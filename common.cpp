
#include <llvm_transcoder.h>
//
#include "common.h"

Common *c = nullptr;  // init in main window creator or audio will not work on
                      // static and general deploy scenarios

void Common::generate_llvm_mt() {
  if (transcoder) delete transcoder;
  transcoder = new llvm_transcoder;

  buff->zero();

  if (compiler.is_ok()) {
    int cnt = 0, n_threads = Thread::getnthreads();

    // scaler list for mt support
    QVector<Scaler<double>> scalers(n_threads);
    int n_values = int(n_channels * secs * sample_rate);
    QVector<float> snd(n_values);

    // calc scaler
    compiler.execute_const();
    transcoder->transcode_mt(compiler, n_threads);

    if (transcoder->ok()) {
      timer.start();

      Thread(iget_n_samples()).run([&](int th, int from, int to) {
        double t = t_inc * from;  // thread timing
        for (int samp = from, isamp = from * n_channels;
             samp < to && control_thread.is_running();
             samp++, t = (t > pi2 * secs) ? 0 : t + t_inc) {
          for (int ch = 0; ch < n_channels; ch++, isamp++) {
            double y = transcoder->evaluate(t, ch, th);
            scalers[th].update(y);
            snd[isamp] = float(y);

            if (++cnt % 50000 == 0) {
              set_progress(100. * cnt / n_values);
              emit progress(
                  QString("%1 %").arg(100. * cnt / n_values, 3, 'f', 0));
            }
          }
        }
      });

      if (control_thread.is_running()) {  // not aborted
        scaler.init(scalers);             // update scaler w/scalers
        scaled_amp = scaler.scale_diff(volume * max_amp());
        float fsa = float(scaled_amp);

        // fill audio buffer w/scaled values buff[i]=scaled_amp * snd[i]
        Thread(n_values).run([&](int i) { buff->set_value(i, fsa * snd[i]); });

        data_completed = true;
        lap = timer.lap();

        emit progress("100 %");
        set_progress(100);
        emit finished();
      }
    }
  }
}

string Common::get_ir_code() { return transcoder->ir_code; }

void Common::generate_mono() {
  buff->zero();

  if (compiler.is_ok()) {
    int cnt = 0;

    // compiler list for mt support

    int n_values = int(n_channels * secs * sample_rate);
    QVector<float> snd(n_values);

    // calc scaler
    compiler.execute(0, 0);

    timer.start();

    for (int samp = 0, isamp = 0; samp < iget_n_samples();
         samp++, t = (t > pi2 * secs) ? 0 : t + t_inc) {
      double t = t_inc * samp;
      for (int ch = 0; ch < n_channels; ch++, isamp++) {
        double y = compiler.execute(t, ch);
        scaler.update(y);
        snd[isamp] = float(y);

        if (++cnt % 50000 == 0) {
          set_progress(100. * cnt / n_values);
          emit progress(QString("%1 %").arg(100. * cnt / n_values, 3, 'f', 0));
        }
      }
    }

    scaled_amp = scaler.scale_diff(volume * max_amp());
    float fsa = float(scaled_amp);

    // fill audio buffer w/scaled values buff[i]=scaled_amp * snd[i]
    Thread(n_values).run([&](int i) { buff->set_value(i, fsa * snd[i]); });

    data_completed = true;
    lap = timer.lap();

    emit progress("100 %");
    set_progress(100);
    emit finished();
  }
}
