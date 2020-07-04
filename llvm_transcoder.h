#pragma once

#include <MusicFreq.h>
#include <VSL_compiler.h>
#include <llvm.h>
#include <math.h>

#include <sstream>

using namespace std;

extern "C" {

double wave(double amp, double hz, double phase) {
  return amp * sin(hz + phase);
}

double wave2(double amp, double hz) { return wave(amp, hz, 0); }

double wave1(double hz) { return wave(1, hz, 0); }

double note(double note_n, double oct)  // note(note#, octave)
{
  if (note_n > 12 || note_n < 0) note_n = 0;
  if (oct < 0 || oct > 10) oct = 0;
  //          s2 = sin(x * NoteOct2Freq(int(s2), int(s1)));
  return MusicFreq::NoteOct2Freq(int(note_n), int(oct));
}

double tone(double note, double oct)  // tone(note#, octave)
{
  if (note > 12 || note < 0) note = 0;
  if (oct < 0 || oct > 10) oct = 0;
  return MusicFreq::NoteOct2Freq(int(note), int(oct));
}

double note_const(double t, double note, double oct) {
  return sin(MusicFreq::NoteOct2Freq(note, oct) * t);
}

double sec(double x) { return 2 * M_PI * x; }

double osc(double hz) { return wave1(hz); }

double lap(double t, double s1, double s2) {
  return (s2 <= s1) ? 0 : (t >= s1 * (2 * M_PI)) && (t <= s2 * (2 * M_PI));
}

double hz2oct(double freq, double oct)  // hz2oct(hz,oct)
{
  return MusicFreq::FreqInOctave(int(freq), int(oct));
}

double factorial(double f) {
  if (f <= 0)
    return 1;
  else
    return f * factorial(f - 1);
}

double saw(double x, double alpha1) {
  static double tn, x1, a1;
  double ret, fm = x;

  if (alpha1 != a1) {
    tn = tan(alpha1 * M_PI / 180.);
    x1 = M_PI / (1 + tn * tn);
  }

  if (x > TWO_PI) fm = fmod(x, TWO_PI);
  if (fm <= x1)
    ret = fm * tn;
  else if (fm <= TWO_PI - x1)
    ret = (-1 / tn) * (fm - M_PI);
  else
    ret = tn * (fm - TWO_PI);
  return ret;
}

// saw(x)
// evaluate between 0..2PI, return -1..+1
double saw1(double x) {
  double ret, fm = x;
  if (x > TWO_PI) fm = fmod(x, TWO_PI);
  if (fm <= M_PI_2)
    ret = fm / M_PI_2;
  else if (fm <= M_PI3_2)
    ret = -2 * (fm / M_PI - 1);
  else
    ret = fm / M_PI_2 - 4;
  return ret;
}

double yinyang(double x, double t) {
  auto k = 6. * M_PI;
  return sin(t * x) * sin(x / (t + k));
}

double freq_mesh(double base, double slope, double islope, double n, double t,
                 FreqMesh *fm) {
  if (!fm->is_init()) fm->set_param(base, slope, uint(islope), uint(n));

  return fm->gen_sample(t);
}

double sequence(double ini, double end, double n, void *_stack_pnt) {
  auto _stack = (VSL_compiler::Stack<double> *)_stack_pnt;

  if (n < _stack->max_stack - 10 && end != ini) {
    if (end < ini) std::swap(ini, end);
    auto inc = (end - ini) / (n - 1);
    for (auto i = ini; i < end; i += inc) _stack->push(i);
    _stack->push(end);
  }
  return 0;
}

double backslash(int oper, double t, void *_stack_pnt) {
  auto _stack = (VSL_compiler::Stack<double> *)_stack_pnt;
  auto &sp = _stack->sp;

  double res = 0;
  if (sp > 1) {
    if (oper == TILDE) {
      for (sp--; sp >= 0; sp--) res += sin(t * (*_stack)[sp].d);
    } else {
      res = (*_stack)[-1].d;
      for (sp -= 2; sp >= 0; sp--) {
        switch (oper) {
          case PLUS:
            res += (*_stack)[sp].d;
            break;
          case MINUS:
            res -= (*_stack)[sp].d;
            break;
          case MULT:
            res *= (*_stack)[sp].d;
            break;
          case DIV:
            res /= (*_stack)[sp].d;
            break;
        }
      }
    }
    (*_stack)[0] = res;
    sp = 1;
  }
  return res;
}

//  "magneticring"
}

class JIT {
 public:
  LLVMContext Context;
  IRBuilder<> Builder = IRBuilder<>(Context);
  unique_ptr<Module> module;
  ExecutionEngine *engine = nullptr;
  EngineBuilder *engine_builder = nullptr;
  string error_msg, ir_code;

 public:
  JIT(string module_name = "test") {
    InitializeNativeTarget();
    InitializeNativeTargetAsmPrinter();
    InitializeNativeTargetAsmParser();

    module = make_unique<Module>(module_name, Context);
  }
  ~JIT() {
    if (engine) delete engine;
    if (engine_builder) delete engine_builder;
  }

  bool create_engine() {  //  create execution engine
    engine_builder = new EngineBuilder(move(module));
    engine_builder->setErrorStr(&error_msg).setEngineKind(EngineKind::JIT);

    // puts(error_msg.c_str());
    engine = engine_builder->create();

    return engine != nullptr;
  }

  void dump_code() {  // Output the bitcode file to stdout
    WriteBitcodeToFile(*module.get(), outs());
  }

  void show_code() {
    module->print(outs(), nullptr);  // print IR
  }

  string get_ir_code() {
    string str;

    raw_string_ostream os(str);

    os << *module;
    os.flush();

    return str;
  }

  Module *get_module() { return module.get(); }

  intptr_t get_function(string name) {
    return engine->getFunctionAddress(name);
  }

  vector<Type *> arg_types(vector<Value *> args) {
    vector<Type *> arg_tys;
    for (Value *v : args) arg_tys.push_back(v->getType());
    return arg_tys;
  }

  vector<Type *> arg_int(int n) {
    vector<Type *> vt;
    for (int i = 0; i < n; i++) vt.push_back(int32_type());
    return vt;
  }

  vector<Type *> arg_fp(int n) {
    vector<Type *> vt;
    for (int i = 0; i < n; i++) vt.push_back(fp_type());
    return vt;
  }

  vector<Type *> arg_dbl(int n) {
    vector<Type *> vt;
    for (int i = 0; i < n; i++) vt.push_back(dbl_type());
    return vt;
  }

  vector<Type *> arg_ptr(int n) {
    vector<Type *> vt;
    for (int i = 0; i < n; i++) vt.push_back(ptr_type());
    return vt;
  }

  int64_t get_int(Value *v) { return dyn_cast<ConstantInt>(v)->getSExtValue(); }
  double get_dbl(Value *v) {
    return dyn_cast<ConstantFP>(v)->getValueAPF().convertToDouble();
  }
  double get_flt(Value *v) {
    return dyn_cast<ConstantFP>(v)->getValueAPF().convertToFloat();
  }

  Type *int32_type() { return Type::getInt32Ty(Context); }
  Type *int_type() { return Type::getInt32Ty(Context); }
  Type *fp_type() { return Type::getFloatTy(Context); }
  Type *dbl_type() { return Type::getDoubleTy(Context); }
  Type *void_type() { return Type::getVoidTy(Context); }
  Type *ptr_type() { return Type::getInt64Ty(Context); }

  Value *init_int(string name = "", int i = 0) {
    auto v = ConstantInt::get(Type::getInt32Ty(Context), i);
    v->setName(name);
    return v;
  }
  Value *init_fp(string name = "", float f = 0) {
    auto v = ConstantFP::get(Type::getFloatTy(Context), f);
    v->setName(name);
    return v;
  }
  Value *init_dbl(string name = "", double d = 0) {
    auto v = ConstantFP::get(Type::getDoubleTy(Context), d);
    v->setName(name);
    return v;
  }

  Value *init_ptr(string name = "", void *ptr = nullptr) {
    auto *addr = ConstantInt::get(Type::getInt64Ty(Context), u_int64_t(ptr));
    Value *v = ConstantExpr::getIntToPtr(
        addr, PointerType::getUnqual(Type::getInt64Ty(Context)));

    v->setName(name);
    return v;
  }

  GlobalVariable *gl_dbl(string id, double d, bool is_constant) {
    module->getOrInsertGlobal(id, Builder.getDoubleTy());
    GlobalVariable *gv = module->getNamedGlobal(id);
    gv->setLinkage(GlobalValue::ExternalLinkage);
    gv->setInitializer(ConstantFP::get(Type::getDoubleTy(Context), d));
    gv->setConstant(is_constant);
    return gv;
  }

  GlobalVariable *gl_array_dbl(string id, int dim) {
    module->getOrInsertGlobal(id, ArrayType::get(dbl_type(), dim));
    GlobalVariable *gv = module->getNamedGlobal(id);
    gv->setLinkage(GlobalValue::ExternalLinkage);

    ArrayRef<Constant *> vd;

    auto _vec = vd.vec();
    for (int i = 0; i < dim; i++) _vec.push_back(0);

    gv->setInitializer(ConstantArray::get(ArrayType::get(dbl_type(), dim), vd));
    return gv;
  }

  Value *index_pointer(GlobalVariable *array, Value *index) {
    return Builder.CreateGEP(
        array, {ConstantInt::get(Context, APInt(64, 0, true)), index});
  }
};

// reads generated code and transcodes to llvm ir
class llvm_transcoder : public JIT {
  using func_type = double(double);  // y=chan_func(x), single thread
  vector<func_type *> func;

  using func_type_mt = double(double, int);  // y=chan_func(x, nth), multithread
  vector<func_type_mt *> func_mt;

  map<int, Function *> funcs;

  VSL_compiler::Stack<double> *_stack = nullptr;

  class Stack {
   public:
    vector<Value *> stack;

    Value *get(int sp_offset = 0) {  // get item _stack[-sp_offset]
      return sp_offset == 0 ? stack.back()
                            : *::prev(stack.end(), sp_offset + 1);
    }

    Value *top() { return stack.back(); }
    Value *prev() { return get(1); }  // stack[sp-1]

    void push(Value *v) { stack.push_back(v); }
    void pop(int n_times = 1) {
      for (int i = 0; i < n_times; i++) stack.pop_back();
    }

    Value *bin_oper(Instruction *instr) {
      pop(2);
      push(instr);
      return instr;
    }

    Value *bin_oper(Value *v) {
      pop(2);
      push(v);
      return v;
    }

    Value *tri_oper(Value *instr) {
      pop(3);
      push(instr);
      return instr;
    }

    Value *unit_oper(Value *instr) {
      pop();
      push(instr);
      return instr;
    }

    vector<Value *> get_n_items(int n) {
      vector<Value *> vv;
      if (n <= size())
        for (int i = 0; i < n; i++) {
          vv.push_back(top());
          pop();
        }
      std::reverse(vv.begin(), vv.end());
      return vv;
    }
    int size() { return stack.size(); }

    void clear() { stack.clear(); }
  };

  Stack stack;

  char *code;
  bool err = false;
  int pc = 0;
  VSL_compiler::Compiler *compiler = nullptr;

  map<string, GlobalVariable *> numerical_values;
  int n_threads = 0;  // single thread
  map<string, Function *> local_funcs;

  int get_iparam(int n = 0) {
    return *reinterpret_cast<int *>(code + pc + 1 + n * sizeof(int));
  }
  double get_dparam(int n = 0) {
    return *reinterpret_cast<double *>(code + pc + 1 + n * sizeof(int));
  }

 public:
  llvm_transcoder() : JIT("expr_main") {}

  int generate(Value *t_value, Value *nth_value) {
    int sym = code[pc];
    double _d;

    switch (sym) {
      case PUSH_CONST:
        _d = get_dparam();
        pc += sizeof(double) + 1;

        stack.push(init_dbl("constant", _d));
        break;
      case PUSH_T:
        pc++;
        stack.push(t_value);
        break;
      case PUSH_ID: {
        int nid = get_iparam();
        string id = compiler->get_ids(nid);
        auto tp = compiler->get_type(nid);
        pc += sizeof(int) + 1;
        if (n_threads == 0)  // single thread
          stack.push(Builder.CreateLoad(
              dbl_type(), numerical_values[compiler->get_ids(nid)]));
        else {  // mt-> index numerical_values[id][nth_value]
          if (numerical_values[compiler->get_ids(nid)]->isConstant())
            stack.push(Builder.CreateLoad(
                dbl_type(), numerical_values[compiler->get_ids(nid)]));
          else
            stack.push(Builder.CreateLoad(index_pointer(
                numerical_values[compiler->get_ids(nid)], nth_value)));
        }

      } break;

      case RET:  // ret int(n_params)
        pc += sizeof(int) + 1;
        Builder.CreateRet(stack.top());
        break;

      case FUNC: {  // last argument is 't, nth'
        int addr = get_iparam(), n_params = get_iparam(1);

        auto tf = compiler->get_by_address(addr);
        string id = tf.sid;

        auto params = stack.get_n_items(n_params);
        params.push_back(t_value);
        params.push_back(nth_value);

        stack.push(Builder.CreateCall(local_funcs[id], params));

        pc += sizeof(int) * 2 + 1;
      } break;

      case PARAM: {
        int n_param = get_iparam(), i_func = get_iparam(1);

        string id = compiler->get_ids(i_func);
        stack.push(local_funcs[id]->arg_begin() + n_param);  // 't' is last
        pc += 2 * sizeof(int) + 1;
      } break;

      case POP: {
        int nid = get_iparam();

        if (n_threads == 0) {  // single thread
          Builder.CreateStore(stack.top(),
                              numerical_values[compiler->get_ids(nid)]);
        } else {  // mt -> index numerical_values[id][nth_value]
          if (numerical_values[compiler->get_ids(nid)]->isConstant())
            Builder.CreateStore(stack.top(),
                                numerical_values[compiler->get_ids(nid)]);
          else
            Builder.CreateStore(
                stack.top(),
                index_pointer(numerical_values[compiler->get_ids(nid)],
                              nth_value));
        }
        stack.pop();
        pc += sizeof(int) + 1;
      }

      break;

      case NEG:
        pc++;
        if (stack.size() >= 1)
          stack.unit_oper(Builder.CreateFNeg(stack.top()));
        else
          err = true;
        break;

      case PLUS:
        pc++;
        if (stack.size() >= 2)
          stack.bin_oper(Builder.CreateFAdd(stack.top(), stack.prev()));
        else
          err = true;
        break;

      case MINUS:
        pc++;
        if (stack.size() >= 2)
          stack.bin_oper(Builder.CreateFSub(stack.prev(), stack.top()));
        else
          err = true;
        break;

      case MULT:
        pc++;
        if (stack.size() >= 2)
          stack.bin_oper(Builder.CreateFMul(stack.top(), stack.prev()));
        else
          err = true;
        break;
      case DIV:
        pc++;
        if (stack.size() >= 2)
          stack.bin_oper(Builder.CreateFDiv(stack.prev(), stack.top()));
        else
          err = true;
        break;

      case POWER:
        if (stack.size() >= 2) {
          stack.bin_oper(Builder.CreateCall(create_func(sym),
                                            {stack.prev(), stack.top()}));
        } else
          err = true;
        break;

      case EQ:
      case NE:
      case GT:
      case LT:
      case LE:
      case GE: {
        static map<int, CmpInst::Predicate> op_map{
            {EQ, CmpInst::FCMP_OEQ}, {NE, CmpInst::FCMP_ONE},
            {GT, CmpInst::FCMP_OGT}, {LT, CmpInst::FCMP_OLT},
            {LE, CmpInst::FCMP_OLE}, {GE, CmpInst::FCMP_OGE}};

        if (stack.size() >= 2) {
          stack.bin_oper(Builder.CreateUIToFP(
              Builder.CreateFCmp(op_map[sym], stack.prev(), stack.top()),
              fp_type()));

        } else
          err = true;
      } break;

      case SWAVE:  // wave(amp(2), hz(1), phase(top))
        pc++;
        if (stack.size() >= 3) {                                   // 3 arg func
          Value *hzt = Builder.CreateFMul(t_value, stack.get(1));  // hzt=hz*t
          stack.tri_oper(Builder.CreateCall(create_func(sym),
                                            {stack.get(2), hzt, stack.top()}));

        } else
          err = true;
        break;

      case SWAVE2:  // wave(amp, hz)
        pc++;
        if (stack.size() >= 2) {                                  // 2 arg func
          Value *hzt = Builder.CreateFMul(t_value, stack.top());  // hzt=hz*t
          stack.bin_oper(
              Builder.CreateCall(create_func(sym), {stack.prev(), hzt}));
        } else
          err = true;
        break;

      case OSC:
      case SWAVE1:  // wave(hz)
        pc++;
        if (stack.size() >= 1) {                                  // 1 arg func
          Value *hzt = Builder.CreateFMul(t_value, stack.top());  // hzt=hz*t
          stack.unit_oper(Builder.CreateCall(create_func(sym), {hzt}));
        } else
          err = true;
        break;

      case YINYANG:  // yy(x)
        pc++;
        stack.unit_oper(
            Builder.CreateCall(create_func(sym), {stack.top(), t_value}));
        break;

      case TONE:  // 2 parameter funcs
      case NOTE:
      case LAP:
      case HZ2OCT:
        pc++;
        if (stack.size() >= 2) {  // 2 arg func
          stack.bin_oper(Builder.CreateCall(create_func(sym),
                                            {stack.prev(), stack.top()}));
        } else
          err = true;
        break;

      case NOTE_CONST:  // note_const(t, note, oct)
      {
        Value *note = init_dbl("note", get_iparam(0)),
              *oct = init_dbl("oct", get_iparam(1));
        stack.push(Builder.CreateCall(create_func(sym), {t_value, note, oct}));
      }
        pc += 2 * sizeof(int) + 1;
        break;

      case FREQ_MESH: {  // *(base, slope, islope, n)
        int n_mesh = get_iparam();
        pc += sizeof(int) + 1;

        if (stack.size() >= 4) {
          FreqMesh *fm = &compiler->freq_mesh[n_mesh];
          auto params = stack.get_n_items(4);
          params.push_back(t_value);
          params.push_back(init_ptr("mesh", fm));

          stack.push(Builder.CreateCall(create_func(sym), params));

        } else {
          err = true;
        }
      } break;

      case SEQUENCE: {  // (ini, end, n) works on compiler _stack
        pc++;
        if (stack.size() >= 3) {
          auto params = stack.get_n_items(3);
          params.push_back(init_ptr("_stack", _stack));
          Builder.CreateCall(create_func(sym), params);
        }
      } break;

      case BACKSLASH: {  // \{}operator on compiler _stack: (int oper, double t,
                         // void *_stack_pnt)
        pc++;
        // has compiler _stack? -> sequence
        ArrayRef<Value *> params{init_int("oper", code[pc]), t_value,
                                 init_ptr("_stack", _stack)};
        stack.push(Builder.CreateCall(create_func(sym), params));

        if (stack.size() > 1) {
          int oper = code[pc];

          if (oper == TILDE) {
            while (stack.size() > 1)
              stack.bin_oper(
                  Builder.CreateCall(create_func(SWAVE1), stack.top()));
          } else {
            auto res = stack.top();
            while (stack.size() > 1) {
              switch (oper) {
                case PLUS:
                  res = stack.bin_oper(Builder.CreateFAdd(res, stack.prev()));
                  break;
                case MINUS:
                  res = stack.bin_oper(Builder.CreateFSub(res, stack.prev()));
                  break;
                case MULT:
                  res = stack.bin_oper(Builder.CreateFMul(res, stack.prev()));
                  break;
                case DIV:
                  res = stack.bin_oper(Builder.CreateFDiv(res, stack.prev()));
                  break;
              }
            }
          }
        }

        pc++;
      } break;

      case FACT:
      case FSIN:
      case FCOS:
      case FTAN:
      case FEXP:
      case FLOG:
      case FLOG10:
      case FINT:
      case FSQRT:
      case FASIN:
      case FACOS:
      case FATAN:
      case FABS:

      case SEC:

      case ABS:
      case SAW:
      case SAW1:

      case MAGNETICRING:
        pc++;
        if (stack.size() >= 1) {  // 1 arg func
          stack.unit_oper(Builder.CreateCall(create_func(sym), {stack.top()}));
        } else
          err = true;

        break;

      default:
        err = true;
        break;
    }

    return sym;
  }

  // multi threaded version
  void transcode_mt(VSL_compiler &vsl_compiler, int n_threads) {
    this->n_threads = n_threads;

    err = false;
    compiler = vsl_compiler.get_compiler();
    _stack = &vsl_compiler._stack;
    code = compiler->get_code();

    // global const / let vars -> get values from Table_Values
    numerical_values.clear();
    for (auto &tv : compiler->get_numericals()) {
      if (tv.is_const) {
        numerical_values[tv.sid] =
            gl_dbl(tv.sid, tv.get_double(), true);  // constant
      } else
        numerical_values[tv.sid] = gl_array_dbl(
            tv.sid, n_threads);  // let: global double array one item per thread
    }

    // let's, if any
    Function *let_func = nullptr;

    stack.clear();
    if (vsl_compiler.has_lets()) {
      let_func =  // let_func(double t, int n_thread)
          Function::Create(
              FunctionType::get(void_type(), {dbl_type(), int_type()}, false),
              Function::InternalLinkage, "let_func", get_module());
      Builder.SetInsertPoint(
          BasicBlock::Create(Context, "entry_let", let_func));

      auto blk_let = vsl_compiler.get_block_let();

      Value *t = let_func->arg_begin(), *nth = let_func->arg_begin() + 1;
      for (pc = blk_let.first; pc < blk_let.second && !err;) {
        generate(t, nth);  // t, nth
      }
      Builder.CreateRetVoid();
    }

    // func's
    stack.clear();
    local_funcs.clear();

    if (vsl_compiler.has_funcs()) {
      for (auto &tf : vsl_compiler.get_funcs()) {
        // pass 'double t, int nth_value' as last argument (+2)
        vector<Type *> params = arg_dbl(tf.n_params + 1);
        params.push_back(int_type());

        local_funcs[tf.sid] =
            Function::Create(FunctionType::get(dbl_type(), params, false),
                             Function::ExternalLinkage, tf.sid, get_module());

        // new_fp_func(tf.sid, tf.n_params + 2);
        Builder.SetInsertPoint(BasicBlock::Create(
            Context, "entry_local_func_" + tf.sid, local_funcs[tf.sid]));

        Value *t = local_funcs[tf.sid]->arg_begin() + tf.n_params,
              *nth = local_funcs[tf.sid]->arg_begin() + tf.n_params + 1;
        while (generate(t, nth) != RET && !err &&
               pc < compiler->get_code_size())
          ;
      }
    }

    // expression code
    stack.clear();
    for (int chan = 0; chan < compiler->get_chans() && !err; chan++) {
      // main callable double expr_result = double chan_func#ch(double x, int
      // nth)
      Function *chan_function = Function::Create(
          FunctionType::get(dbl_type(), {dbl_type(), int_type()}, false),
          Function::ExternalLinkage, "chan_func" + itostr(chan), get_module());
      Builder.SetInsertPoint(BasicBlock::Create(
          Context, "entry_chan" + itostr(chan), chan_function));

      Value *t = chan_function->arg_begin(),
            *nth = chan_function->arg_begin() + 1;

      if (chan == 0 && let_func)  // call let_func(t, nth) only once on chan 0
        Builder.CreateCall(let_func, {t /* t */, nth} /*  nth */);

      for (pc = vsl_compiler.get_block_code(chan).first;
           pc < vsl_compiler.get_block_code(chan).second && !err;) {
        generate(t, nth);
      }

      if (stack.size() == 1 && !err) {
        Builder.CreateRet(stack.top());

      } else
        err = true;

      stack.clear();
    }

    if (!err) {
      ir_code = get_ir_code();
      //      puts(ir_code.c_str());

      if (!create_engine()) {
        err = true;
      } else {  // set the 'c' callable func that evaluates the rpn expression
        func.clear();
        func_mt.clear();

        for (int chan = 0; chan < compiler->get_chans(); chan++) {
          func.push_back(reinterpret_cast<func_type *>(
              ((void *)get_function("chan_func" + itostr(chan)))));

          func_mt.push_back(reinterpret_cast<func_type_mt *>(
              ((void *)get_function("chan_func" + itostr(chan)))));
        }
      }
    }
  }

  // double name(double di(n_params))
  Function *dbl_func(string name, int n_params) {
    Function *fnc = Function::Create(
        FunctionType::get(dbl_type(), arg_dbl(n_params), false),
        Function::ExternalLinkage, name, get_module());
    fnc->setCallingConv(CallingConv::C);
    return fnc;
  }

  Function *gen_func(string name, vector<Type *> params) {
    Function *fnc =
        Function::Create(FunctionType::get(dbl_type(), params, false),
                         Function::ExternalLinkage, name, get_module());
    fnc->setCallingConv(CallingConv::C);
    return fnc;
  }

  Function *dbl_func_ptr(string name, int n_params) {
    auto params = arg_dbl(n_params);
    params.push_back(ptr_type());

    Function *fnc =
        Function::Create(FunctionType::get(dbl_type(), params, false),
                         Function::ExternalLinkage, name, get_module());
    fnc->setCallingConv(CallingConv::C);
    return fnc;
  }

  Function *create_func(int sym) {  // create internal func on demand by 'sym'
    static const vector<string> fname = {
        "sin",   "cos",  "tan",    "exp",         "log",  "log10", "floor",
        "sqrt",  "asin", "acos",   "atan",        "fabs", "wave",  "wave1",
        "wave2", "tone", "note",   "sec",         "osc",  "fabs",  "saw",
        "saw1",  "lap",  "hz2oct", "magneticring"};

    if (funcs.find(sym) == funcs.end()) switch (sym) {
        case POWER:
          funcs[sym] = dbl_func("pow", 2);
          break;
        case SWAVE:
          funcs[sym] = dbl_func("wave", 3);
          break;
        case SWAVE2:
          funcs[sym] = dbl_func("wave2", 2);
          break;
        case FACT:
          funcs[sym] = dbl_func("factorial", 1);
          break;
        case YINYANG:
          funcs[sym] = dbl_func("yinyang", 2);
          break;
        case NOTE_CONST:
          funcs[sym] = dbl_func("note_const", 3);
          break;
        case SEQUENCE:
          funcs[sym] = gen_func(
              "sequence", {dbl_type(), dbl_type(), dbl_type(), ptr_type()});
          break;
        case BACKSLASH:  // backslash(int oper, double t, void *_stack_pnt)
          funcs[sym] =
              gen_func("backslash", {int_type(), dbl_type(), ptr_type()});
          break;
        case FREQ_MESH:
          funcs[sym] =
              gen_func("freq_mesh", {dbl_type(), dbl_type(), dbl_type(),
                                     dbl_type(), dbl_type(), ptr_type()});
          break;
        case HZ2OCT:
          funcs[sym] = dbl_func("hz2oct", 2);
          break;

        default:  // 1 arg funcs
          funcs[sym] = dbl_func(fname[sym - FSIN], 1);
          break;
      }

    return funcs[sym];
  }

  // single thread
  inline double evaluate(double x, int ch) { return func[ch](x); }
  // multi thread (mt)
  inline double evaluate(double x, int ch, int nth) {
    return func_mt[ch](x, nth);
  }

  bool ok() { return !err; }
};
