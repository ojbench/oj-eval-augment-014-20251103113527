#include "Evalvisitor.h"
#include <cmath>
#include <iomanip>
#include <sstream>

using std::any;
using std::any_cast;
using std::get;
using std::get_if;
using std::holds_alternative;
using std::make_pair;
using std::optional;
using std::pair;
using std::string;
using std::vector;

// ===== Utilities =====
std::string EvalVisitor::strip_quotes(const std::string& s) {
    if (s.size() >= 2 && ((s.front()=='"' && s.back()=='"') || (s.front()=='\'' && s.back()=='\''))) {
        return s.substr(1, s.size()-2);
    }
    return s;
}

// ===== Env helpers =====
std::unordered_map<std::string, Value>& EvalVisitor::current_locals() {
    if (local_stack_.empty()) {
        static std::unordered_map<std::string, Value> dummy; // shouldn't happen
        return dummy;
    }
    return local_stack_.back();
}

Value EvalVisitor::getVar(const std::string& name) const {
    if (!local_stack_.empty()) {
        const auto& loc = local_stack_.back();
        auto it = loc.find(name);
        if (it != loc.end()) return it->second;
    }
    auto itg = globals_.find(name);
    if (itg != globals_.end()) return itg->second;
    return std::monostate{}; // None by default
}

void EvalVisitor::setVar(const std::string& name, const Value& v) {
    if (!local_stack_.empty()) {
        auto& loc = local_stack_.back();
        auto it = loc.find(name);
        if (it != loc.end()) { it->second = v; return; }
    }
    globals_[name] = v;
}

// ===== Conversions =====
bool EvalVisitor::as_bool_strict(const Value& v) {
    if (holds_alternative<bool>(v)) return get<bool>(v);
    throw std::runtime_error("not a bool");
}

bool EvalVisitor::to_bool(const Value& v) {
    if (holds_alternative<std::monostate>(v)) return false;
    if (holds_alternative<bool>(v)) return get<bool>(v);
    if (holds_alternative<long long>(v)) return get<long long>(v)!=0;
    if (holds_alternative<double>(v)) return get<double>(v)!=0.0;
    if (holds_alternative<std::string>(v)) return !get<std::string>(v).empty();
    return false;
}

long long EvalVisitor::to_int(const Value& v) {
    if (holds_alternative<long long>(v)) return get<long long>(v);
    if (holds_alternative<bool>(v)) return get<bool>(v) ? 1 : 0;
    if (holds_alternative<double>(v)) return static_cast<long long>(get<double>(v)); // truncate toward 0
    if (holds_alternative<std::string>(v)) {
        try {
            return std::stoll(get<std::string>(v));
        } catch (...) { return 0; }
    }
    return 0;
}

double EvalVisitor::to_float(const Value& v) {
    if (holds_alternative<double>(v)) return get<double>(v);
    if (holds_alternative<long long>(v)) return static_cast<double>(get<long long>(v));
    if (holds_alternative<bool>(v)) return get<bool>(v) ? 1.0 : 0.0;
    if (holds_alternative<std::string>(v)) {
        try {
            return std::stod(get<std::string>(v));
        } catch (...) { return 0.0; }
    }
    return 0.0;
}

std::string EvalVisitor::to_string(const Value& v) {
    if (holds_alternative<std::monostate>(v)) return "None";
    if (holds_alternative<bool>(v)) return get<bool>(v) ? "True" : "False";
    if (holds_alternative<long long>(v)) return std::to_string(get<long long>(v));
    if (holds_alternative<double>(v)) {
        std::ostringstream oss; oss.setf(std::ios::fixed); oss<<std::setprecision(6)<<get<double>(v); return oss.str();
    }
    if (holds_alternative<std::string>(v)) return get<std::string>(v);
    return "";
}

// ===== Binary ops =====
static bool is_int(const Value& v){return holds_alternative<long long>(v);}
static bool is_float(const Value& v){return holds_alternative<double>(v);}
static bool is_num(const Value& v){return is_int(v)||is_float(v);}
static bool is_str(const Value& v){return holds_alternative<std::string>(v);}

Value EvalVisitor::add(const Value& a, const Value& b) {
    if (is_str(a) && is_str(b)) return get<string>(a) + get<string>(b);
    if (is_num(a) && is_num(b)) {
        if (is_float(a) || is_float(b)) return to_float(a) + to_float(b);
        return get<long long>(a) + get<long long>(b);
    }
    // bool participates like int
    if (holds_alternative<bool>(a) || holds_alternative<bool>(b)) {
        if (is_str(a) || is_str(b)) throw std::runtime_error("bad +");
        double fa = to_float(a), fb = to_float(b);
        if (!std::floor(fa) && !std::floor(fb)) {}
        if ((is_float(a) || is_float(b))) return fa+fb;
        return static_cast<long long>(to_int(a) + to_int(b));
    }
    throw std::runtime_error("unsupported +");
}

Value EvalVisitor::sub(const Value& a, const Value& b) {
    if (is_num(a) && is_num(b)) {
        if (is_float(a) || is_float(b)) return to_float(a) - to_float(b);
        return get<long long>(a) - get<long long>(b);
    }
    return static_cast<long long>(to_int(a) - to_int(b));
}

Value EvalVisitor::mul(const Value& a, const Value& b) {
    if (is_str(a) && is_num(b)) {
        long long n = to_int(b); if (n<=0) return std::string("");
        std::string res; res.reserve(get<string>(a).size()*std::min<long long>(n, 100000));
        for (long long i=0;i<n;++i) res += get<string>(a);
        return res;
    }
    if (is_str(b) && is_num(a)) return mul(b,a);
    if (is_num(a) && is_num(b)) {
        if (is_float(a) || is_float(b)) return to_float(a) * to_float(b);
        return get<long long>(a) * get<long long>(b);
    }
    return static_cast<long long>(to_int(a) * to_int(b));
}

Value EvalVisitor::truediv(const Value& a, const Value& b) {
    return static_cast<double>(to_float(a) / to_float(b));
}

static long long floordiv_ll(long long a, long long b){
    if (b==0) return 0; // undefined
    long long q = a / b; long long r = a % b;
    if ((r!=0) && ((a<0) != (b<0))) --q; // floor toward -inf
    return q;
}

Value EvalVisitor::floordiv(const Value& a, const Value& b) {
    if (is_num(a) && is_num(b) && !is_float(a) && !is_float(b)) {
        return floordiv_ll(get<long long>(a), get<long long>(b));
    }
    // mixed -> use double then floor
    double q = std::floor(to_float(a) / to_float(b));
    return static_cast<long long>(q);
}

Value EvalVisitor::mod(const Value& a, const Value& b) {
    long long ia = to_int(a), ib = to_int(b);
    long long q = floordiv_ll(ia, ib);
    return ia - q * ib;
}

Value EvalVisitor::unary_plus(const Value& a) {
    if (is_num(a)) return a;
    return static_cast<long long>(to_int(a));
}

Value EvalVisitor::unary_minus(const Value& a) {
    if (is_num(a)) {
        if (is_float(a)) return -get<double>(a);
        return -get<long long>(a);
    }
    return static_cast<long long>(-to_int(a));
}

int EvalVisitor::cmp(const Value& a, const Value& b) {
    if (is_str(a) && is_str(b)) {
        const auto& sa = get<string>(a); const auto& sb = get<string>(b);
        if (sa<sb) return -1; if (sa>sb) return 1; return 0;
    }
    // try numeric
    if (is_num(a) || is_num(b) || holds_alternative<bool>(a) || holds_alternative<bool>(b)) {
        if (is_float(a) || is_float(b)) {
            double da = to_float(a), db = to_float(b);
            if (da<db) return -1; if (da>db) return 1; return 0;
        } else {
            long long ia = to_int(a), ib = to_int(b);
            if (ia<ib) return -1; if (ia>ib) return 1; return 0;
        }
    }
    // None comparisons -> only == and !=, others false
    if (is_none(a) && is_none(b)) return 0;
    return -2; // signal incomparable for ordering
}

// ===== Statement visitors =====
std::any EvalVisitor::visitFile_input(Python3Parser::File_inputContext* ctx) {
    for (auto st : ctx->stmt()) {
        visit(st);
    }
    return {};
}

std::any EvalVisitor::visitStmt(Python3Parser::StmtContext* ctx) {
    if (ctx->simple_stmt()) return visit(ctx->simple_stmt());
    if (ctx->compound_stmt()) return visit(ctx->compound_stmt());
    return {};
}

std::any EvalVisitor::visitSimple_stmt(Python3Parser::Simple_stmtContext* ctx) {
    return visit(ctx->small_stmt());
}

std::any EvalVisitor::visitSmall_stmt(Python3Parser::Small_stmtContext* ctx) {
    if (ctx->expr_stmt()) return visit(ctx->expr_stmt());
    if (ctx->flow_stmt()) return visit(ctx->flow_stmt());
    return {};
}

std::any EvalVisitor::visitSuite(Python3Parser::SuiteContext* ctx) {
    if (ctx->simple_stmt()) {
        visit(ctx->simple_stmt());
        return {};
    }
    // NEWLINE INDENT (stmt)+ DEDENT
    for (auto st : ctx->stmt()) {
        visit(st);
    }
    return {};
}

std::any EvalVisitor::visitCompound_stmt(Python3Parser::Compound_stmtContext* ctx) {
    if (ctx->if_stmt()) return visit(ctx->if_stmt());
    if (ctx->while_stmt()) return visit(ctx->while_stmt());
    if (ctx->funcdef()) return visit(ctx->funcdef());
    return {};
}

std::any EvalVisitor::visitIf_stmt(Python3Parser::If_stmtContext* ctx) {
    int nConds = static_cast<int>(ctx->test().size());
    for (int i=0;i<nConds;++i) {
        if (to_bool(any_cast<Value>(visit(ctx->test(i))))) {
            visit(ctx->suite(i));
            return {};
        }
    }
    // else
    if (ctx->ELSE()) {
        visit(ctx->suite(nConds));
    }
    return {};
}

std::any EvalVisitor::visitWhile_stmt(Python3Parser::While_stmtContext* ctx) {
    while (to_bool(any_cast<Value>(visit(ctx->test())))) {
        try {
            visit(ctx->suite());
        } catch (const ContinueSignal&) {
            continue;
        } catch (const BreakSignal&) {
            break;
        } catch (const ReturnSignal&) {
            throw; // propagate
        }
    }
    return {};
}

std::any EvalVisitor::visitFlow_stmt(Python3Parser::Flow_stmtContext* ctx) {
    if (ctx->break_stmt()) return visit(ctx->break_stmt());
    if (ctx->continue_stmt()) return visit(ctx->continue_stmt());
    if (ctx->return_stmt()) return visit(ctx->return_stmt());
    return {};
}

std::any EvalVisitor::visitBreak_stmt(Python3Parser::Break_stmtContext* /*ctx*/) {
    throw BreakSignal{};
}

std::any EvalVisitor::visitContinue_stmt(Python3Parser::Continue_stmtContext* /*ctx*/) {
    throw ContinueSignal{};
}

std::any EvalVisitor::visitReturn_stmt(Python3Parser::Return_stmtContext* ctx) {
    ReturnSignal rs;
    if (ctx->testlist()) {
        auto anyVals = visit(ctx->testlist());
        if (auto p = any_cast<Values>(&anyVals)) rs.values = *p;
    }
    throw rs;
}

std::any EvalVisitor::visitAugassign(Python3Parser::AugassignContext* /*ctx*/) { return {}; }

// ===== Function def =====
std::any EvalVisitor::visitFuncdef(Python3Parser::FuncdefContext* ctx) {
    std::string name = ctx->NAME()->getText();
    FunctionDef f;
    // parameters
    if (auto pc = ctx->parameters(); pc && pc->typedargslist()) {
        auto tl = pc->typedargslist();
        int n = static_cast<int>(tl->tfpdef().size());
        f.param_names.reserve(n);
        f.default_tests.assign(n, nullptr);
        // Parameters appear interleaved with optional ASSIGN test; we iterate by commas via indices
        int ti = 0; // index into test() defaults
        int ni = 0; // index into names
        for (int i=0;i<n;++i) {
            auto p = tl->tfpdef(i);
            f.param_names.push_back(p->NAME()->getText());
            // If there are fewer ASSIGN than names so far, try to attach defaults from the tail.
            // Simpler: check if (i < tl->ASSIGN().size()) but this is not aligned when some earlier have no default.
            // We use heuristic: defaults are aligned to the last k parameters. We'll attach from the end.
        }
        // Attach defaults from the end
        int dcnt = static_cast<int>(tl->test().size());
        for (int di=0; di<dcnt; ++di) {
            // default belongs to param index: last dcnt params
            int param_index = n - dcnt + di;
            if (param_index >= 0 && param_index < n) f.default_tests[param_index] = tl->test(di);
        }
    }
    f.suite = ctx->suite();
    functions_[name] = f;
    return {};
}

std::any EvalVisitor::visitParameters(Python3Parser::ParametersContext* /*ctx*/) { return {}; }
std::any EvalVisitor::visitTypedargslist(Python3Parser::TypedargslistContext* /*ctx*/) { return {}; }
std::any EvalVisitor::visitTfpdef(Python3Parser::TfpdefContext* /*ctx*/) { return {}; }

// ===== Expr stmt (assign / augassign / expr) =====
std::any EvalVisitor::visitExpr_stmt(Python3Parser::Expr_stmtContext* ctx) {
    int nLists = static_cast<int>(ctx->testlist().size());
    if (ctx->augassign()) {
        // lhs op= rhs, lhs must be single name
        std::string lhsName = ctx->testlist(0)->getText();
        auto rhsValsAny = visit(ctx->testlist(1));
        auto rhsVals = any_cast<Values>(rhsValsAny);
        Value rhs = rhsVals.empty()? Value{} : rhsVals[0];
        Value cur = getVar(lhsName);
        auto aug = ctx->augassign();
        Value nv;
        if (aug->ADD_ASSIGN()) nv = add(cur, rhs);
        else if (aug->SUB_ASSIGN()) nv = sub(cur, rhs);
        else if (aug->MULT_ASSIGN()) nv = mul(cur, rhs);
        else if (aug->DIV_ASSIGN()) nv = truediv(cur, rhs);
        else if (aug->IDIV_ASSIGN()) nv = floordiv(cur, rhs);
        else if (aug->MOD_ASSIGN()) nv = mod(cur, rhs);
        else nv = rhs;
        setVar(lhsName, nv);
        return {};
    }

    if (nLists == 1) {
        // expression statement, evaluate and discard
        visit(ctx->testlist(0));
        return {};
    }
    // chain assignment: a=b=c or tuple
    // Evaluate rightmost once
    auto rhsAny = visit(ctx->testlist(nLists-1));
    auto rhsVals = any_cast<Values>(rhsAny);
    for (int i=0;i<nLists-1;++i) {
        auto tl = ctx->testlist(i);
        int ln = static_cast<int>(tl->test().size());
        if (ln == 0) continue;
        if (ln == 1) {
            std::string name = tl->test(0)->getText();
            if (rhsVals.size() <= 1) {
                setVar(name, rhsVals.empty()? Value{} : rhsVals[0]);
            } else {
                // pack? choose first
                setVar(name, rhsVals[0]);
            }
        } else {
            // multiple names
            if (static_cast<int>(rhsVals.size()) == ln) {
                for (int k=0;k<ln;++k) {
                    std::string name = tl->test(k)->getText();
                    setVar(name, rhsVals[k]);
                }
            } else if (!rhsVals.empty()) {
                for (int k=0;k<ln;++k) {
                    std::string name = tl->test(k)->getText();
                    setVar(name, rhsVals[0]);
                }
            } else {
                for (int k=0;k<ln;++k) {
                    std::string name = tl->test(k)->getText();
                    setVar(name, Value{});
                }
            }
        }
    }
    return {};
}

// ===== Expressions =====
std::any EvalVisitor::visitTest(Python3Parser::TestContext* ctx) {
    return visit(ctx->or_test());
}

std::any EvalVisitor::visitOr_test(Python3Parser::Or_testContext* ctx) {
    auto parts = ctx->and_test();
    if (parts.empty()) return Value{};
    if (parts.size() == 1) return visit(parts[0]);
    bool res = to_bool(any_cast<Value>(visit(parts[0])));
    for (size_t i=1;i<parts.size();++i) {
        if (res) { res = true; break; }
        res = to_bool(any_cast<Value>(visit(parts[i])));
    }
    return Value{res};
}

std::any EvalVisitor::visitAnd_test(Python3Parser::And_testContext* ctx) {
    auto parts = ctx->not_test();
    if (parts.empty()) return Value{};
    if (parts.size() == 1) return visit(parts[0]);
    bool res = to_bool(any_cast<Value>(visit(parts[0])));
    for (size_t i=1;i<parts.size();++i) {
        if (!res) { res = false; break; }
        res = to_bool(any_cast<Value>(visit(parts[i])));
    }
    return Value{res};
}

std::any EvalVisitor::visitNot_test(Python3Parser::Not_testContext* ctx) {
    if (ctx->NOT()) {
        bool v = to_bool(any_cast<Value>(visit(ctx->not_test())));
        return Value{!v};
    }
    return visit(ctx->comparison());
}

std::any EvalVisitor::visitComparison(Python3Parser::ComparisonContext* ctx) {
    auto exprs = ctx->arith_expr();
    auto ops = ctx->comp_op();
    if (exprs.size()==1) return visit(exprs[0]);
    std::vector<Value> vals; vals.reserve(exprs.size());
    for (auto e: exprs) vals.push_back(any_cast<Value>(visit(e)));
    bool ok = true;
    for (size_t i=0;i<ops.size();++i) {
        int c = cmp(vals[i], vals[i+1]);
        auto op = ops[i];
        bool r = false;
        if (op->LESS_THAN()) r = (c==-1);
        else if (op->GREATER_THAN()) r = (c==1);
        else if (op->EQUALS()) r = (c==0);
        else if (op->GT_EQ()) r = (c==1 || c==0);
        else if (op->LT_EQ()) r = (c==-1 || c==0);
        else /*NOT_EQ_2*/ r = (c!=0);
        if (c==-2) { // incomparable
            if (op->EQUALS()) r=false; else if (op->NOT_EQ_2()) r=true; else r=false;
        }
        ok = ok && r;
        if (!ok) break;
    }
    return Value{ok};
}

std::any EvalVisitor::visitArith_expr(Python3Parser::Arith_exprContext* ctx) {
    auto terms = ctx->term();
    if (terms.empty()) return Value{};
    Value acc = any_cast<Value>(visit(terms[0]));
    for (size_t i=1;i<terms.size();++i) {
        auto op = ctx->addorsub_op(i-1);
        Value rhs = any_cast<Value>(visit(terms[i]));
        if (op->ADD()) acc = add(acc, rhs);
        else acc = sub(acc, rhs);
    }
    return acc;
}

std::any EvalVisitor::visitTerm(Python3Parser::TermContext* ctx) {
    auto facts = ctx->factor();
    if (facts.empty()) return Value{};
    Value acc = any_cast<Value>(visit(facts[0]));
    for (size_t i=1;i<facts.size();++i) {
        auto op = ctx->muldivmod_op(i-1);
        Value rhs = any_cast<Value>(visit(facts[i]));
        if (op->STAR()) acc = mul(acc, rhs);
        else if (op->DIV()) acc = truediv(acc, rhs);
        else if (op->IDIV()) acc = floordiv(acc, rhs);
        else acc = mod(acc, rhs);
    }
    return acc;
}

std::any EvalVisitor::visitFactor(Python3Parser::FactorContext* ctx) {
    if (ctx->atom_expr()) return visit(ctx->atom_expr());
    Value v = any_cast<Value>(visit(ctx->factor()));
    if (ctx->ADD()) return unary_plus(v);
    if (ctx->MINUS()) return unary_minus(v);
    return v;
}

std::any EvalVisitor::visitAtom_expr(Python3Parser::Atom_exprContext* ctx) {
    auto at = ctx->atom();
    auto tr = ctx->trailer();
    if (!tr) return visit(at);
    // function call: atom must be NAME
    if (!at->NAME()) return Value{};
    std::string fname = at->NAME()->getText();

    // Build args (positional and keyword)
    std::vector<std::pair<std::optional<std::string>, Value>> args; args.clear();
    if (tr->arglist()) {
        for (auto arg : tr->arglist()->argument()) {
            if (arg->ASSIGN()) {
                std::string an = arg->test(0)->getText();
                Value av = any_cast<Value>(visit(arg->test(1)));
                args.emplace_back(an, av);
            } else {
                Value av = any_cast<Value>(visit(arg->test(0)));
                args.emplace_back(std::nullopt, av);
            }
        }
    }

    // Builtins
    if (fname == "print" || fname == "int" || fname == "float" || fname == "str" || fname == "bool") {
        return call_builtin(fname, args);
    }

    // User-defined function
    auto it = functions_.find(fname);
    if (it == functions_.end()) {
        // calling an undefined name -> None
        return Value{};
    }
    const FunctionDef& f = it->second;
    std::unordered_map<std::string, Value> frame;
    frame.reserve(f.param_names.size());

    // First positional
    size_t posi = 0;
    for (auto& pr : args) {
        if (!pr.first) {
            if (posi < f.param_names.size()) {
                frame[f.param_names[posi]] = pr.second;
                ++posi;
            }
        }
    }
    // Then keywords
    for (auto& pr : args) {
        if (pr.first) {
            auto itn = std::find(f.param_names.begin(), f.param_names.end(), *pr.first);
            if (itn != f.param_names.end()) {
                frame[*itn] = pr.second;
            }
        }
    }
    // Fill defaults
    for (size_t i=0;i<f.param_names.size();++i) {
        const std::string& pn = f.param_names[i];
        if (frame.find(pn) == frame.end()) {
            if (i < f.default_tests.size() && f.default_tests[i]) {
                frame[pn] = any_cast<Value>(visit(f.default_tests[i]));
            } else {
                frame[pn] = Value{}; // None
            }
        }
    }

    // Execute in new scope where only parameters are local
    local_stack_.push_back(std::move(frame));
    Value retv = Value{};
    try {
        visit(f.suite);
    } catch (const ReturnSignal& rs) {
        if (rs.values.empty()) retv = Value{};
        else if (rs.values.size()==1) retv = rs.values[0];
        else retv = rs.values[0];
    }
    local_stack_.pop_back();
    return retv;
}

std::any EvalVisitor::visitTrailer(Python3Parser::TrailerContext* /*ctx*/) { return {}; }

std::any EvalVisitor::visitAtom(Python3Parser::AtomContext* ctx) {
    if (auto n = ctx->NAME()) {
        return getVar(n->getText());
    }
    if (auto num = ctx->NUMBER()) {
        std::string t = num->getText();
        if (t.find('.') != std::string::npos) {
            try { return Value{std::stod(t)}; } catch (...) { return Value{0.0}; }
        } else {
            try { return Value{std::stoll(t)}; } catch (...) { try { return Value{std::stod(t)}; } catch (...) { return Value{0LL}; } }
        }
    }
    if (ctx->NONE()) return Value{};
    if (ctx->TRUE()) return Value{true};
    if (ctx->FALSE()) return Value{false};
    if (ctx->format_string()) return visit(ctx->format_string());
    if (ctx->OPEN_PAREN()) {
        if (ctx->test()) return visit(ctx->test());
        return Value{};
    }
    // STRING (+)
    std::string res;
    for (auto s : ctx->STRING()) {
        res += strip_quotes(s->getText());
    }
    return Value{res};
}

#include "Python3Lexer.h"
#include "Python3Parser.h"

std::any EvalVisitor::visitFormat_string(Python3Parser::Format_stringContext* ctx) {
    // Rebuild from raw text and evaluate {...} by parsing the substring as a testlist
    std::string raw = ctx->getText(); // like f"..." or f'...'
    if (raw.size() < 3) return Value{std::string("")};
    size_t pos = 0;
    if (raw[pos] == 'f' || raw[pos] == 'F') ++pos;
    if (pos >= raw.size()) return Value{std::string("")};
    char quote = raw[pos++]; // ' or "
    if (pos >= raw.size()) return Value{std::string("")};
    // content excludes the trailing quote
    std::string content = raw.substr(pos, raw.size()-pos-1);
    std::string out;
    for (size_t i=0;i<content.size();) {
        if (content[i] == '{') {
            if (i+1 < content.size() && content[i+1] == '{') {
                out.push_back('{'); i+=2; continue;
            }
            // extract balanced {...}
            int depth = 1; size_t j = i+1;
            while (j < content.size() && depth > 0) {
                if (content[j] == '{') ++depth;
                else if (content[j] == '}') --depth;
                ++j;
            }
            size_t end = j; // position after '}'
            std::string expr = content.substr(i+1, end - (i+1) - 1); // between braces
            // Evaluate expr as testlist
            antlr4::ANTLRInputStream input(expr);
            Python3Lexer lexer(&input);
            antlr4::CommonTokenStream tokens(&lexer);
            tokens.fill();
            Python3Parser parser(&tokens);
            auto tree = parser.testlist();
            auto anyVals = visit(tree);
            if (auto p = any_cast<Values>(&anyVals)) {
                if (!p->empty()) out += to_string((*p)[0]);
            }
            i = end; // continue after '}'
        } else if (content[i] == '}') {
            if (i+1 < content.size() && content[i+1] == '}') {
                out.push_back('}'); i+=2; continue;
            }
            // lone '}' shouldn't happen; skip
            ++i;
        } else {
            out.push_back(content[i]); ++i;
        }
    }
    return Value{out};
}

std::any EvalVisitor::visitTestlist(Python3Parser::TestlistContext* ctx) {
    Values vs; vs.reserve(ctx->test().size());
    for (auto t : ctx->test()) vs.push_back(any_cast<Value>(visit(t)));
    return vs;
}

std::any EvalVisitor::visitArglist(Python3Parser::ArglistContext* /*ctx*/) { return {}; }
std::any EvalVisitor::visitArgument(Python3Parser::ArgumentContext* /*ctx*/) { return {}; }

// ===== Builtins =====
Value EvalVisitor::call_builtin(const std::string& name, const std::vector<std::pair<std::optional<std::string>, Value>>& args) {
    if (name == "print") {
        // print all args separated by space, then newline
        bool first = true;
        for (auto& pr : args) {
            const Value& v = pr.second;
            if (!first) std::cout << ' ';
            std::cout << to_string(v);
            first = false;
        }
        std::cout << '\n';
        return Value{};
    }
    if (name == "int") {
        if (!args.empty()) return Value{to_int(args[0].second)};
        return Value{0LL};
    }
    if (name == "float") {
        if (!args.empty()) return Value{to_float(args[0].second)};
        return Value{0.0};
    }
    if (name == "str") {
        if (!args.empty()) return Value{to_string(args[0].second)};
        return Value{std::string("")};
    }
    if (name == "bool") {
        if (!args.empty()) return Value{to_bool(args[0].second)};
        return Value{false};
    }
    return Value{};
}
