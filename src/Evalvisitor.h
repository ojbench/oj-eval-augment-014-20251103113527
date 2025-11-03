#pragma once
#ifndef PYTHON_INTERPRETER_EVALVISITOR_H
#define PYTHON_INTERPRETER_EVALVISITOR_H

#include <any>
#include <iostream>
#include <optional>
#include <string>
#include <unordered_map>
#include <utility>
#include <variant>
#include <vector>

#include "Python3ParserBaseVisitor.h"

// A simple dynamic value type: None | bool | int64 | double | string
using Value = std::variant<std::monostate, bool, long long, double, std::string>;
using Values = std::vector<Value>;

struct ReturnSignal {
    Values values; // empty -> None, size==1 -> single
};
struct BreakSignal {};
struct ContinueSignal {};

struct FunctionDef {
    std::vector<std::string> param_names;
    // default_tests.size() == param_names.size(); nullptr means no default
    std::vector<Python3Parser::TestContext*> default_tests;
    Python3Parser::SuiteContext* suite = nullptr;
};

class EvalVisitor : public Python3ParserBaseVisitor {
public:
    EvalVisitor() = default;

    // Top-level entry
    std::any visitFile_input(Python3Parser::File_inputContext* ctx) override;
    std::any visitFuncdef(Python3Parser::FuncdefContext* ctx) override;
    std::any visitParameters(Python3Parser::ParametersContext* ctx) override;
    std::any visitTypedargslist(Python3Parser::TypedargslistContext* ctx) override;
    std::any visitTfpdef(Python3Parser::TfpdefContext* ctx) override;


    // Statements
    std::any visitStmt(Python3Parser::StmtContext* ctx) override;
    std::any visitSimple_stmt(Python3Parser::Simple_stmtContext* ctx) override;
    std::any visitSmall_stmt(Python3Parser::Small_stmtContext* ctx) override;
    std::any visitExpr_stmt(Python3Parser::Expr_stmtContext* ctx) override;
    std::any visitAugassign(Python3Parser::AugassignContext* ctx) override;
    std::any visitFlow_stmt(Python3Parser::Flow_stmtContext* ctx) override;
    std::any visitBreak_stmt(Python3Parser::Break_stmtContext* ctx) override;
    std::any visitContinue_stmt(Python3Parser::Continue_stmtContext* ctx) override;
    std::any visitReturn_stmt(Python3Parser::Return_stmtContext* ctx) override;
    std::any visitCompound_stmt(Python3Parser::Compound_stmtContext* ctx) override;
    std::any visitIf_stmt(Python3Parser::If_stmtContext* ctx) override;
    std::any visitWhile_stmt(Python3Parser::While_stmtContext* ctx) override;
    std::any visitSuite(Python3Parser::SuiteContext* ctx) override;

    // Expressions
    std::any visitTest(Python3Parser::TestContext* ctx) override;
    std::any visitOr_test(Python3Parser::Or_testContext* ctx) override;
    std::any visitAnd_test(Python3Parser::And_testContext* ctx) override;
    std::any visitNot_test(Python3Parser::Not_testContext* ctx) override;
    std::any visitComparison(Python3Parser::ComparisonContext* ctx) override;
    std::any visitArith_expr(Python3Parser::Arith_exprContext* ctx) override;
    std::any visitTerm(Python3Parser::TermContext* ctx) override;
    std::any visitFactor(Python3Parser::FactorContext* ctx) override;
    std::any visitAtom_expr(Python3Parser::Atom_exprContext* ctx) override;
    std::any visitTrailer(Python3Parser::TrailerContext* ctx) override;
    std::any visitAtom(Python3Parser::AtomContext* ctx) override;
    std::any visitFormat_string(Python3Parser::Format_stringContext* ctx) override;
    std::any visitTestlist(Python3Parser::TestlistContext* ctx) override;
    std::any visitArglist(Python3Parser::ArglistContext* ctx) override;
    std::any visitArgument(Python3Parser::ArgumentContext* ctx) override;

private:
    // Environment management
    std::unordered_map<std::string, Value> globals_;
    std::vector<std::unordered_map<std::string, Value>> local_stack_;
    std::unordered_map<std::string, FunctionDef> functions_;

    std::unordered_map<std::string, Value>& current_locals();
    bool in_function() const { return !local_stack_.empty(); }

    // Variable access consistent with problem's scope rules
    Value getVar(const std::string& name) const;
    void setVar(const std::string& name, const Value& v);

    // Builtins
    Value call_builtin(const std::string& name, const std::vector<std::pair<std::optional<std::string>, Value>>& args);

    // Helpers: conversions and operators
    static bool is_none(const Value& v) { return std::holds_alternative<std::monostate>(v); }
    static bool as_bool_strict(const Value& v);
    static bool to_bool(const Value& v);
    static long long to_int(const Value& v);
    static double to_float(const Value& v);
    static std::string to_string(const Value& v);

    static Value add(const Value& a, const Value& b);
    static Value sub(const Value& a, const Value& b);
    static Value mul(const Value& a, const Value& b);
    static Value truediv(const Value& a, const Value& b);
    static Value floordiv(const Value& a, const Value& b);
    static Value mod(const Value& a, const Value& b);
    static Value unary_plus(const Value& a);
    static Value unary_minus(const Value& a);

    static int cmp(const Value& a, const Value& b); // like Python compare after conversions; string only with string

    // Utilities
    static std::string strip_quotes(const std::string& s);
};

#endif // PYTHON_INTERPRETER_EVALVISITOR_H
