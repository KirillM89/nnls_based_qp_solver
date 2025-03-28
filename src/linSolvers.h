#ifndef LINSOLVERS_H
#define LINSOLVERS_H
#include "types.h"
#include "utils.h"
namespace QP_NNLS {
class ILinSolver {
    // Interface for linear solver
    // [M s] * [M_T s_T] * y = - gamma * s
    // methods Add() and Delete() calls when corresponding constraint
    // adds/deletes to/from active set
    // Solve() calls in place where the problem has to be solved
public:
    virtual ~ILinSolver() = default;
    virtual bool Add(unsg_t indx) = 0;
    virtual bool Delete(unsg_t indx) = 0;
    virtual void SetGamma(double gamma) = 0;
    virtual void Reset() {};
    virtual void Set(bool rejectSingular) {};
    virtual const LinSolverOutput& Solve() = 0;
protected:
    ILinSolver() = default;
};

class CumulativeSolver: public ILinSolver {
    // Add / Delete methods constructs linear system
    // Solve() solves pre-constructed linear system
public:
    CumulativeSolver() = delete;
    CumulativeSolver(const matrix_t& M, const std::vector<double>& s);
    virtual ~CumulativeSolver() override = default;
    virtual bool Add(unsg_t indx) override;
    virtual bool Delete(unsg_t indx) override;
    virtual void SetGamma(double gamma) override {this->gamma = gamma;}
protected:
    const unsg_t nConstraints;
    unsg_t nVariables;
    unsg_t nActive;
    double gamma;
    std::vector<bool> activeSet;
    const matrix_t& M;
    const std::vector<double>& s;
    LinSolverOutput output;

};

class CumulativeLDLTSolver: public ILinSolver {
    // Solve linear system using custom LDLT decomposition
public:
    CumulativeLDLTSolver() = delete;
    CumulativeLDLTSolver(const matrix_t& M, const std::vector<double>& s, bool rejectSingular);
    virtual ~CumulativeLDLTSolver() override = default;
    const LinSolverOutput& Solve() override;
    virtual bool Add(unsg_t indx) override;
    virtual bool Delete(unsg_t indx) override;
    virtual void SetGamma(double gamma) override {this->gamma = -gamma;}
    virtual void Set(bool rejectSingular) override { this->rejectSingular = rejectSingular;}
protected:
    bool rejectSingular;
    const double zeroTol = g_GetMachineEps();
    double gamma;
    LDLT ldlt;
    unsigned int ndzero;
    const std::size_t maxSize;
    const std::vector<double>& S;
    std::vector<double> forward;
    std::vector<double> backward;
    LinSolverOutput output;
};

}
#endif // LINSOLVERS_H
