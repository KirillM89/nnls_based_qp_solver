#ifndef QP_SOLVERS_H
#define QP_SOLVERS_H
#include "qp_types.h"
#include "decorators.h"
#include "api.h"
#include "timer.h"
namespace QP_SOLVERS {
class iQpSolver {
public:
    virtual ~iQpSolver() = default;
    virtual bool SetProblem(const DenseQpProblem& problem) = 0;
    virtual void Solve(SolverOutput& output) = 0;
    virtual const std::string& GetInitMsg() = 0;
protected:
    iQpSolver() = default;
};

class QLD : public iQpSolver {
public:
    QLD() = default;
    virtual ~QLD() override = default;
    bool SetProblem(const DenseQpProblem& problem) override;
    void Solve(SolverOutput& output) override;
    const std::string& GetInitMsg() override {return initMsg;}
protected:
    std::vector<double> FlattenByColums(const matrix_t& m);
    int nVars = 0;
    int nConstraints = 0;
    int nEqConstraints = 0;
    const double machinePrecision = 1e-16;
    std::string initMsg;
    matrix_t H;
    std::vector<double> hFlat;
    std::vector<double> aFlat;
    std::vector<double> c;
    std::vector<double> b;
    std::vector<double> lw;
    std::vector<double> up;
    matrix_t A;

};

class NNLS : public iQpSolver {
public:
    NNLS() = default;
    virtual ~NNLS() override = default;
    bool SetProblem(const DenseQpProblem& problem) override;
    void Solve(SolverOutput& output) override;
    const std::string& GetInitMsg() override {return initMsg;}
protected:
    QP_NNLS::QPNNLSDense solver;
    std::string initMsg;
    TIMER::wcTimer timer;

};

class DAQP: public iQpSolver {
public:
    DAQP() = default;
    virtual ~DAQP() override = default;
    bool SetProblem(const DenseQpProblem& problem) override;
    void Solve(SolverOutput& output) override;
    const std::string& GetInitMsg() override {return initMsg;}
protected:
    const double defaultBoundValue = 1.0e19;
    const double zeroTol = 1.0e-10;
    const int efSoftOptimal = 2;
    const int efOptimal = 1;
    const int efInfeasible = -1;
    const int efCycling = -2;
    const int efUnbounded = -3;
    const int efIterations = -4;
    const int efNonConvex = -5;
    const int efInitSet = -6;

    QP_NNLS::QPNNLSDense solver;
    std::string initMsg;
    std::size_t nVariables = 0;
    std::size_t nConstraints = 0;
    std::size_t nEqConstraints = 0;
    std::vector<double> c;
    std::vector<double> b;
    std::vector<double> h;
    std::vector<double> a;
    std::vector<double> lb;
    std::vector<double> ub;
    std::vector<int> constraintsType;
    DAQPSettings settings;

};

}
#endif // QP_SOLVERS_H
