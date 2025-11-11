#ifndef SOLVER_COMPARATOR_H
#define SOLVER_COMPARATOR_H
#include <memory>
#include <vector>
#include "qp_solvers.h"
#include "data_writer.h"

bool g_isSame(double a, double b, double eps = 1.0e-14);
namespace SOLVER_COMPARATOR {
class LogGenerator {
public:
    LogGenerator() = default;
    ~LogGenerator() = default;
    void SetLogFile(const std::string& file);
    void Dump(const std::string& problemName, const QP_SOLVERS::QpMetrics& output);
    void Err(const std::string& err);
private:
     FMT_WRITER::FmtWriter logger;
     std::string file;
};

class QualityChecker {
public:
    QualityChecker() = default;
    ~QualityChecker() = default;
    void Check(const QP_SOLVERS::SolverOutput& output,
               const QP_SOLVERS::DenseQpProblem& problem);
    const QP_SOLVERS::QpMetrics& GetMetrics() { return metrics;}
private:
    const QP_SOLVERS::SolverOutput* output = nullptr;
    const QP_SOLVERS::DenseQpProblem* problem = nullptr;
    const double zeroPlus = 1.0e-20;
    const double zeroMinus = -zeroPlus;
    double xHx = 0.0;
    double cTx = 0.0;
    double bTL = 0.0;
    double lTL = 0.0;
    double uTL = 0.0;
    unsigned int nEqConstraints = 0;
    unsigned int nIEqConstraints = 0;
    unsigned int nVariables = 0;
    QP_SOLVERS::QpMetrics metrics;
    void CheckPrimalFeasibility();
    void CheckDualFeasibility();
    void CheckDualityGap();
};

class QpProblemSolver {
public:
    QpProblemSolver() = default;
    ~QpProblemSolver() = default;
    void SetQpSolver(QP_SOLVERS::QpSolverType solverType);
    bool SetProblem(const QP_SOLVERS::DenseQpProblem& problem);
    void Solve(QP_SOLVERS::SolverOutput& output);
    void CheckQuality(const QP_SOLVERS::SolverOutput& output,
                      QP_SOLVERS::QpMetrics& metrics);
protected:
    std::string msg;
    QP_SOLVERS::DenseQpProblem problem;
    std::unique_ptr<QP_SOLVERS::iQpSolver> solver;
    QualityChecker qlChecker;
};

class ProblemCreator {
public:
    ProblemCreator() = default;
    ~ProblemCreator() = default;
    bool ConfigureFromPath(const std::string& filePath);
    const QP_SOLVERS::DenseQpProblem& GetProblem();
    std::string errMsg;
private:
    QP_SOLVERS::DenseQpProblem problem;
};

void TestQpSolver(QP_SOLVERS::QpSolverType qpType, unsigned char cType, unsigned char num);
}
#endif // SOLVER_COMPARATOR_H
