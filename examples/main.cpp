#include <cassert>
#include <cmath>
#include <iostream>
#include "decorators.h"
#include "tools.h"

// Test problem:
// Number of variables = 12
// Number of constraints = 3
// Number of equality constraints = 1
// H - identity matrix

double GetMachineEps() {
    double eps = 1.0;
    while (eps + 1.0 > 1.0) {
        eps *= 0.5;
    }
    return 2 * eps;
}

int main() {
    // check correctness of input and baseline
    const std::string fileH = "C:/Users/m00829527/nnls_based_qp_solver/examples/data/ata_1.txt";
    const std::string fileA = "C:/Users/m00829527/nnls_based_qp_solver/examples/data/low_rank_1.txt";
    const std::string logF = "C:/Users/m00829527/nnls_based_qp_solver/examples/data/log1.txt";
    auto H = readMatrix(fileH);
    auto A = readMatrix(fileA);

    std::size_t nRows = H.size();
    std::size_t nCols = H.front().size();
    assert(nRows == nCols);
    std::vector<double> b(A.size(), 0.1);
    std::vector<double> c(nCols, 0.0);
    for (int i = 0; i < nCols; ++i) {
        for (int j = 0; j < A.size(); ++j) {
           c[i] -= A[j][i] * b[j];
        }
    }
    using namespace QP_NNLS;
    DenseQPProblem problem;
    problem.H = std::move(H);
    problem.c = std::move(c);
    problem.A = std::vector<std::vector<double>>(1, std::vector<double>(nCols, 0.0));
    A.front().front() = -1.0;
    problem.b = {0.0};
    problem.up = std::vector<double>(nCols, 10.0);
    problem.lw = std::vector<double>(nCols, 0.05);
    Settings settings;     // default settings
    settings.logLevel = 5u;
    QPNNLSDense solver;    // create solver instance
    std::unique_ptr<Callback> cb = std::make_unique<Callback1>(logF);
    solver.SetCallback(std::move(cb));
    solver.Init(settings);
    // set problem
    // method SetProblem() preprocesses the input problem.
    // Returns true if all the preprocessing procedures passed without numerical problems, false otherwise
    // If SetProblem() returns "true" the properties of the problem can be seen using method GetInitStatus()
    // QP_NNLS::InitStageStatus initStatus = solver.GetInitStatus()
    // "SUCCESS" means that problem has positive definite H
    // "D_Z" means positive semidefinite
    // "D_ZN" or "D_N"  - indefinite or negative definite
    // Solver can't solve indefinite or negative definite problems but can solve positive semidefinite problems,
    // in this case H correction will take place
    std::string initMsg;
    if (!solver.SetProblem(problem)) {
        QP_NNLS::InitStageStatus initStatus = solver.GetInitStatus();
        if (initStatus == QP_NNLS::InitStageStatus::D_Z) {
            initMsg = "LDLT: zero D components";
        } else if (initStatus == QP_NNLS::InitStageStatus::D_N) {
            initMsg = "LDLT: negative D components";
        } else if (initStatus == QP_NNLS::InitStageStatus::D_ZN) {
            initMsg = "LDLT: zero and negative D components";
        }
        std::cout << initMsg << std::endl;
        return 1;
    }
    solver.Solve();
    QP_NNLS::SolverOutput nnlsOutput = solver.GetOutput();
    bool success = true;
    if (nnlsOutput.dualExitStatus == QP_NNLS::DualLoopExitStatus::INFEASIBILITY) {
        std::cout << "infeasibility" << std::endl;
        success = false;
    } else if (nnlsOutput.dualExitStatus == QP_NNLS::DualLoopExitStatus::ITERATIONS) {
        std::cout << "iterations" << std::endl;
    } else if (nnlsOutput.dualExitStatus == QP_NNLS::DualLoopExitStatus::ALL_DUAL_POSITIVE)  {
        std::cout << "all dual positive" << std::endl;
    } else if (nnlsOutput.dualExitStatus == QP_NNLS::DualLoopExitStatus::FULL_ACTIVE_SET) {
        std::cout <<  "full active set" << std::endl;
    } else {
        abort();
    }
    return 0;
}

