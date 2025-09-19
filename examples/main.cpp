#include <cassert>
#include <cmath>
#include "decorators.h"
#include "tools.h"

namespace SMOKE_PROBLEM {
// Test problem:
// Number of variables = 12
// Number of constraints = 3
// Number of equality constraints = 1
// H - identity matrix
using namespace QP_NNLS;
matrix_t H;
const std::vector<double> c = {2.17299e6, -13162.7, 1.36145e6, 557966.0, 1.46127e6, 6.23872e6, 24368.4, -24368.4, -14567.7, -9040.55, 46731.1, 36536.2};
const matrix_t A = {{1.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0},
                    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0},
                    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0}};
const std::vector<double> b = {1.0e-7, 1.0e-7, 1.0e-7};
const std::vector<double> lw = {-10.00553, -5.9893,  -0.0107047,  0.367332,  0.0107047,  0.0107047, -0.989295,  -10.0,  -1.0e-8,  -1.0e-8,  -1.0e-8,  1.0e-8};
const std::vector<double> up = {9998.4,    0.0107047, 9998.69,    9998.73,   0.189295,    0.189295,  0.0107047,  10.0,   50.0,    50.0,     50.0,     50.0};
}

namespace SMOKE_PROBLEM_BASELINE {
//solution
const std::vector<double> x = {-10.00553000112996, 0.01070470000195201, -0.01070469967089594, 9.994825300178494, 0.01070470036938787,
                        0.01070469990372658, -0.9892949999994021, 9.99789335764945e-08, 50.00000000000912, 50.00000000000366,
                        -9.989889804273847e-09, 9.989889804273847e-09};
const std::vector<double> dualA = {-557975.9948253002, 582344.3948252, 0.0};
const std::vector<double> dualL = {1615003.999644698, 0.0, 803473.9944700002, 0.0, 903294.0158794001, 5680744.0158794, 24367.410705,
                        0.0, 0.0, 0.0, 46731.09999999001, 36536.20000000999};
const std::vector<double> dualU = {0.0, 13162.68929529999, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 14517.69999999998, 8990.549999999992, 0.0, 0.0};
const double cost = -17299352.22259864;
const std::size_t nIterations = 12;
}
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
    /*std::vector<std::vector<double>> ATA(nCols, std::vector<double>(nCols, 0.0));
    for (std::size_t i = 0; i < nCols; ++i) {
        for (std::size_t j = 0; j < nCols; ++j) {
            for (std::size_t k = 0; k < nRows; ++k) {
                ATA[i][j] += A[k][j] * A[k][i];
            }
        }
    }*/
    std::vector<double> b(nRows, 0.1);
    std::vector<double> c(nCols, 0.0);
    for (int i = 0; i < nCols; ++i) {
        for (int j = 0; j < nRows; ++j) {
           c[i] -= A[j][i] * b[j];
        }
    }
    using namespace QP_NNLS;
    DenseQPProblem problem;
    problem.H = std::move(H);
    problem.c = std::move(c);
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

    // check output
    /*
    const double eps = sqrt(GetMachineEps());
    if (success) {
        if ((std::fabs(nnlsOutput.cost - SMOKE_PROBLEM_BASELINE::cost) >= eps) ||
            nnlsOutput.nDualIterations != SMOKE_PROBLEM_BASELINE::nIterations) {
            std::cout << "FAILED" << std::endl;
            return 1;
        }
        for (std::size_t v = 0; v < nVariables; ++v) {
            if ((std::fabs(nnlsOutput.x[v] - SMOKE_PROBLEM_BASELINE::x[v]) >= eps) ||
                (std::fabs(nnlsOutput.lambdaUp[v] - SMOKE_PROBLEM_BASELINE::dualU[v]) >= eps) ||
                (std::fabs(nnlsOutput.lambdaLw[v] - SMOKE_PROBLEM_BASELINE::dualL[v]) >= eps)) {
                    std::cout << "FAILED" << std::endl;
                    return 1;
            }
        }
        for (std::size_t c = 0; c < nConstraints; ++c) {
            if (std::fabs(nnlsOutput.lambda[c] - SMOKE_PROBLEM_BASELINE::dualA[c]) >= eps) {
                std::cout << "FAILED" << std::endl;
                return 1;
            }
        }
        std::cout << "PASSED" << std::endl;
        return 0;
    } else {
        std::cout << "FAILED" << std::endl;
        return 1;
    }*/
}
