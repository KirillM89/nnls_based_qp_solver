#include <cassert>
#include <cmath>
#include "decorators.h"

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
    // check correctness of input and baseline data
    const std::size_t nVariables = 12;
    assert(SMOKE_PROBLEM::c.size() == nVariables);
    assert(SMOKE_PROBLEM::lw.size() == nVariables);
    assert(SMOKE_PROBLEM::up.size() == nVariables);
    assert(SMOKE_PROBLEM_BASELINE::x.size() == nVariables);
    assert(SMOKE_PROBLEM_BASELINE::dualL.size() == nVariables);
    assert(SMOKE_PROBLEM_BASELINE::dualU.size() == nVariables);
    SMOKE_PROBLEM::H.resize(nVariables, std::vector<double>(nVariables, 0.0));
    for (std::size_t v = 0; v < nVariables; ++v) {
        assert( SMOKE_PROBLEM::lw[v] < SMOKE_PROBLEM::up[v]);
        SMOKE_PROBLEM::H[v][v] = 1.0;
    }
    const std::size_t nConstraints = 3;
    assert(SMOKE_PROBLEM::A.size() == nConstraints);
    assert(SMOKE_PROBLEM::b.size() == nConstraints);
    assert(SMOKE_PROBLEM_BASELINE::dualA.size() == nConstraints);
    for (const auto& constraint: SMOKE_PROBLEM::A) {
        assert(constraint.size() == nVariables);
    }
    // initialize problem
    QP_NNLS::DenseQPProblem nnlsProblem;
    nnlsProblem.H = SMOKE_PROBLEM::H;
    nnlsProblem.A = SMOKE_PROBLEM::A;
    nnlsProblem.c = SMOKE_PROBLEM::c;
    nnlsProblem.b = SMOKE_PROBLEM::b;
    nnlsProblem.lw = SMOKE_PROBLEM::lw;
    nnlsProblem.up = SMOKE_PROBLEM::up;
    nnlsProblem.nEqConstraints = 1; // first nEqConstraints rows in A are equality constraints
    QP_NNLS::Settings settings;     // default settings
    QP_NNLS::QPNNLSDense solver;    // create solver instance
    solver.Init(settings);          // apply settings
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
    if (!solver.SetProblem(nnlsProblem)) {
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
    }
}
