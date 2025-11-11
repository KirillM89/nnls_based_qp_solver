#include "qp_solvers.h"
#include "qld.h"
#include "decorators.h"
#include <cmath>
#include <cassert>
namespace QP_SOLVERS {
bool QLD::SetProblem(const DenseQpProblem& problem)  {
    nVars = static_cast<int>(problem.H.size());
    nConstraints = static_cast<int>(problem.A.size());
    nEqConstraints = problem.nEqConstraints;
    if (nEqConstraints > nConstraints) {
        return false;
    }
    A = matrix_t(nConstraints, std::vector<double>(nVars, 0.0));
    for (int i = 0; i < nConstraints; ++i) {
        for (int j = 0; j < nVars; ++j) {
            A[i][j] = -problem.A[i][j];
        }
    }
    aFlat = FlattenByColums(A);
    hFlat = FlattenByColums(problem.H);
    H = problem.H;
    c = problem.c;
    b = problem.b;
    lw = problem.lw;
    up = problem.up;
    return true;
}
void QLD::Solve(SolverOutput& output) {
    const double machinePrecision = 1e-16;
    int mPlus2N = nConstraints + 2 * nVars;
    int iOut = 6;
    int zero = 0;
    int status = 0;
    int workArraySize = static_cast<int>(1.5 * nVars  * nVars)  + 10 * nVars  + 2 * nConstraints + 1;
    std::vector<int> intArray(nVars, 0.0);
    std::vector<double> workArray(workArraySize);
    std::vector<double> lambda(mPlus2N, 0.0);
    std::vector<double> x(nVars, 0.0);
    // The following asks QP-solver to use LDL-decomposition (square-root-free Cholesky decomposition)
    // C = L*D*L^T, where L - lower triangular, D - diagonal
    intArray[0] = 1;  // use LDL-decomposition
    TIMER::wcTimer timer;
    timer.Start();
    ql0001_(&nConstraints, &nEqConstraints, &nConstraints, &nVars,
            &nVars, &mPlus2N, hFlat.data(), c.data(), aFlat.data(), b.data(), lw.data(), up.data(),
            x.data(), lambda.data(), &iOut, &status, &zero, workArray.data(), &workArraySize, intArray.data(),
            &nVars, &machinePrecision);
    timer.Stop();
    output.time = timer.TimeFromFirstStart();
    double cost = 0.0;
    const double half = 0.5;
    for (std::size_t i = 0; i < nVars; ++i) {
        for (std::size_t j = 0; j <= i; ++j) {
            cost += H[i][j] * x[i] * x[j] * (i == j ? half : 1.0);
        }
        cost += c[i] * x[i];
    }
    const std::string qpReturnMessages[7] = {"", "Termination after too many iterations: 40*(N+M)",
    "Insufficient accuracy", "Inconsistency, division by zero", "Numerical instability",
    "Incorrect epsilon or working/int array sizes", "Inconsistent constraint"};
    if (status > 100) {
        output.msg =  qpReturnMessages[6] + std::to_string(status - 100);
        output.status = SolverRetStatus::INFEASIBILITY;
    } else if (status > 0 && status < 6) {
        output.msg =  qpReturnMessages[status];
        output.status = SolverRetStatus::NUM_PROBLEMS;
    } else if (status == 0) {
        output.msg.clear();
        output.status = SolverRetStatus::SUCCESS;
    }
    if (status == 0) {
        output.primal = x;
        output.cost = cost;
        output.dualA = std::vector<double>(lambda.begin(), lambda.begin() + nConstraints);
        output.dualL = std::vector<double>(lambda.begin() + nConstraints, lambda.begin() + nConstraints + nVars);
        output.dualU = std::vector<double>(lambda.begin() + nConstraints + nVars, lambda.end());
    }
}

std::vector<double> QLD::FlattenByColums(const matrix_t& mat) {
    std::size_t m = mat.size();
    std::size_t n = m ? mat[0].size() : 0;
    std::vector<double> matrixByColumns(m * n);
    for (std::size_t j = 0; j < n; ++j) {
        for (std::size_t i = 0; i < m; ++i) {
            matrixByColumns[m * j + i] = mat[i][j];
        }
    }
    return matrixByColumns;
}

bool NNLS::SetProblem(const DenseQpProblem& problem)  {
    QP_NNLS::Settings settings;
    timer.Start();
    solver.Init(settings);
    timer.Stop();
    initMsg.clear();
    QP_NNLS::DenseQPProblem nnlsProblem;
    nnlsProblem.H = problem.H;
    nnlsProblem.A = problem.A;
    nnlsProblem.c = problem.c;
    nnlsProblem.b = problem.b;
    nnlsProblem.lw = problem.lw;
    nnlsProblem.up = problem.up;
    nnlsProblem.nEqConstraints = problem.nEqConstraints;
    if (!solver.SetProblem(nnlsProblem)) {
        QP_NNLS::InitStageStatus initStatus = solver.GetInitStatus();
        if (initStatus == QP_NNLS::InitStageStatus::D_Z) {
            initMsg = "LDLT: zero D components";
        } else if (initStatus == QP_NNLS::InitStageStatus::D_N) {
            initMsg = "LDLT: negative D components";
        } else if (initStatus == QP_NNLS::InitStageStatus::D_ZN) {
            initMsg = "LDLT: zero and negative D components";
        }
    }
    return initMsg.empty();
}
void NNLS::Solve(SolverOutput& output) {
    timer.Start();
    solver.Solve();
    timer.Stop();
    output.time = timer.TimeFromFirstStart();
    timer.Reset();
    QP_NNLS::SolverOutput nnlsOutput = solver.GetOutput();
    output.msg.clear();
    bool success = false;
    if (!initMsg.empty()) {
        output.msg = initMsg;
        output.status = SolverRetStatus::PREPROCESSING;
    } else if (nnlsOutput.dualExitStatus == QP_NNLS::DualLoopExitStatus::INFEASIBILITY) {
        output.msg = "infeasibility";
        output.status = SolverRetStatus::INFEASIBILITY;
    } else if (nnlsOutput.dualExitStatus == QP_NNLS::DualLoopExitStatus::ITERATIONS) {
        output.msg = "iterations";
        output.status = SolverRetStatus::ITERATIONS;
        success = true;
    } else if (nnlsOutput.dualExitStatus == QP_NNLS::DualLoopExitStatus::UNKNOWN){
        output.msg = "unknown exit status";
    } else if (nnlsOutput.dualExitStatus == QP_NNLS::DualLoopExitStatus::ALL_DUAL_POSITIVE)  {
        success = true;
        output.msg = "all dual positive";
        output.status = SolverRetStatus::SUCCESS;
    } else if (nnlsOutput.dualExitStatus == QP_NNLS::DualLoopExitStatus::FULL_ACTIVE_SET) {
        success = true;
        output.msg = "full active set";
        output.status = SolverRetStatus::SUCCESS;
    } else {
        abort();
    }
    if (success) {
        output.nIterations = nnlsOutput.nDualIterations;
        output.primal = nnlsOutput.x;
        output.dualA = nnlsOutput.lambda;
        output.dualL = nnlsOutput.lambdaLw;
        output.dualU = nnlsOutput.lambdaUp;
        output.cost = nnlsOutput.cost;
    }
}

bool DAQP::SetProblem(const DenseQpProblem& problem) {
    nVariables = problem.H.size();
    nEqConstraints = problem.nEqConstraints;
    assert((problem.A.size() - nEqConstraints) % 2 == 0);
    nConstraints = (problem.A.size() - nEqConstraints) / 2 + nEqConstraints;
    h = std::vector<double>(nVariables * nVariables);
    a = std::vector<double>(nVariables * nConstraints);
    c = problem.c;
    const std::size_t nTotalConstraints = nConstraints + nVariables;
    lb = std::vector<double>(nTotalConstraints);
    ub = std::vector<double>(nTotalConstraints);
    constraintsType = std::vector<int>(nTotalConstraints ,0);
    for (std::size_t i = 0; i < std::max(nVariables, nConstraints); i++) {
        if (i < nVariables) {
            lb[i] = std::fmax(problem.lw[i], -defaultBoundValue);
            ub[i] = std::fmin(problem.up[i], defaultBoundValue);
        }
        for (std::size_t j = 0; j < nVariables; j++) {
            if (i < nVariables) {
                h[i * nVariables + j] = problem.H[i][j];
            }
            if (i < nEqConstraints) {
                a[i * nVariables + j] = problem.A[i][j];
            } else if (i < nConstraints) {
                a[i * nVariables + j] = problem.A[2 * i - nEqConstraints][j];
            }
        }
    }
    for (std::size_t i = 0; i < nConstraints; ++i) {
        const std::size_t indx = nVariables + i;
        if (i < nEqConstraints) {
            ub[indx] = problem.b[i];
            lb[indx] = ub[indx];
            constraintsType[indx] = 5;
        } else {
            lb[indx] = -problem.b[2 * i - nEqConstraints + 1];
            ub[indx] = problem.b[2 * i - nEqConstraints];
        }
    }
    daqp_default_settings(&settings);
    settings.iter_limit = 10000;   /*10000 denotes total number of iterations*/
    settings.cycle_tol = 1000;     /*1000 denotes allowed number of iterations without progress before terminating*/
    settings.progress_tol = 1e-12; /*1e-12 denotes minimum change in objective function to consider it progress*/
    return true;
}

void DAQP::Solve(SolverOutput& output) {
    double sol[nVariables];
    double lam[nConstraints + nVariables];
    for (int i = 0; i < nConstraints + nVariables; ++i) {
        if (i <  nVariables) {
            sol[i] = 0.0;
        }
        lam[i] = 0.0;
    }
    DAQPResult result;
    result.x = sol;   // primal variable
    result.lam = lam; // dual variable
    DAQPProblem qp = {nVariables, nConstraints + nVariables, nVariables, h.data(), c.data(),
                      a.data(), ub.data(), lb.data(), constraintsType.data()};
    TIMER::wcTimer timer;
    timer.Start();
    daqp_quadprog(&result, &qp, &settings);
    timer.Stop();
    output.time = timer.TimeFromFirstStart();
    if (result.exitflag == efInfeasible) {
        output.status = SolverRetStatus::INFEASIBILITY;
        output.msg = "infeasible";
    } else if (result.exitflag == efCycling) {
        output.status = SolverRetStatus::CYCLING;
        output.msg = "cycling";
    } else if (result.exitflag == efUnbounded) {
        output.status = SolverRetStatus::INFEASIBILITY;
        output.msg = "unbounded";
    } else if (result.exitflag == efNonConvex) {
        output.status = SolverRetStatus::NONCONVEX;
        output.msg = "nonconvex";
    } else if (result.exitflag == efInitSet) {
        output.status = SolverRetStatus::OTHER;
        output.msg = "init active set";
    } else if (result.exitflag == -4 || result.exitflag == 1 || result.exitflag == 2) {
        //success
        output.status = SolverRetStatus::SUCCESS;
        output.msg.clear();
        output.primal = std::vector<double>(result.x, result.x + nVariables);
        output.dualA = std::vector<double>(2 * nConstraints - nEqConstraints, 0.0);
        output.dualL = std::vector<double>(nVariables, 0.0);
        output.dualU = std::vector<double>(nVariables, 0.0);
        for (std::size_t i = 0; i < nVariables; ++i) {
            if (std::fabs(result.lam[i]) >= zeroTol && (std::fabs(output.primal[i] - ub[i]) <= zeroTol)) {
                output.dualU[i] = result.lam[i];
            } else {
                output.dualL[i] = -result.lam[i];
            }
        }
        for (std::size_t i = 0; i < nEqConstraints; ++i) {
            output.dualA[i] = result.lam[i + nVariables];
        }
        for (std::size_t i = nEqConstraints; i < nConstraints; ++i) {
            double ax = 0.0;
            for (std::size_t j = 0; j < nVariables; ++j) {
                ax += a[i * nVariables + j] * result.x[j];
            }
            const std::size_t indx = nVariables + i;
            if (std::fabs(result.lam[indx]) > zeroTol && (std::fabs(ax - ub[indx]) <= zeroTol)) {
                output.dualA[2 * i - nEqConstraints] = result.lam[indx];
            } else {
                output.dualA[2 * i - nEqConstraints + 1] = -result.lam[indx];
            }
        }
        output.nIterations = result.iter;
        output.cost = result.fval;
    } else {
        output.msg = "abnormal exit flag";
    }
}

} //namespace QP_SOLVERS
