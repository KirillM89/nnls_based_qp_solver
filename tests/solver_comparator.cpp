#include <cmath>
#include "solver_comparator.h"
#include "TxtParser.h"
#include "configuration.h"
namespace SOLVER_COMPARATOR {
bool g_isSame(double a, double b, double eps = 1.0e-14) {
    if (eps < 0.0) {
        eps = -eps;
    }
    const double diff = a - b;
    return ((-eps <= diff) && (diff <= eps));
}
void LogGenerator::SetLogFile(const std::string& file) {
    this->file = file;
    logger.SetFile(file);
}
void LogGenerator::Dump(const std::string& problemName, const QP_SOLVERS::QpMetrics& result) {
    using namespace FMT_WRITER;
    using namespace FMT_WRITER;
    if (logger.LineNumber() % 20 == 0) {
        logger.Write("test name",
                     "time",
                     "matrix def",
                     "n variables",
                     "n constraints",
                     "n eq constraints",
                     "n iterations",
                     "status",
                     "max pr infsb ineq cnstr",
                     "violated value",
                     "n violated",
                     "max pr infsb eq cnstr",
                     "violated value",
                     "n violated",
                     "max pr infsb bnds",
                     "violated value",
                     "n violated",
                     "primal cost",
                     "duality gap",
                     "max dual infsb",
                     "max negative dual");
        logger.NewLine();
    }
    logger.Write(problemName,
                 result.time,
                 result.mDef,
                 result.nVariables,
                 result.nConstraints,
                 result.nEqConstraints,
                 result.nIterations);
    if (!result.errMsg.empty()) {
        logger.Write(result.errMsg);
    } else {
        logger.Write("solved",
                     result.maxPrInfsbIEq,
                     result.violatedIEq,
                     result.nPrInfsbIEq,
                     result.maxPrInfsbEq,
                     result.violatedEq,
                     result.nPrInfsbEq,
                     result.maxPrInfsbB,
                     result.violatedB,
                     result.nPrInfsbB,
                     result.primalCost,
                     result.dualityGap,
                     result.maxDlInfsb,
                     result.maxNegDl);
    }
    logger.NewLine();

}
void LogGenerator::Err(const std::string& err) {
   logger.Write(err);
   logger.NewLine();
}
/////////////////////
void QualityChecker::Check(const QP_SOLVERS::SolverOutput &output, const QP_SOLVERS::DenseQpProblem &problem) {
    this->output = &output;
    this->problem = &problem;
    assert(output.primal.size() == problem.H.size());
    if (problem.twoSided) {


    } else {
        assert(output.dualA.size() == problem.A.size());
        assert(problem.nEqConstraints <= output.dualA.size());
        nIEqConstraints = static_cast<unsigned int>(problem.A.size())
                - problem.nEqConstraints;
    }
    assert(output.dualL.size() == output.primal.size());
    assert(output.dualU.size() == output.primal.size());
    nVariables = static_cast<unsigned int>(output.primal.size());
    nEqConstraints = problem.nEqConstraints;
    metrics.nEqConstraints = nEqConstraints;
    metrics.nConstraints = problem.A.size();
    metrics.nVariables = nVariables;
    metrics.primalCost = 0.0;
    metrics.time = output.time;
    metrics.nIterations = output.nIterations;
    for (std::size_t i = 0; i < nVariables; ++i) {
        metrics.primalCost += problem.c[i] * output.primal[i] + 0.5 * problem.H[i][i] * output.primal[i] * output.primal[i];
        for (std::size_t j = 0; j < i; ++j) {
            metrics.primalCost += (problem.H[i][j] * output.primal[i] * output.primal[j]);
        }
    }
    CheckPrimalFeasibility();
    CheckDualFeasibility();
    CheckDualityGap();
}

void QualityChecker::CheckPrimalFeasibility() {
    //Ax <= b, l <= x <= u
    //constraints
    double maxInfsblIEq = 0.0;
    double maxInfsblEq = 0.0;
    unsigned int nViolatedIEq = 0;
    unsigned int nViolatedEq = 0;
    double violatedIEq = 0.0;
    double violatedEq = 0.0;
    if (problem->twoSided) {

    } else {
        for (std::size_t i = 0; i < nEqConstraints + nIEqConstraints ; ++i) {
            double constraint = -(problem->b[i]);
            for (std::size_t j = 0; j < nVariables; ++j) {
                constraint += problem->A[i][j] * output->primal[j];
            }
            if (i < nEqConstraints) {
                if (!g_isSame(constraint, 0.0)) {
                    ++nViolatedEq;
                    const double cAbs = std::fabs(constraint);
                    if (cAbs > maxInfsblEq) {
                        violatedEq = problem->b[i];
                        maxInfsblEq = cAbs;
                    }
                }
            } else if (constraint > zeroPlus) {
                ++nViolatedIEq;
                if (constraint > maxInfsblIEq) {
                    violatedIEq = problem->b[i];
                    maxInfsblIEq = constraint;
                }
            }
        }
    }
    //bounds
    double maxInfsblB = 0.0;
    unsigned int nViolatedB = 0;
    double violatedB = 0.0;
    for (std::size_t i = 0; i < nVariables; ++i) {
        assert(!g_isSame(problem->up[i], problem->lw[i]));
        double shiftU = output->primal[i] - problem->up[i];
        double shiftL = problem->lw[i] - output->primal[i];
        //assert(shiftU <= zeroPlus || shiftL <= zeroPlus);
        if (shiftU > zeroPlus) {
            ++nViolatedB;
        } else if (shiftL > zeroPlus) {
            ++nViolatedB;
        }
        double shift = std::fmax(shiftL, shiftU);
        if (shift > maxInfsblB) {
            violatedB = shiftL > shiftU ? problem->lw[i] : problem->up[i];
            maxInfsblB = shift;
        }
    }
    metrics.violatedIEq = violatedIEq;
    metrics.violatedEq = violatedEq;
    metrics.violatedB = violatedB;
    metrics.maxPrInfsbIEq = maxInfsblIEq;
    metrics.maxPrInfsbEq = maxInfsblEq;
    metrics.maxPrInfsbB = maxInfsblB;
    metrics.nPrInfsbIEq = nViolatedIEq;
    metrics.nPrInfsbB = nViolatedB;
    metrics.nPrInfsbEq = nViolatedEq;
}

void QualityChecker::CheckDualFeasibility() {
    // KKT conditions:
    // H * x + c + A_T * lambda - lambdaL + lambdaUp = 0 - dual feasibility
    // lambda, lambdaL, lambdaU >= 0
    // dualityGap = 0 <=> lambda * ( A * x - b) = 0, lambdaL * (x - l) = 0, lambdaU * (x - u) = 0
    // A * x <= b, x <= u, x >= l -- checks in ComputePrInfeasibility
    // dual feasibility
    double maxDualFsb = 0.0;
    unsigned int nViolated = 0;
    std::vector<double> dualFsb = problem->c;
    //H * x + A_T * lambda
    unsigned int nConstraints = nEqConstraints + nIEqConstraints;
    xHx = 0.0;
    for (std::size_t iv = 0; iv < nVariables; ++iv) {
        for (std::size_t j = 0; j < std::max(nVariables, nConstraints); ++j) {
            if (j < nVariables) {
                dualFsb[iv] += problem->H[iv][j] * output->primal[j];
                xHx += output->primal[iv] * problem->H[iv][j] * output->primal[j];
            }
            if (j < nConstraints) {
                dualFsb[iv] += problem->A[j][iv] * output->dualA[j];
            }
        }
    }
    // xTx, cTx, bTl, lTl, uTl are duality gap components
    cTx = 0.0;
    for (std::size_t i = 0; i < nVariables; ++i) {
        cTx += problem->c[i] * output->primal[i];
        dualFsb[i] += output->dualU[i] - output->dualL[i];
        const double violation = std::fabs(dualFsb[i]);
        if (violation > zeroPlus) {
            maxDualFsb = std::fmax(violation, maxDualFsb);
            ++nViolated;
        }
    }
    metrics.maxDlInfsb = maxDualFsb;
    metrics.nDlInfsb = nViolated;
    //complementary slackness
    double maxNegDual = 0.0;
    unsigned int nNegDual = 0;
    bTL = 0.0;
    for (std::size_t i = 0; i < nConstraints; ++i) {
        bTL += problem->b[i] * output->dualA[i];
        if (i >= nEqConstraints && output->dualA[i] < zeroMinus) {
            maxNegDual = std::fmax(maxNegDual, -(output->dualA[i]));
            ++nNegDual;
        }
    }
    lTL = 0.0;
    uTL = 0.0;
    for (std::size_t i = 0; i < nVariables; ++i) {
        lTL += output->dualL[i] * problem->lw[i];
        uTL += output->dualU[i] * problem->up[i];
        if (output->dualL[i] < zeroMinus) {
            maxNegDual= std::fmax(maxNegDual, -(output->dualL[i]));
            ++nNegDual;
        }
        if (output->dualU[i] < zeroMinus) {
            maxNegDual = std::fmax(maxNegDual, -(output->dualU[i]));
            ++nNegDual;
        }
    }
    metrics.maxNegDl = maxNegDual;
    metrics.nNegDl = nViolated;
}
void QualityChecker::CheckDualityGap() {
    metrics.dualityGap  = xHx + cTx + bTL - lTL + uTL;
}
//////////////////////
void QpProblemSolver::SetQpSolver(QP_SOLVERS::QpSolverType solverType) {
    if (solverType == QP_SOLVERS::QpSolverType::NNLS) {
        solver = std::make_unique<QP_SOLVERS::NNLS>();
    } else if (solverType == QP_SOLVERS::QpSolverType::QLD) {
        solver = std::make_unique<QP_SOLVERS::QLD>();
    } else if (solverType == QP_SOLVERS::QpSolverType::DAQP) {
        solver = std::make_unique<QP_SOLVERS::DAQP>();
    }
}
bool QpProblemSolver::SetProblem(const QP_SOLVERS::DenseQpProblem& problem) {
    this->problem = problem;
    if (!solver->SetProblem(problem)) {
        msg = solver->GetInitMsg();
        return false;
    }
    return true;
}
void QpProblemSolver::Solve(QP_SOLVERS::SolverOutput& output){
    solver->Solve(output);
}
void QpProblemSolver::CheckQuality(const QP_SOLVERS::SolverOutput& output,
                                   QP_SOLVERS::QpMetrics& metrics) {
    qlChecker.Check(output, problem);
    metrics = qlChecker.GetMetrics();
}
/////////////////////
bool ProblemCreator::ConfigureFromPath(const std::string& filePath) {
    errMsg.clear();
    TXT_QP_PARSER::DenseProblemFormatter reader;
    const QP_NNLS::DenseQPProblem pr = reader.PrepareProblem(filePath);
    TXT_QP_PARSER::ProblemFormatterStatus retStatus = reader.GetRetStatus();
    if (retStatus.status == false) {
        errMsg = retStatus.errMsg;
        return false;
    } else {
        problem.nEqConstraints = pr.nEqConstraints;
        problem.A = pr.A;
        problem.H = pr.H;
        problem.c = pr.c;
        problem.b = pr.b;
        problem.lw = pr.lw;
        problem.up = pr.up;
        problem.twoSided = false;
        return true;
    }
}
const QP_SOLVERS::DenseQpProblem& ProblemCreator::GetProblem() {
    return problem;
}
/////////////////////
void TestQpSolver(QP_SOLVERS::QpSolverType qpType, unsigned char cType, unsigned char num) {

    std::vector<std::string> filesIEq = {"HS21" ,"HS35", "HS76", "HS268", "HS118", "QPTEST", "S268", "ZECEVIC2",  //small
                                      "KSIP", "QISRAEL", "PRIMALC1", "PRIMALC2", "PRIMALC5", "PRIMALC8", "PRIMAL1", // med
                                      "PRIMAL2", "PRIMAL3", "MOSARQP2" }; //large
    std::vector<std::string> filesEq = {"HS35MOD", "TAME", "HS51", "HS52", "HS53", "GENHS28", "LOTSCHD", "DUALC1", "DUALC2", "QAFIRO", "DUALC5", "DUALC8", //small
                                        "DUAL1", "DUAL2", "DUAL4", "QPCBLEND", "QSHARE2B", "CVXQP2_S", "QADLITTL", "DUAL3", "CVXQP1_S", "CVXQP3_S", "QSCAGR7",
                                        "DPKLO1", "QPCBOEI2", "QRECIPE", "VALUES",  //medium
                                        "QSC205", "QSHARE1B", "QBRANDY", "QBEACONF", "QE226", "QGROW7", "QBORE3D", "QSCORPIO", "QCAPRI", "QFORPLAN", "QPCBOEI1",
                                        "QSCFXM1", "QBANDM", "QPCSTAIR", "QSTAIR", "QSCTAP1", "QSCAGR25", //large
                                        "QGROW15", "QSCSD1", "GOULDQP3", "GOULDQP2", "QETAMACR", "QFFFFF80", "QGROW22", "CVXQP2_M", "QSCFXM2", "CVXQP1_M", "CVXQP3_M" }; // > 10MB

    //std::vector<std::string> filesEq = {"GENHS28"};
    const std::string logFile = TST_CONFIG::LOG_DIR + "report_" + std::to_string(num) + ".txt";
    QpProblemSolver qpSolverWrapper;
    qpSolverWrapper.SetQpSolver(qpType);
    ProblemCreator pCreator;
    QP_SOLVERS::SolverOutput output;
    QP_SOLVERS::QpMetrics metrics;
    LogGenerator logger;
    logger.SetLogFile(logFile);
    const std::vector<std::string>& files = cType == 1 ? filesIEq : filesEq;
    const bool noEqC = cType == 1 ? true : false;
    for (const auto& file: files) {
        const std::string& txtQpRoot = noEqC ?  TST_CONFIG::DENSE_ONLY_INQ_CONSTR_PATH
                                              : TST_CONFIG::DENSE_INQ_AND_EQ_CONSTR_PATH;
        if (pCreator.ConfigureFromPath(txtQpRoot + file + ".txt")) {
            if (!qpSolverWrapper.SetProblem(pCreator.GetProblem())) {
                logger.Err(file + " " + pCreator.errMsg);
            } else {
                qpSolverWrapper.Solve(output);
                bool solved = (output.status == QP_SOLVERS::SolverRetStatus::SUCCESS) ||
                          (output.status == QP_SOLVERS::SolverRetStatus::ITERATIONS);
                if (!solved) {
                    logger.Err(output.msg + file);
                } else {
                    qpSolverWrapper.CheckQuality(output, metrics);
                    logger.Dump(file, metrics);
                }
            }
        } else {
            logger.Err("configuration error: " + file);
        }
    }
}

} //namespace SOLVER_COMPARATOR
