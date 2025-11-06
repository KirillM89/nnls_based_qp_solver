
#include <cmath>
#include <algorithm>
#include <iostream>
#include "core.h"
#define BTOL 1.0e10

namespace QP_NNLS {
Core::Core():
    timer(std::make_unique<wcTimer>()) //wall-clock timer
{}

void Core::Init(const Input& problem) {
    nPVariables = static_cast<unsg_t>(problem.H.size());
    nPConstraints = static_cast<unsg_t>(problem.A.size());
    nEqConstraints = problem.nEqConstraints;
    nVariables = nPVariables;
    nConstraints = nPConstraints + 2 * nPVariables;
    if (config.largeBoundsPenalty) {
        nPVariables += 1;
        nConstraints += 1;
    }
    linSolverTimes.clear();
    mEps = g_GetMachineEps();
    minEl = 10.0 * sqrt(mEps);
    scaleFactorDB = 1.0;
    initStatus = InitStageStatus::SUCCESS;
    nDualIterations = static_cast<unsg_t>((2 * nPVariables + nPConstraints) * 1.5);
    gamma = 1.0;
    singularIndices.clear();
}
void Core::Set(const Configuration& config) {
    this->config = config;
}
void Core::SetCallback(Callback* callback) {
    uCallback = callback;
}
bool Core::InitProblem(const Input &problem) {
    if (!PrepareNNLS(problem)) {
        return false;
    }
    if (uCallback) {
        uCallback->SetLogLevel(config.logLevel);
        uCallback->Init();
    }
    linSolverTimes.reserve(nDualIterations);
    SetInitData(problem);
    return true;
}
void Core::CheckFactorization(const matrix_t& ld, const matrix_t& H, const std::vector<int>& pmt) {
    //ld : L - is in low diagonal part, D - is on diagonal
    //H = P * L * D * L_T * P_T
    const std::size_t nV = H.size();
    matrix_t ldt(ld);
    for (std::size_t i = 0; i < nV; ++i) {
        ldt[i][i] = 1.0;
        for (std::size_t j = i + 1; j < nV; ++j) {
            ldt[i][j] = 0.0;
        }
    }
    //apply permutations in reverse order
    for (int v = nV - 1; v >= 0; --v) {
        if (pmt[v] != -1) {
            std::swap(ldt[v], ldt[pmt[v]]);
        }
    }
    matrix_t lByD(ldt);
    for (std::size_t i = 0; i < nV; ++i) {
        for (std::size_t j = 0; j < nV; ++j) {
            lByD[i][j] *= ld[j][j];
        }
    }
    matrix_t HF(H);
    M1M2T(lByD, ldt, HF);
    const double errTol = 1.0e-7;
    for (std::size_t i = 0; i < nV; ++i) {
        for (std::size_t j = 0; j < nV; ++j) {
            const double diff = H[i][j] - HF[i][j];
            if (std::fabs(diff) >= errTol) {
                std::cout << "FACTORIZATION ERROR: " << diff << std::endl;
            }
        }
    }
}
matrix_t Core::ComputeLDLT(const matrix_t& H) {
    matrix_t ld(H);
    const std::size_t nV = H.size();
    std::vector<int> pmt(nV, -1.0); // must be filled with -1.0 by default
    std::size_t nZr, nNg; //number of zero and negative diagonal elements in matrix D of LDL_T decomposition
    const bool inPlaceDZeroCor = false;
    unsigned char rStat  = InPlaceLdlt(ld, pmt, nZr, nNg, true, inPlaceDZeroCor);
    //continue if D has only positive values or zero values with option zeroDCor==true
    if ((rStat == 1u && !config.posDefCorrection) || rStat > 1u) {
        if (rStat == 1u) {
            initStatus = InitStageStatus::D_Z;
        } else if (rStat == 2u) {
            initStatus = InitStageStatus::D_N;
        } else if (rStat == 3u) {
            initStatus = InitStageStatus::D_ZN;
        }
    } else {
       // if (settings.checkFactorization) {
       //     CheckFactorization(ld, H, pmt);
       // }
        output.isPositiveDefinite = false;
        const double minChol = 2.0 * mEps; // min diagonal value
        const double maxCondNumInv = sqrt(mEps);
        const double maxEgVal = ld[0][0];
        const double minEgVal = maxCondNumInv * maxEgVal;
        for (unsg_t v = 0; v < nV; ++v) {
            if (config.posDefCorrection && ld[v][v] < minEgVal) {  // replace all values which violate condition number
                ld[v][v] = std::fmax(minChol, maxCondNumInv *  minEgVal);  // (CondNum / maxCondNum) * minEgVal = (maxEgVal / (minEgVal * maxCondNum)) * minEgVal = maxEgVal / maxCondNum
                output.isPositiveDefinite = false;
                ws.dCorrected[v] = true;
            }
            ws.Chol[v] = 1.0 / sqrt(ld[v][v]);
        }
        InvertL(ld); // L_-1 will be saved in upper diagonal part
                     // lower diagonal part will be filled with zeros
                     // diagonal will be filled with 1.0
        for (int v = nV - 1; v >= 0; --v) {
            if (pmt[v] != -1) {
                std::swap(ld[v], ld[pmt[v]]); //P * L_-T
            }
        }
    }
    return ld;
}
void Core::ComputeLS4Constraints(const matrix_t& ld,  const Input& problem) {
    const std::size_t nV = config.largeBoundsPenalty ? nVariables - 1 : nVariables;
    const std::size_t nC = config.largeBoundsPenalty ? nConstraints -  2 * nV - 1 : nConstraints - 2 * nV;
    for (std::size_t c = 0; c < nC; ++c) {
        double b = problem.b[c];
        int stat = 1;
        bool scale = false;
        if (config.largeBoundsPenalty) {
            AddExtraComponent(c, b, stat);
        } else {
            if (b < -BTOL) {
                b = -BTOL;
            } else if (b > BTOL) {
                b = BTOL;
            } else {
                scale = true;
            }
            ws.s[c] = b;
        }
        double norm2 = ws.M[c].back() * ws.M[c].back();
        for (std::size_t v = 0; v < nV; ++v) {
            const std::size_t iBnd = nC + 2 * v;
            if (c == 0) {
                const std::vector<double>& minLd = -ld[v];
                std::copy(ld[v].begin(), ld[v].end(), ws.M[iBnd].begin());   // M part corresponding to bounds
                std::copy(minLd.begin(), minLd.end(), ws.M[iBnd + 1].begin());
            }
            for (std::size_t j = 0; j < nV; ++j) {
                if (c == 0) {
                    ws.v[v] += ld[j][v] * problem.c[j] * ws.Chol[v]; // v = Chol_-T * c_new =  Chol_-T * (P * L_-T)_T * c
                }
                ws.M[c][v] += problem.A[c][j] * ld[j][v]; // A_new = A * (P * L_-T)
            }
            ws.M[c][v] *= ws.Chol[v];
            ws.s[c] += ws.M[c][v] * ws.v[v];
            norm2 += ws.M[c][v] * ws.M[c][v];
        }
        ws.lambda[c] = norm2;
        UpdateScaleFactor(c);
    }
}
void Core::ComputeLS4Bounds(const Input& problem) {
    const std::size_t nV = config.largeBoundsPenalty ? nVariables - 1 : nVariables;
    const std::size_t nC = config.largeBoundsPenalty ? nConstraints -  2 * nV - 1 : nConstraints - 2 * nV;
    for (unsg_t c = 0; c < nV; ++c) {
        bool scaleU = false;
        bool scaleL = false;
        double norm2 = 0.0;
        double bUp = problem.up[c];
        double bLw = -problem.lw[c];
        const std::size_t iBnd = nC + 2 * c;
        int statUp = 1, statLw = 1;
        if (config.largeBoundsPenalty) {
            AddExtraComponent(iBnd, bUp, statUp);
            AddExtraComponent(iBnd + 1, bLw, statLw);
        } else {
            if (bUp > BTOL) {
                bUp = BTOL;
            } else if (bUp < -BTOL) {
                bUp = -BTOL;
            } else {
                scaleU = true;
            }
            if (bLw > BTOL) {
                bLw = BTOL;
            } else if (bLw < -BTOL) {
                bLw = -BTOL;
            } else {
                scaleL = true;
            }
            ws.s[iBnd] = bUp;
            ws.s[iBnd  + 1] = bLw;
        }
        for (unsg_t v = 0; v < nVariables; ++v) {
            ws.M[iBnd][v] *= ws.Chol[v];
            ws.M[iBnd + 1][v] *= ws.Chol[v];
            ws.s[iBnd] += ws.M[iBnd][v] * ws.v[v];
            ws.s[iBnd + 1] += ws.M[iBnd  + 1][v] * ws.v[v];
            norm2 += ws.M[iBnd][v] * ws.M[iBnd][v];
        }
        ws.lambda[iBnd] = norm2;
        UpdateScaleFactor(iBnd);
        ws.lambda[iBnd + 1] = norm2;
        UpdateScaleFactor(iBnd + 1);
    }
}

void Core::Scale() { 
    if (isSame(scaleFactorDB, 1.0)) {
        //scaleFactorDB = 10. * mEps;
    }
    scaleFactorDB = sqrt(scaleFactorDB);
    config.origPrimalFsb *= scaleFactorDB;
    const double minS = g_GetMachineEps();
    //s scaling and ortogonalization
    for (std::size_t r = 0; r < nConstraints; ++r) {
        if (r < nVariables) {
            ws.v[r] *= scaleFactorDB;
        }
        ws.s[r] *= scaleFactorDB;
        double fullNorm = ws.lambda[r] + ws.s[r] * ws.s[r];
        if (!isSame(fullNorm, 0.0)) {
            double scaleCoef = 1.0 / sqrt(fullNorm);
            ws.lambda[r] = scaleCoef;
            ws.s[r] *= scaleCoef;
            double maxM = std::numeric_limits<double>::min();
            unsg_t iMaxM = 0;
            for (std::size_t c = 0; c < nVariables; ++c) {
                ws.M[r][c] *= scaleCoef;
                if (ws.M[r][c] > maxM) {
                    maxM = ws.M[r][c];
                    iMaxM = c;
                }
            }
            ws.M[r][iMaxM] += minS * ws.M[r][iMaxM];
        }       
    }
}
void Core::AllocateWs() {
    ws.primal = std::vector<double>(nConstraints, 0.0);     // primal vars for nnls problem
    ws.dual = std::vector<double>(nConstraints, 0.0);       // dual vars for nnls problem
    ws.zp = std::vector<double>(nConstraints, 0.0);         // aux array for compuatation of primal components on current active set
    ws.lambda = std::vector<double>(nConstraints, 1.0);     // dual vars for orig problem
    ws.s = std::vector<double>(nConstraints, 0.0);          // vector s = b + M * v
    ws.aux = std::vector<double>(nConstraints, 0.0);        // auxilary array
    ws.x = std::vector<double>(nVariables, 0.0);            // primal vars for orig problem
    ws.v = std::vector<double>(nVariables, 0.0);            // vector v = L_-1 * c
    ws.Chol = std::vector<double>(nVariables, 1.0);         // invert Choletsky factor
    if (config.posDefCorrection) {
        ws.dCorrected = std::vector<bool>(nVariables, false);
    }
    ws.M = matrix_t(nConstraints, std::vector<double>(nVariables, 0.0)); // M = A * L_-1
    ws.activeConstraints.clear();
    ws.linEqConstraints.clear();
    ws.addHistory.clear();
    for (unsg_t i = 0; i < nEqConstraints; ++i) {
        ws.linEqConstraints.insert(i);
    }
}

bool Core::PrepareNNLS(const Input &problem) {
    Init(problem);
    AllocateWs();
    timer->Start();
    const matrix_t ld = ComputeLDLT(problem.H);
    if (initStatus != InitStageStatus::SUCCESS) {
        return false;
    }
    ComputeLS4Constraints(ld, problem);
    ComputeLS4Bounds(problem);
    if (config.largeBoundsPenalty) {
        ws.M.back().back() = -1.0;
    }
    Scale();
    if (uCallback) {
        TimeInterval(uCallback->initData.tM);
    }
    SetLinearSolver();
    for (auto indx : ws.linEqConstraints) {
        AddToActiveSet(indx);
    }
    if (!ws.activeConstraints.empty()) {
        SolvePrimal();
        ws.primal =  ws.zp;
    }
    return true;
}
void Core::UpdateScaleFactor(std::size_t iConstraint) {
    //compute critical scale factor for s. This need to
    //prevent the appearance of too small [M s] matrix components during ortogonalization
    //real scale factor must be >= critical
    //TODO: prevent appearance of too big components !
    if (!isSame(ws.s[iConstraint], 0.0)) {
        const double sp2 = ws.s[iConstraint] * ws.s[iConstraint];
        const double norm2 = ws.lambda[iConstraint];
        for (std::size_t i = 0; i < nVariables; ++i) {
            const double m = ws.M[iConstraint][i];
            const double diff = m * m - minEl * norm2;
            //not positive diff means that it's impossible to set current m to value >= minEl
            if (isSame(m, 0.0) || diff <= 0.0) {
                continue;
            }
            double rat = diff / (minEl * sp2);
            if (rat < minEl) {
                scaleFactorDB = 0.5 * mEps; // mEps because real scale factor is sqrt(scaleFactorDB)
            } else {
                scaleFactorDB = std::fmin(scaleFactorDB, rat); // minEl > mEps
            }
        }
    }
}
void Core::AddExtraComponent(unsg_t indx, double bound, int& stat) {
    if (bound >= BTOL) {
        ws.M[indx].back() = -1.0;
        ws.s[indx] = 0.0;
        stat = 1;
    } else if (bound <= -BTOL) {
        ws.M[indx].back() = 1.0;
        ws.s[indx] = 0.0;
        stat = -1;
    } else {
        ws.M[indx].back() = 0.0;
        ws.s[indx] = bound;
        stat = 0;
    }
}
void Core::TimeInterval(std::string& buf) {
    TimeIntervals tIntervals;
    timer->toIntervals(timer->Ticks(), tIntervals);
    buf = std::to_string(tIntervals.minutes) + " min " +
          std::to_string(tIntervals.sec) + " sec " +
          std::to_string(tIntervals.ms) + " ms " +
          std::to_string(tIntervals.mus) + " mus";
}
bool Core::OrigInfeasible() {
    MultTransp(ws.M, ws.primal, ws.activeConstraints, ws.x); // M_T * primal
    styGamma = gamma + DotProduct(ws.s, ws.primal, ws.activeConstraints);
    rsNorm = DotProduct(ws.x, ws.x) + styGamma * styGamma;
    return rsNorm < config.nnlsResidNormFsb;
}

bool Core::FullActiveSet() {
    return (static_cast<unsg_t>(ws.activeConstraints.size()) == nConstraints);
}
void Core::ComputeDualVariable() {
    Mult(ws.M, ws.x, ws.dual); // M * M_T * primal,  M_T * primal is saved in ws.x
    for (unsg_t i = 0 ; i < nConstraints; ++i) {
        ws.dual[i] += styGamma * ws.s[i];
    }
}
bool Core::SkipCandidate(unsg_t indx) {
    if (config.rejectSingular && singularIndices.find(indx) != singularIndices.end()) {
        return true;
    } else {
        return false;
    }
}
void Core::AddToActiveSet(unsg_t indx) {
    if (indx != nConstraints) {
        ws.activeConstraints.insert(indx);
        ws.addHistory.push_back(indx);
        lSolver->Add(indx);
    }
}
void Core::RmvFromActiveSet(unsg_t indx) {
    if (ws.linEqConstraints.find(indx) == ws.linEqConstraints.end()) {
        ws.activeConstraints.erase(indx);
        lSolver->Delete(indx);
    }
}
bool Core::IsCandidateForNewActive(unsg_t indx, double toCompare, bool skip) {
    bool res = false;
    const double dl = ws.dual[indx];
    if (ws.activeConstraints.find(indx) != ws.activeConstraints.end() || SkipCandidate(indx)) {
        res = false;
    } else if ((dl < dualTolerance && dl < toCompare)) {
        newActiveIndex = indx;
        res = true;
    }
    return res;
}
unsg_t Core::SelectNewActiveComponent() {
    //Algorith guarantees that all dual on active set must be nonegative
    //but they can be negative due to numerical errors. Such components will be set to zero
    double newActive = std::numeric_limits<double>::max();
    newActiveIndex = nConstraints; //default value
    //check inactive components
    for (unsg_t i = nEqConstraints; i < nConstraints; ++i) {
        if (IsCandidateForNewActive(i, newActive)) {
            newActive = ws.dual[i];
        }
    }
    return newActiveIndex;
}

unsg_t Core::SolvePrimal() {
    timer->Ticks();
    const std::size_t nActive = ws.activeConstraints.size();
    lSolver->Set(config.rejectSingular);
    lSolver->SetGamma(gamma);
    const LinSolverOutput& lsOutput = lSolver->Solve();
    assert(lsOutput.indices.size() == ws.activeConstraints.size());
    for (auto indx : lsOutput.indices) {
        assert(ws.activeConstraints.find(indx) != ws.activeConstraints.end());
    }
    std::fill(ws.zp.begin(), ws.zp.end(), 0.0);
    if (lsOutput.nDNegative == 0 || !config.rejectSingular) {
        ws.negativeZp.clear();
        std::size_t i = 0;
        if (nActive > 0) {
            for (auto indx: lsOutput.indices) {
                ws.zp[indx] = lsOutput.solution[i];
                if (indx >= nEqConstraints && ws.zp[indx] <= 1.0e-30) {
                    ws.negativeZp.insert(indx);
                }
                ++i;
            }
        }
    }
    linSolverTimes.emplace_back();
    linSolverTimes.back().us = timer->Ticks();
    linSolverTimes.back().nConstraints = nActive;
    // TODO: check quality
    return lsOutput.nDNegative;
}

bool Core::MakeLineSearch() {
    double minStep = 1.0;
    bool nonZeroStep = false;
    std::size_t iBlocking = nConstraints;
    std::size_t nBlocking = 0;
    const double strictZeroTol = 1.0e-30;
    // step = 0: primal_i = 0 for any zp_i
    // but all the primal except primal for new active set index must be strict positive
    // so step can be zero only if zp for new active component is nonpositive it's blocking component
    for (auto indx: ws.negativeZp) {
        const double primal = ws.primal[indx];
        std::cout << "dit: " << dualIteration << " pit:" << primalIteration << " idx:" << indx << " x:" << primal << " y:" << ws.zp[indx] << std::endl;
        const double denominator = primal - ws.zp[indx];
        if (!isSame(denominator, 0.0, strictZeroTol)) { // zp < 0
            double rat = primal / denominator;
            assert(rat >= 0.0);
            if (!isSame(primal, 0.0, strictZeroTol)) {
                minStep = std::fmin(minStep, rat);
                nonZeroStep = true;
            } else {
                iBlocking = indx;
                ++nBlocking;
            }
        } else { // zp == 0
            assert(isSame(primal, 0.0, strictZeroTol));
            iBlocking = indx;
            ++nBlocking;
        }
    }
    assert(nBlocking <= 1);
    std::cout << "step " << minStep << std::endl;
    if (nonZeroStep) {
        // At least for one component LS step was found
        gammaCorrection = 0.0;
        for (unsg_t i = 0; i < nConstraints; ++i) {
            if (i != iBlocking) {
                //don't shift blocking component
                ws.primal[i] += minStep * (ws.zp[i] - ws.primal[i]);
                if (i >= nEqConstraints) {
                    ws.primal[i] = std::fmax(ws.primal[i], 0.0); // correct possible numerical errors
                }
            }
            //remove constraints corresponding to zero or neg primal from active set
            if (ws.primal[i] <= 1.0e-30) {
                if (i >= nEqConstraints) {
                    if (ws.activeConstraints.find(i) != ws.activeConstraints.end()) {
                        gammaCorrection += std::fabs(ws.s[i]);
                    }
                    RmvFromActiveSet(i);
                }
            }
        }
        UpdateGammaOnPrimalIteration();
    } else {
        //cycling occurs: example case QFFFFF80
        RmvFromActiveSet(iBlocking);
    }
    return true;
}
int Core::UpdatePrimal() {
    int res = 0;
    if (SolvePrimal() > 0) {
        res |= SINGULARITY;
    }
    if (res == 0 || !config.rejectSingular) {
        if (!ws.negativeZp.empty()) {
            singularIndices.clear();
            if (!MakeLineSearch()) {
                res |= LINE_SEARCH_FAILED;
            }
        } else {
            std::copy(ws.zp.begin(), ws.zp.end(), ws.primal.begin());
        }
    }
    return res;
}

void Core::UnscaleD() {
    if (scaleFactorDB < 1.0) {
        config.origPrimalFsb /= scaleFactorDB;
    }
}
void Core::UpdateGammaOnPrimalIteration() {
    if (config.gammaUpdate == true) {
        gamma = std::fabs(gamma - gammaCorrection);
    }
}
void Core::UpdateGammaOnDualIteration() {
    if (config.gammaUpdate == true) {
        gamma += std::fabs(ws.s[newActiveIndex]);
    }
}

void Core::ComputeOrigSolution() {
    const double sty = DotProduct(ws.s, ws.primal, ws.activeConstraints);
    const double lambdaTerm = -1.0 / (gamma + sty);
    std::copy(ws.lambda.begin(), ws.lambda.end(), ws.aux.begin()); //save [M s] norms
    for (unsg_t i = 0; i < nConstraints; ++i) {
        ws.lambda[i] = lambdaTerm * ws.primal[i];
    }
    MultTransp(ws.M, ws.lambda, ws.activeConstraints, ws.x);
    const double invScaleFactor = 1.0 / scaleFactorDB;
    for (unsg_t i = 0; i < nVariables; ++i) {
        ws.x[i] = (ws.x[i] - ws.v[i]) * ws.Chol[i];
    }
    //compute cost
    //f = 0.5 * x_T * H * x + c_T * x =
    //  = 0.5 * x_T * Ch * Ch * x + v_T * Ch * x
    cost = 0.0;
    for (unsg_t i = 0; i < nPVariables; ++i) {
        double cx = ws.x[i] / ws.Chol[i];
        cost += (0.5 * cx + ws.v[i]) * cx;
    }
    cost *= invScaleFactor * invScaleFactor;
    for (unsg_t i = 0; i < nConstraints; ++i) {
        ws.lambda[i] *= (-ws.aux[i] * invScaleFactor);
    }
    //recompute x
    for (std::size_t i = 0; i < nPVariables; ++i) {
        ws.aux[i] = 0.0;
        for (std::size_t j = 0; j < nPVariables; ++j) {
            ws.aux[i] += (ws.M[nPConstraints + 2 * i][j] * ws.x[j] / ws.Chol[j]);
        }
        ws.aux[i] /= (ws.aux[nPConstraints + 2 * i] * scaleFactorDB);
    }
    std::copy(ws.aux.begin(), ws.aux.begin() + nPVariables, ws.x.begin());
}


void Core::FillOutput() {
    output.dualExitStatus = static_cast<unsigned char>(dualExitStatus);
    output.primalExitStatus = static_cast<unsigned char>(primalExitStatus);
    if (dualExitStatus != DualLoopExitStatus::INFEASIBILITY){
        output.x = std::vector<double>(ws.x.begin(), ws.x.begin() + nPVariables);
        output.lambdaC.resize(nPConstraints, 0.0);
        output.lambdaLw.resize(nPVariables, 0.0);
        output.lambdaUp.resize(nPVariables, 0.0);
        for (std::size_t i = 0; i < nPConstraints; ++i) {
            output.lambdaC[i] = ws.lambda[i];
        }
        for (std::size_t i = 0; i < nPVariables; ++i) {
            output.lambdaUp[i] = ws.lambda[nPConstraints + 2 * i];
            output.lambdaLw[i] = ws.lambda[nPConstraints + 2 * i + 1];
        }
        output.cost = cost;
    }
    output.nIterations = dualIteration;
}
void Core::SetInitData(const Input &problem) {
    if (config.logLevel >= 1u) {
        uCallback->initData.nVariables= nPVariables;
        uCallback->initData.nConstraints = nPConstraints;
        uCallback->initData.nEqConstraints = nEqConstraints;
        if (config.logLevel >= 2u) {
            uCallback->initData.scaleDB = scaleFactorDB;
            if (config.logLevel >= 3u) {
                uCallback->initData.Chol = &ws.Chol;
                uCallback->initData.CholInv = &ws.Chol;
                uCallback->initData.M = &ws.M;
                uCallback->initData.s = &ws.s;
                uCallback->initData.c = &const_cast<std::vector<double>&>(problem.c);
                uCallback->initData.b = &const_cast<std::vector<double>&>(problem.b);
            }
        }
        uCallback -> ProcessData(1);
    }
}
void Core::SetIterationData() {
    if (config.logLevel >= 2u) {
        uCallback->iterData.iteration = dualIteration;
        uCallback->iterData.newIndex = newActiveIndex;
        uCallback->iterData.gamma = gamma;
        uCallback->iterData.dualTol = dualTolerance;
        uCallback->iterData.rsNorm = rsNorm;
        if (config.logLevel >= 3u) {
            uCallback->iterData.activeSet = &ws.activeConstraints;
            uCallback->iterData.primal = &ws.primal;
            uCallback->iterData.dual = &ws.dual;
            uCallback->iterData.zp = &ws.zp;
        }
        uCallback->ProcessData(2);
    }
}

void Core::SetFinalData() {
    if (config.logLevel >= 1u) {
        uCallback->finalData.dualStatus = static_cast<unsigned char>(dualExitStatus);
        uCallback->finalData.primalStatus = static_cast<unsigned char>(primalExitStatus);
        uCallback->finalData.nIterations = nDualIterations;
        if (dualExitStatus != DualLoopExitStatus::INFEASIBILITY) {
            uCallback->finalData.cost = output.cost;
            uCallback->finalData.x = &output.x;
            uCallback->finalData.lambda = &output.lambdaC;
            uCallback->finalData.lambdaLw = &output.lambdaLw;
            uCallback->finalData.lambdaUp = &output.lambdaUp;
        }
        uCallback->ProcessData(3);
    }
}
void Core::SetLinearSolver() {
    if (config.linSolverType == 2) {
   //     lSolver = std::make_unique<CumulativeEGNSolver>(ws.M, ws.s);
    } else if (config.linSolverType == 1) {
   //     lSolver = std::make_unique<MssCumulativeSolver>(ws.M, ws.s);
    } else {
        lSolver = std::make_unique<CumulativeLDLTSolver>(ws.M, ws.s, config.rejectSingular);
    }
}
void Core::Solve() {
    dualExitStatus = DualLoopExitStatus::UNKNOWN;
    primalExitStatus = PrimalLoopExitStatus::DIDNT_STARTED;
    rsNorm = std::numeric_limits<double>::max();
    cost = std::numeric_limits<double>::max();
    gamma = 1.0;
    newActiveIndex = nConstraints;
    dualIteration = 0;
    while (dualIteration < nDualIterations) {
        if (OrigInfeasible()) {
            dualExitStatus = DualLoopExitStatus::INFEASIBILITY;
            break;
        }
        if (FullActiveSet()) {
            dualExitStatus = DualLoopExitStatus::FULL_ACTIVE_SET;
            break;
        }
        ComputeDualVariable();
        dualTolerance = -styGamma * config.origPrimalFsb; // primal feasiblility was scaled in DB scaling
        SelectNewActiveComponent();
        if(newActiveIndex == nConstraints) { //set to nConstraints in not found
            //Proccess singular components after all nonsingular
            if (config.rejectSingular) {
                config.rejectSingular = false;
                singularIndices.clear();
                continue;
            }
            dualExitStatus = DualLoopExitStatus::ALL_DUAL_POSITIVE;
            break;
        }
        UpdateGammaOnDualIteration();
        AddToActiveSet(newActiveIndex);
        primalIteration = 0;
        primalExitStatus = PrimalLoopExitStatus::UNKNOWN;
        while (primalIteration < nPrimalIterations) {
            if (ws.activeConstraints.empty()) {
                primalExitStatus = primalIteration == 0 ? PrimalLoopExitStatus::EMPTY_ACTIVE_SET_ON_ZERO_ITERATION :
                                                          PrimalLoopExitStatus::EMPTY_ACTIVE_SET;
                break;
            }
            int prStat = UpdatePrimal();
            const bool success = (prStat == 0) ||((prStat == SINGULARITY)
                    && !config.rejectSingular);
            const bool rejectSingular = (prStat == SINGULARITY) && config.rejectSingular;
            if (success) {
                if (ws.negativeZp.empty()) {
                    primalExitStatus = PrimalLoopExitStatus::ALL_PRIMAL_POSITIVE;
                    break;
                } else {
                    ++primalIteration;
                }
            } else if (rejectSingular) {
                primalExitStatus = PrimalLoopExitStatus::SINGULAR_MATRIX;
                RmvFromActiveSet(newActiveIndex);
                singularIndices.insert(newActiveIndex);
                break;
            } else {
                primalExitStatus = PrimalLoopExitStatus::LINE_SEARCH_FAILED;
                break;
            }
        }
        if (primalIteration >= nPrimalIterations) {
            primalExitStatus = PrimalLoopExitStatus::ITERATIONS;
        }
        SetIterationData();
        ++dualIteration;
    }

    if (dualIteration >= nDualIterations) {
        dualExitStatus = DualLoopExitStatus::ITERATIONS;
    }
    if (dualExitStatus == DualLoopExitStatus::ALL_DUAL_POSITIVE ||
        dualExitStatus == DualLoopExitStatus::FULL_ACTIVE_SET ||
        dualExitStatus == DualLoopExitStatus::ITERATIONS) {
        SolvePrimal();
        ws.primal = std::move(ws.zp);
    }
    if (OrigInfeasible()) {
        dualExitStatus = DualLoopExitStatus::INFEASIBILITY;
    }
    if (dualExitStatus != DualLoopExitStatus::INFEASIBILITY) {
        ComputeOrigSolution();
        UnscaleD();
    }
    FillOutput();
    SetFinalData();
}

}
