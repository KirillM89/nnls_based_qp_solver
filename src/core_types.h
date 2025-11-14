#ifndef CORE_TYPES_H
#define CORE_TYPES_H
#include <list>
#include <string>
#include <deque>
#include <unordered_set>
#include <set>
#include "types.h"
namespace QP_NNLS {
enum class LinSolverType {
    CUMULATIVE_LDLT = 0,
    //CUMULATIVE_EG_LDLT,
    DYNAMIC_LDLT,
    MSS1,
};

enum class InitStageStatus {
    SUCCESS = 0,
    D_Z = 1,
    D_N = 2,
    D_ZN = 3,
};

enum class PrimalLoopExitStatus {
    EMPTY_ACTIVE_SET = 0,
    ALL_PRIMAL_POSITIVE,
    ITERATIONS,
    EMPTY_ACTIVE_SET_ON_ZERO_ITERATION,
    SINGULAR_MATRIX,
    DIDNT_STARTED,
    LINE_SEARCH_FAILED,
    ZERO_STEP,
    UNKNOWN
};

enum class DualLoopExitStatus {
    ALL_DUAL_POSITIVE = 0,
    FULL_ACTIVE_SET = 1,
    ITERATIONS = 2,
    INFEASIBILITY = 3,
    UNKNOWN
};


struct LinSolverOutput {
    bool emptyInput = false;
    unsg_t nDNegative = 0; // number of d<=0 in LDLT
    std::vector<double> solution;
    std::list<unsg_t> indices;
};

struct ActiveSetUpdateSettings {
    int rptInterval = 0;
    bool rejectSingular = false;
};

struct CoreSettings : public ActiveSetUpdateSettings  {
    LinSolverType linSolverType = LinSolverType::CUMULATIVE_LDLT;
    unsg_t nDualIterations = 1000;
    unsg_t nPrimalIterations = 100;
    unsg_t logLevel = 0u;
    double nnlsResidNormFsb = 1.0e-12; // infeasibility criterion
    double origPrimalFsb = 1.0e-6;     // Ax - b <= origPrimalFsb
    double nnlsPrimalZero = 1.0e-30;   // x==0 if |x| < nnlsPrimalZero
    double prLtZero = 1.0e-30;         // x < 0 if x < prLtZero
    bool gammaUpdate = true;
    bool largeBoundsPenalty = false;
    bool zeroDCor = true;
    bool checkFactorization = false;
};

struct Settings : public CoreSettings {
//RESERVED
};

struct DenseQPProblem {
    // 0.5 * x_T * H * x + c_T * x
    // Ax <= b; lw <= x <= up
    // First nEqConstraints in A are equality constraints
    matrix_t H;
    matrix_t A;
    std::vector<double> b;
    std::vector<double> c;
    std::vector<double> up;
    std::vector<double> lw;
    unsg_t nEqConstraints = 0;
};

struct LinSolverTime {
    unsigned long long us;
    unsg_t nConstraints;
};

struct SolverOutput {
    DualLoopExitStatus dualExitStatus;
    PrimalLoopExitStatus primalExitStatus;
    unsg_t mDef;
    unsg_t nDualIterations;
    unsg_t nVariables;
    unsg_t nConstraints;
    unsg_t nEqConstraints;
    double maxViolation;
    double dualityGap;
    double cost;
    std::vector<double> x;
    std::vector<double> lambda;
    std::vector<double> lambdaLw;
    std::vector<double> lambdaUp;
    std::vector<double> violations;
};
} //

#endif // CORE_TYPES_H
