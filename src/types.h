
#ifndef NNLS_QP_SOLVER_TYPES_H
#define NNLS_QP_SOLVER_TYPES_H
#include <vector>
#include <list>
#include <string>
#include <deque>
#include <unordered_set>
#include <set>
#include "timers.h"
namespace QP_NNLS {
using matrix_t = std::vector<std::vector<double>>;
using unsg_t = unsigned int;

double g_GetMachineEps();

enum class LinSolverType {
    CUMULATIVE_LDLT = 0,
    CUMULATIVE_EG_LDLT,
    DYNAMIC_LDLT,
    MSS1,
};

enum class InitStageStatus {
    SUCCESS = 0,
    D_Z,
    D_N,
    D_ZN,
};

struct LinSolverOutput {
    bool emptyInput = false;
    unsg_t nDNegative = std::numeric_limits<unsg_t>::max(); // number of d<=0 in LDLT
    std::vector<double> solution;
    std::list<unsg_t>  indices;
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

enum class DualLoopExitStatus {
    ALL_DUAL_POSITIVE = 0,
	FULL_ACTIVE_SET,
	ITERATIONS,
	INFEASIBILITY,
	UNKNOWN
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

struct LinSolverTime {
    ticks_t us;
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

} // namespace QP_NNLS
#endif
