
#ifndef NNLS_QP_SOLVER_TYPES_H
#define NNLS_QP_SOLVER_TYPES_H

#include <vector>
#define FP_PREC 1
#if FP_PREC == 0
#define fp_t float
#elif FP_PREC == 1
#define fp_t double
#endif
namespace QP_NNLS {
static constexpr unsigned char LIN_SOLVER_TYPE_DEFAULT = 0;
static constexpr unsigned char LOG_LEVEL_DEFAULT = 0;
static constexpr bool USE_LARGE_BOUNDS_PENALTY_DEFAULT = false;
static constexpr bool USE_POS_DEF_CORRECTION_DEFAULT = false;
static constexpr bool USE_GAMMA_UPDATE_DEFAULT = true;
#if FP_PREC == 1
static constexpr fp_t NNLS_RESID_TOL_DEFAULT = 1.0e-12;
static constexpr fp_t PRIMAL_FSB_DEFAULT = 1.0e-6;
#elif FP_PREC == 0
static constexpr fp_t NNLS_RESID_TOL_DEFAULT = 1.0e-5;
static constexpr fp_t PRIMAL_FSB_DEFAULT = 1.0e-5;
#endif
using matrix_t = std::vector<std::vector<fp_t>>;
using unsg_t = unsigned int;
struct Configuration {
    Configuration():
        linSolverType(LIN_SOLVER_TYPE_DEFAULT), logLevel(LOG_LEVEL_DEFAULT),
        largeBoundsPenalty(USE_LARGE_BOUNDS_PENALTY_DEFAULT),
        posDefCorrection(USE_POS_DEF_CORRECTION_DEFAULT), gammaUpdate(USE_GAMMA_UPDATE_DEFAULT),
#if FP_PREC == 1
        nnlsResidNormFsb(NNLS_RESID_TOL_DEFAULT), origPrimalFsb(PRIMAL_FSB_DEFAULT)
#elif FP_PREC == 0
        nnlsResidNormFsb(NNLS_RESID_TOL_DEFAULT), origPrimalFsb(PRIMAL_FSB_DEFAULT)
#endif
    {}
    unsigned char linSolverType;
    unsigned char logLevel;
    bool largeBoundsPenalty;
    bool posDefCorrection;
    bool gammaUpdate;
    bool rejectSingular = false;
    fp_t nnlsResidNormFsb; // infeasibility criterion
    fp_t origPrimalFsb;    // Ax - b <= origPrimalFsb
};

// 0.5 * x_T * H * x + c_T * x
// Ax <= b; lw <= x <= up
// First nEqConstraints in A are equality constraints
struct Input {
    matrix_t H;
    matrix_t A;
    std::vector<fp_t> b;
    std::vector<fp_t> c;
    std::vector<fp_t> lw;
    std::vector<fp_t> up;
    unsg_t nEqConstraints = 0;
};

struct Output {
    bool isPositiveDefinite;
    unsigned char exitStatus;
    unsg_t nIterations;
    fp_t cost;
    std::vector<fp_t> x;
    std::vector<fp_t> lambdaC;
    std::vector<fp_t> lambdaLw;
    std::vector<fp_t> lambdaUp;
};


} // namespace QP_NNLS
#endif
