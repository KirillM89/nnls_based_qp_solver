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
    std::vector<fp_t> solution;
    std::list<unsg_t> indices;
};

struct LinSolverTime {
    unsigned long long us;
    unsg_t nConstraints;
};

} //

#endif // CORE_TYPES_H
