#ifndef QP_TYPES_H
#define QP_TYPES_H
#include <string>
#include <vector>
#include <limits>
namespace QP_SOLVERS {
using matrix_t = std::vector<std::vector<double>>;
enum class QpSolverType {
    NNLS = 0,
    QLD,
    DAQP,
    UNKNOWN,
};
enum class SolverRetStatus {
    SUCCESS = 0,
    PREPROCESSING,
    INFEASIBILITY,
    CYCLING,
    ITERATIONS,
    NUM_PROBLEMS,
    NONCONVEX,
    OTHER,
};

struct SolverOutput {
    std::size_t nIterations;
    std::string msg;
    SolverRetStatus status;
    double cost;
    unsigned long long time;
    std::vector<double> primal;
    std::vector<double> dualA;
    std::vector<double> dualL;
    std::vector<double> dualU;
};

struct QpMetrics {
   QpMetrics() {
        Reset();
    }
    double dualityGap;    //duality gap
    double maxPrInfsbB;   //maximum primal infsb for bounds
    double maxPrInfsbIEq; //maximum primal infsb for ineq constraints
    double maxPrInfsbEq;  //maximum primal infsb for eq constraints
    double maxNegDl;      //maximun negative dual variable for ineq. constraints
    double maxDlInfsb;    //dual infeasibility
    double primalCost;    //cost of primal problem
    double violatedIEq;   //ineq. constraint with max primal violation
    double violatedB;     //bound with max primal violation
    double violatedEq;    //eq. constraint with max primal violation
    unsigned int mDef;    //0 - positive def, 1 - pos. semidef, 2 - indef.
    unsigned int nConstraints;     //number of ineq. constraints
    unsigned int nEqConstraints;   //number of eq. constraints
    unsigned int nVariables;       //number of variables
    unsigned int nPrInfsbB;        //number of bounds with primal violation
    unsigned int nPrInfsbIEq;      // -//- ineq. constraints -//-
    unsigned int nPrInfsbEq;       // -//- eq. constraints -//-
    unsigned int nNegDl;           //number of negative dual for bounds and ineq. constraints
    unsigned int nDlInfsb;         //number of nonzero components of KKT for dual problem
    unsigned int nIterations;      //number of iterations
    unsigned long long time;       //computation time in milliseconds
    std::string errMsg;
    void Reset() {
        dualityGap = std::numeric_limits<double>::max();
        maxPrInfsbIEq = std::numeric_limits<double>::max();
        maxPrInfsbB = std::numeric_limits<double>::max();
        maxNegDl = std::numeric_limits<double>::max();
        maxDlInfsb = std::numeric_limits<double>::max();
        primalCost =  std::numeric_limits<double>::max();
        violatedIEq = 0.0;
        violatedB = 0.0;
        mDef = 0;
        nConstraints = 0;
        nEqConstraints = 0;
        nVariables = 0;
        nPrInfsbIEq = std::numeric_limits<unsigned int>::max();
        nPrInfsbB = std::numeric_limits<unsigned int>::max();
        nNegDl = std::numeric_limits<unsigned int>::max();
        nDlInfsb = std::numeric_limits<unsigned int>::max();
        nIterations = 0;
        time = 0;
        errMsg.clear();
    }
};

struct DenseQpProblem {
    //min 0.5 * x_T * H * x + c_T * x
    //if twoSided == true
    // s.t. bL <= A * x <= bU
    //      lw <= x <= Up
    //      size of b == 2 * size of A
    //      b[2 * i] = bU
    //      b[2 * i + 1] = bL
    //      for eq Constrainst bL = bU
    //if twoSided == false
    // s.t. A * x <= b
    //      lw <= x <= up
    //      size of b == size of A
    // size of lw, up == size of H
    bool twoSided = false;
    matrix_t H;
    matrix_t A;
    std::vector<double> c;
    std::vector<double> b;
    std::vector<double> lw;
    std::vector<double> up;
    std::size_t nEqConstraints;
};

struct DenseProblem {
    matrix_t H;
    matrix_t A;
    std::vector<double> c;
    std::vector<double> b;
    std::vector<double> eqA;
    std::vector<double> lwA;
    std::vector<double> upA;
    std::vector<double> lw;
    std::vector<double> up;
};
}
#endif // QP_TYPES_H
