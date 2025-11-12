#ifndef NNLS_CORE_H
#define NNLS_CORE_H
#include <memory>
#include <map>
#include "core_types.h"
#include "linSolvers.h"
#include "timers.h"
#include "callback.h"
namespace QP_NNLS {
class Core {
    struct WorkSpace {
        // na - number of constraints
        // nv - number of variables
        // nc = na + 2 * nv
        // total size = nc x nv + 8 * nc + 3 * nv ~ na * nv + 2 * nv^2
        std::vector<fp_t> s;      // nc
        std::vector<fp_t> zp;     // nc
        std::vector<fp_t> primal; // nc
        std::vector<fp_t> dual;   // nc
        std::vector<fp_t> lambda; // nc
        std::vector<fp_t> aux;    // nc
        std::vector<fp_t> x;      // nv
        std::vector<fp_t> v;      // nv
        std::vector<fp_t> Chol;   // nv
        std::vector<bool> dCorrected; // nv
        std::set<unsg_t> activeConstraints;    // max nc
        std::set<unsg_t> linEqConstraints;     // max nc - 2 * nv
        std::unordered_set<unsg_t> negativeZp; // max nc
        matrix_t M;  // nc x nv
        std::deque<unsg_t> addHistory;
    };

public:
    Core();
    ~Core() = default;
    void Set(const Configuration& config);
    void ResetProblem();
    void SetCallback(Callback* callback);
    bool InitProblem(const Input& problem);
    void Solve();
    const Output& GetOutput() { return output; }
    InitStageStatus GetInitStatus() { return initStatus; }
private:
    enum LinSolverStatus {
        SINGULARITY = 0x01,
        LINE_SEARCH_FAILED = 0x02,
    };
    unsg_t nVariables = 0;
    unsg_t nConstraints = 0;
    unsg_t nLwBounds = 0;
    unsg_t nUpBounds = 0;
    unsg_t nPVariables = 0;
    unsg_t nPConstraints = 0;
    unsg_t nEqConstraints = 0;
    unsg_t newActiveIndex = 0;
    unsg_t rptInterval = 0;
    unsg_t nDualIterations = 0;
    unsg_t nPrimalIterations = 0;
    std::unordered_set<unsg_t> singularIndices;
    unsg_t dualIteration = 0;
    unsg_t primalIteration = 0;
    DualLoopExitStatus dualExitStatus;
    PrimalLoopExitStatus primalExitStatus;
    fp_t mEps;
    fp_t minEl;
    fp_t gamma;
    fp_t gammaCorrection;
    fp_t styGamma;
    fp_t scaleFactorDB;
    fp_t rsNorm;
    fp_t newActive;
    fp_t dualTolerance;
    fp_t cost;
    Configuration config;
    WorkSpace ws;
    std::unique_ptr<iTimer> timer;
    Callback* uCallback = nullptr;
    std::unique_ptr<ILinSolver> lSolver;
    Output output;
    InitStageStatus initStatus;
    std::vector<LinSolverTime> linSolverTimes;
    bool PrepareNNLS(const Input& problem);
    bool OrigInfeasible();
    bool FullActiveSet();
    bool SkipCandidate(unsg_t indx);
    bool MakeLineSearch();
    bool IsCandidateForNewActive(unsg_t index, fp_t toCompare, bool skip = true);
    void TimeInterval(std::string& buf);
    void UnscaleD();
    void ComputeDualVariable();
    void UpdateGammaOnPrimalIteration();
    void UpdateGammaOnDualIteration();
    void AddToActiveSet(unsg_t indx);
    void RmvFromActiveSet(unsg_t indx);
    void ComputeOrigSolution();
    void FillOutput();
    void SetInitData(const Input &problem);
    void SetLinearSolver();
    void SetIterationData();
    void SetFinalData();
    void AllocateWs();
    void UpdateScaleFactor(std::size_t iConstraint);
    void AddExtraComponent(unsg_t indx, fp_t bound, int& status);
    void ComputeLS4Constraints(const matrix_t& ld, const Input& problem);
    void ComputeLS4Bounds(const matrix_t& ld, const Input& problem);
    void Scale();
    void Init(const Input& problem);
    matrix_t ComputeLDLT(const matrix_t& H);
    unsg_t SelectNewActiveComponent();
    unsg_t SolvePrimal();
    int UpdatePrimal();
    //auxilary methods
    void CheckFactorization(const matrix_t& ld, const matrix_t& H, const std::vector<int>& pmt);
};
}
#endif // CORE_H
