#ifndef NNLS_CORE_H
#define NNLS_CORE_H
#include <memory>
#include <map>
#include "types.h"
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
        std::vector<double> s;      // nc
        std::vector<double> zp;     // nc
        std::vector<double> primal; // nc
        std::vector<double> dual;   // nc
        std::vector<double> lambda; // nc
        std::vector<double> aux;    // nc
        std::vector<double> x;      // nv
        std::vector<double> v;      // nv
        std::vector<double> Chol;   // nv
        std::vector<bool> dCorrected; // nv
        std::set<unsigned int> activeConstraints;    // max nc
        std::set<unsigned int> linEqConstraints;     // max nc - 2 * nv
        std::unordered_set<unsigned int> negativeZp; // max nc
        matrix_t M;  // nc x nv
        std::deque<unsg_t> addHistory;
    };

    enum PrimalRetStatus {
        SINGULARITY = 0x01,
        LINE_SEARCH_FAILED = 0x02
    };

public:
    Core();
    ~Core() = default;
    void Set(const CoreSettings& settings);
    void ResetProblem();
    void SetCallback(std::unique_ptr<Callback> callback);
    bool InitProblem(const DenseQPProblem& problem);
    void Solve();
    const SolverOutput& GetOutput() { return output; }
    InitStageStatus GetInitStatus() { return initStatus; }
private:
    unsg_t nVariables;
    unsg_t nConstraints;
    unsg_t nPVariables;
    unsg_t nPConstraints;
    unsg_t nEqConstraints;
    unsg_t newActiveIndex;
    unsg_t rptInterval;
    std::unordered_set<unsg_t> singularIndices;
    unsg_t dualIteration;
    unsg_t primalIteration;
    DualLoopExitStatus dualExitStatus;
    PrimalLoopExitStatus primalExitStatus;
    double mEps;
    double minEl;
    double gamma;
    double gammaCorrection;
    double styGamma;
    double scaleFactorDB;
    double rsNorm;
    double newActive;
    double dualTolerance;
    double cost;
    CoreSettings settings;
    WorkSpace ws;
    std::unique_ptr<iTimer> timer;
    std::unique_ptr<Callback> uCallback;
    std::unique_ptr<ILinSolver> lSolver;
    SolverOutput output;
    InitStageStatus initStatus;
    std::vector<LinSolverTime> linSolverTimes;
    bool PrepareNNLS(const DenseQPProblem& problem);
    bool OrigInfeasible();
    bool FullActiveSet();
    bool SkipCandidate(unsg_t indx);
    bool MakeLineSearch();
    bool IsCandidateForNewActive(unsg_t index, double toCompare, bool skip = true);
    void TimeInterval(std::string& buf);
    void UnscaleD();
    void ComputeDualVariable();
    void UpdateGammaOnPrimalIteration();
    void UpdateGammaOnDualIteration();
    void AddToActiveSet(unsg_t indx);
    void RmvFromActiveSet(unsg_t indx);
    void ComputeOrigSolution();
    void FillOutput();
    void SetInitData(const DenseQPProblem &problem);
    void SetLinearSolver();
    void SetIterationData();
    void SetFinalData();
    void AllocateWs();
    void UpdateScaleFactor(std::size_t iConstraint);
    void AddExtraComponent(unsg_t indx, double bound, int& status);
    void ComputeLS4Constraints(const matrix_t& ld, const DenseQPProblem& problem);
    void ComputeLS4Bounds(const DenseQPProblem& problem);
    void Scale();
    void Init(const DenseQPProblem& problem);
    matrix_t ComputeLDLT(const matrix_t& H);
    unsg_t SelectNewActiveComponent();
    unsg_t SolvePrimal();
    int UpdatePrimal();
    //auxilary methods
    void CheckFactorization(const matrix_t& ld, const matrix_t& H, const std::vector<int>& pmt);
};
}
#endif // CORE_H
