#include "linSolvers.h"
#include "utils.h"
#include <cmath>
namespace QP_NNLS {
CumulativeSolver::CumulativeSolver(const matrix_t& M,
                                   const std::vector<double>& s ):
    nConstraints(M.size()),
    nVariables(0),
    nActive(0),
    gamma(1.0),
    M(M),
    s(s)
{
    if (nConstraints > 0) {
        nVariables = M.front().size();
    }
    activeSet.resize(nConstraints, false);
}
bool CumulativeSolver::Add(unsg_t indx) {
    activeSet[indx] = true;
    ++nActive;
    return true;
}
bool CumulativeSolver::Delete(unsg_t indx) {
    if (activeSet[indx]) {
        activeSet[indx] = false;
        if (nActive > 0) {
            --nActive;
        }
    }
    return true;
}
CumulativeLDLTSolver::CumulativeLDLTSolver(const matrix_t& M,
                                           const std::vector<double>& s,
                                           bool rejectSingular):
    rejectSingular(rejectSingular),
    gamma(1.0), ldlt(M, s), ndzero(0),
    maxSize(s.size()), S(s),
    forward(std::vector<double>(maxSize)),
    backward(std::vector<double>(maxSize))
{ }

bool CumulativeLDLTSolver::Add(unsg_t indx) {
    ldlt.Add(indx);
    //ldlt.AddPvt(indx);
    return true;
}
bool CumulativeLDLTSolver::Delete(unsg_t indx) {
    ldlt.Delete(indx);
    return true;
}
const LinSolverOutput& CumulativeLDLTSolver::Solve() {
    const matrix_t& l = ldlt.GetL();
    const std::vector<double>& d = ldlt.GetD();
    const std::list<unsigned int>& rows = ldlt.GetRows();
    output.nDNegative = ldlt.GetNDzero();
    output.indices = rows;
    if (!rejectSingular || output.nDNegative == 0) {
        const std::size_t nr = rows.size();
        output.nDNegative = 0;
        if (nr == 0) {
            output.solution = std::vector<double>(maxSize, 0.0);
            output.indices.clear();
        } else {
            std::size_t i = 0;
            for (auto iAct : rows) {
                double sum = 0.0;
                for (std::size_t j = 0; j < i; ++j) {
                    sum += l[i][j] * forward[j];
                }
                forward[i] = gamma * S[iAct] - sum;
                ++i;
            }
            for (int i = nr - 1; i >= 0; --i) {
                double sum = 0.0;
                for (int j = i + 1; j < nr; ++j) {
                    sum += l[j][i] * d[i] * backward[j];
                }
                if (std::fabs(d[i]) < zeroTol) {
                    backward[i] = 0.0;
                    ++output.nDNegative;
                } else {
                    backward[i] = (forward[i] - sum) / d[i];
                }
            }
            output.solution = backward;
        }
    }
    return output;
}

} //namespace QP_NNLS
