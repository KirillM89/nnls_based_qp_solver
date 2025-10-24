
#ifndef NNLS_QP_SOLVER_UTILS_H
#define NNLS_QP_SOLVER_UTILS_H
#include <cassert>
#include <unordered_set>
#include <unordered_map>
#include <set>
#include <limits>
#include "types.h"
namespace QP_NNLS {
    double g_GetMachineEps();

    void Mult(const matrix_t& M1, const matrix_t& M2, matrix_t& mult); // M1 * M2

    void MultTransp(const matrix_t& M, const std::vector<double>& v, std::vector<double>& MTv); // MTv = M_T * v

    void MultTransp(const matrix_t& M, const std::vector<double>& v, const std::set<unsg_t>& activesetIndices, std::vector<double>& MTv); // M_T * v on active set

    void M1M2T(const matrix_t& M1, const matrix_t& M2, matrix_t& MMT); // MMT = M1 * M2_T

    void Mult(const matrix_t& M, const std::vector<double>& v, std::vector<double>& Mv); // Mv = M * v

    void VSum(const std::vector<double>& v1, const std::vector<double>& v2, std::vector<double>& sum); // sum = v1 + v2

    void VAdd(std::vector<double>& v1, const std::vector<double>& v2); // v1 += v2

    double DotProduct(const std::vector<double>& v1, const std::vector<double>& v2); // <v1,v2>

    double DotProduct(const std::vector<double>& v1, const std::vector<double>& v2,  const std::set<unsg_t>& activesetIndices); // <v1,v2> for active set indices

    void GetIdentityMatrix(int size, matrix_t& M); // set M as identity matrix 

    void InvertHermit(const matrix_t& Chol, matrix_t& Inv); // invert hemitian matrix M using it's Choletsky decomposition M = L * L_T

    void InvertCholetsky(const matrix_t& Chol, matrix_t& Inv); // invert Choletsky decompostion Chol and save results to Inv

    void InvertL(matrix_t& L); // invert matrix L of LDL_T decompostion and save results to right upper part of L

    unsigned char InPlaceLdlt(matrix_t& M,  std::vector<int>& P, size_t& iZero, size_t& iNeg, bool pvt, bool posDefCor);

    matrix_t& operator-(matrix_t& M); // M -> -M in place

    std::vector<double> operator-(const std::vector<double>& v); // return -v

    static inline bool isSame(double cand, double val, double tol = 1.0e-16) {
	assert(tol > 0.0);
	const double diff = cand - val;
	return ((diff >= -tol) && (diff <= tol));
    }

    class LDLT
    {
    public:
        LDLT() = delete;
        LDLT(const matrix_t& M, const std::vector<double>& S);
        virtual ~LDLT() = default;
        unsigned int Compute(const std::set<unsigned int>& activeColumns);
        void Add(std::size_t row, bool delMode = false);
        void AddPvt(std::size_t row, bool delMode = false);
        void Delete(std::size_t row);
        const std::list<unsigned int>& GetRows() { return rows;}
        const matrix_t& GetL() { return L;}
        const std::vector<double>& GetD() { return D;}
        const std::multiset<double, std::greater<double>>& GetPivots() { return pivots;}
        const unsigned int GetNDzero() { return ndzero;}
    private:
        const double dTol = g_GetMachineEps();
        double minD = std::numeric_limits<double>::max();
        const double inf = 1.0e30;
        double d;
        const std::size_t maxSize;
        const std::size_t nX;
        unsigned int ndzero = 0;
        std::size_t curIndex;
        std::size_t actSize;
        const matrix_t& M;
        const std::vector<double>& S;
        matrix_t L;
        std::vector<double> D;
        std::vector<double> norms2;
        std::list<unsigned int> rows;
        std::multiset<double, std::greater<double>> pivots;
        matrix_t cache;
    };

    class MmtLinSolver {
    public:
        MmtLinSolver() = delete;
        MmtLinSolver(const matrix_t& M, const std::vector<double>& S);
        virtual ~MmtLinSolver() = default;
        unsigned int Solve(const std::set<unsigned int>& active, double gamma);
        const std::vector<double>& GetSolution() { return backward;}
    protected:
        const double zeroTol = 1.0e-16;
        unsigned int nDZero;
        const std::size_t maxSize;
        std::size_t curSize;
        double gamma;
        LDLT ldlt;
        const std::vector<double>& S;
        std::vector<double> forward;
        std::vector<double> backward;
        void Forward(const std::set<unsigned int>& active);
        void Backward();
    };

}
#endif
