#ifndef NNLS_TESTS_UTILS_H
#define NNLS_TESTS_UTILS_H

#include <random>
#include <cmath>
#include <memory>
#include <unordered_set>
#include "types.h"
#include "utils.h"
#if FP_PREC == 1
static constexpr double T_ZERO = 1.0e-15;
static constexpr double FP_TOL = 1.0e-7;
static constexpr double FP_REL_TOL = 1.0e-4;
#elif FP_PREC == 0
static constexpr double T_ZERO = 1.0e-6;
static constexpr double FP_TOL = 1.0e-4;
static constexpr double FP_REL_TOL = 1.0e-3;
#endif
static std::random_device rd;  // Will be used to obtain a seed for the random number engine
static std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
using namespace QP_NNLS;
std::vector<double> GenRandomVector(int dim, double lw, double up);
matrix_t GenRandomMatrix(int nRows, int nCols, double lw, double up);
matrix_t GetRandomPosSemidefMatrix(int dim, double lw, double up);
matrix_t genRandomStrictLowerTriangular(int mSize, int iBg);
matrix_t GetIdentity(std::size_t n);
void InvertByGauss(const matrix_t& M, matrix_t& Minv);
bool isNumber(const double& val);
void TestMatrixMult(const matrix_t& m1, const matrix_t& m2, const matrix_t& baseline);
void TestMatrixMultTranspose(const matrix_t& m, const std::vector<double>& v, const std::vector<double> & baseline);
void TestM1M2T(const matrix_t& m1, const matrix_t& m2, const matrix_t& baseline);
void TestLDLT(LDLT& solver, const matrix_t& M,
              const std::vector<double>&S,
              const std::set<unsigned int>& active);
void TestInPlaceLdlt(const matrix_t& M);
double relativeVal(double a, double b);

fp_t RelTol(fp_t v1, fp_t v2);
class LdltTester {
public:
    LdltTester() = default;
    ~LdltTester() = default;
    void Set(const matrix_t& M, const std::vector<double>& S);
    void Add(unsigned int index);
    void AddPvt(unsigned int index);
    void Delete(unsigned int index);
private:
    matrix_t M;
    std::vector<double> S;
    std::unique_ptr<LDLT> solver;
    std::unordered_set<unsigned int> rows;
    std::size_t nR = 0;
    std::size_t nC = 0;
    void Check();
};

class QPProblem {
public: 
    QPProblem() = delete;
    QPProblem(const Input& in, const Output& out, const Configuration& = Configuration{});
    ~QPProblem() = default;
};

#endif

