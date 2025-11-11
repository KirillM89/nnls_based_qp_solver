#include "gtest/gtest.h"
#include "test_utils.h"
#include "decorator.h"

bool isNumber(const double& val) {
	return std::isfinite(val);
}

double relativeVal(double a, double b) {
	const double zeroTol = 1.0e-15;
	if (isSame(a, 0.0, zeroTol) && isSame(b, 0.0, zeroTol)) {
		return zeroTol;
    } else if (isSame(a, 0.0, zeroTol) || isSame(b, 0.0, zeroTol)) {
		return std::fabs(a - b);
    } else {
		double diff = a - b;
		if (isSame(diff, 0.0, 1.0e-8)) {
			return zeroTol;
		}
		return std::fabs(diff / b);
	}
}
void InvertByGauss(const matrix_t& M, matrix_t& Minv) {
    const double pivotZero = 1.0e-7;
    const std::size_t n = M.size(); // M must be square square matrix
    if (n <= 1) {
        Minv = M;
        return;
    }
    // Add the identity matrix at the end of original matrix.
    for (int i = 0; i < n; i++) {
        Minv[i][i] = 1.0;
    }
    matrix_t Mtmp = M;
    // Interchange the rows of matrix,
    // interchanging of rows will start from the last row
    for (std::size_t i = n - 1; i > 0; i--) {
        if (std::fabs(Mtmp[i - 1][0]) < std::fabs(Mtmp[i][0])) {
            Mtmp[i].swap(Mtmp[i - 1]);
            Minv[i].swap(Minv[i - 1]);
        }
    }
    for (std::size_t i = 0; i < n; i++) {
        const double tmp = (std::fabs(Mtmp[i][i]) < pivotZero) ? (1.0 / pivotZero) : (1.0 / Mtmp[i][i]);
        for (std::size_t j = 0; j < n; j++) {
            if (j != i) {
                double temp = Mtmp[j][i] * tmp;
                for (int k = 0; k < n; k++) {
                    Minv[j][k] -= Minv[i][k] * temp;
                    Mtmp[j][k] -= Mtmp[i][k] * temp;
                }
            }
        }
    }
    // Multiply each row by a nonzero integer.
    // Divide row element by the diagonal element
    for (std::size_t i = 0; i < n; i++) {
        const double tmp = (std::fabs(Mtmp[i][i]) < pivotZero) ? (1.0 / pivotZero) : (1.0 / Mtmp[i][i]);
        for (std::size_t j = 0; j < n; j++) {
            Mtmp[i][j] = Mtmp[i][j] * tmp;
            Minv[i][j] = Minv[i][j] * tmp;
        }
    }
}
void TestMatrixMult(const matrix_t& m1, const matrix_t& m2, const matrix_t& baseline) {
    ASSERT_EQ(m1.front().size(), m2.size());
	matrix_t mult(m1.size(), std::vector<double>(m2.front().size(), 1.0)); // non-zero output
	Mult(m1, m2, mult);
	for (int i = 0; i < m1.size(); ++i) {
		for (int j = 0; j < m2.front().size(); ++j) {
			EXPECT_EQ(baseline[i][j], mult[i][j]);
		}
	}
}
void TestMatrixMultTranspose(const matrix_t& m, const std::vector<double>& v, const std::vector<double>& baseline) {
	ASSERT_EQ(m.size(), v.size());
	std::vector<double> res(m.front().size(), 1.0);
    MultTransp(m, v, res);
	for (std::size_t i = 0; i < m.front().size(); ++i) {
		EXPECT_EQ(res[i], baseline[i]);
	}
}

void TestM1M2T(const matrix_t& m1, const matrix_t& m2, const matrix_t& baseline) {
	ASSERT_EQ(m1.front().size(), m2.front().size());
	matrix_t m1m2T(m1.size(), std::vector<double>(m2.size()));
	M1M2T(m1, m2, m1m2T);
	const double tol = 1.0e-10;
	for (int i = 0; i < m1.size(); ++i) {
		for (int j = 0; j < m2.size(); ++j) {
			EXPECT_NEAR(m1m2T[i][j], baseline[i][j], tol);
		}
	}
}

void TestLDLT(LDLT& solver, const matrix_t& M,
              const std::vector<double>&S,
              const std::set<unsigned int>& active) {
    const std::size_t nR = M.size();
    const std::size_t nC = M.front().size();
    const std::size_t nActive = active.size();
    ASSERT_EQ(S.size(), nR);
    ASSERT_LE(nActive, nR);
    solver.Compute(active);
    const matrix_t& l = solver.GetL();
    const std::vector<double>& d = solver.GetD();
    matrix_t mmt(nActive, std::vector<double>(nActive));
    std::size_t i = 0;
    for (unsigned int iA : active) {
        std::size_t j = 0;
        for (unsigned int jA : active) {
            double sum = S[iA] * S[jA];
            for (std::size_t c = 0; c < nC; ++c) {
                sum += M[iA][c] * M[jA][c];
            }
            mmt[i][j++] = sum;
        }
        ++i;
    }
    matrix_t ld(nActive, std::vector<double>(nActive, 0.0));
    for (unsigned int i = 0; i < nActive; ++i) {
        for (unsigned int j = 0; j < nActive; ++j) {
            ld[i][j] = l[i][j] * d[j];
        }
    }
    matrix_t ldl(mmt);
    for (unsigned int i = 0; i < nActive; ++i) {
        for (unsigned int j = 0; j < nActive; ++j) {
            double sum = 0.0;
            for (std::size_t c = 0; c < nActive; ++c) {
                sum += ld[i][c] * l[j][c];
            }
            EXPECT_NEAR(mmt[i][j], sum, 1.0e-7);
        }
    }
}

void TestInPlaceLdlt(const matrix_t& M) {
    const std::size_t n = M.size();
    matrix_t MCP = M;
    std::vector<double> D(n);
    std::vector<int> P(n, -1.0);
    size_t iZero = 2;
    size_t iNeg = 1;
    bool pvt = true;
    bool posDefCor = false;
    InPlaceLdlt(MCP, P, iZero, iNeg, pvt, posDefCor);
    if (!pvt) {
        for (std::size_t i = 0; i < n; ++i) {
            ASSERT_EQ(P[i], -1);
        }
    }
    matrix_t L(MCP);
    std::cout << "zn " << iZero << " " << iNeg << std::endl;
    std::cout << "d: ";
    for (std::size_t i = 0; i < n; ++i) {
        D[i] = L[i][i];
        L[i][i] = 1.0;
        for (std::size_t j = i + 1; j < n; ++j) {
            L[i][j] = 0.0;
        }
    }
    std::cout << std::endl;
    matrix_t LD(n, std::vector<double>(n));
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < n; ++j) {
            LD[i][j] = L[i][j] * D[j];
        }
    }
    matrix_t LDLT(n, std::vector<double>(n));
    M1M2T(LD, L, LDLT);
    //apply permutation
    for (int i = n - 1; i >= 0; --i) {
        if (P[i] != -1) {
            std::swap(LDLT[i], LDLT[P[i]]);
            for (std::size_t j = 0; j < n; ++j) {
                std::swap(LDLT[j][i], LDLT[j][P[i]]);
            }
        }
    }
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < n; ++j) {
            EXPECT_LT(relativeVal(M[i][j], LDLT[i][j]), 1.0e-6) << "res=" << LDLT[i][j] << " baseline="
                    << M[i][j] << " i-j " << i << " " << j;
        }
    }
}

std::vector<double> GenRandomVector(int dim, double lw, double up) {
	std::uniform_real_distribution<> dis(lw, up);
	std::vector<double> v(dim);
	for (int i = 0; i < dim; ++i) {
		v[i] = dis(gen);
	}
	return v;
}
matrix_t GenRandomMatrix(int nRows, int nCols, double lw, double up) {
	std::uniform_real_distribution<> dis(lw, up);
	matrix_t M(nRows, std::vector<double>(nCols));
	for (int i = 0; i < nRows; ++i) {
		for (int j = 0; j < nCols; ++j) {
			M[i][j] = dis(gen);
			if (i == j && std::fabs(M[i][j]) < 1.0e-7) {
				M[i][j] += 1.0e-7;
			}
		}
	}
	return M;
}
matrix_t GetRandomPosSemidefMatrix(int dim, double lw, double up) {
	matrix_t Mrand = GenRandomMatrix(dim, dim, lw, up);
	matrix_t MPS = Mrand;
	M1M2T(Mrand, Mrand, MPS);
	return MPS;
}
matrix_t genRandomStrictLowerTriangular(int mSize, int iBg) {
	const int lw = -100;
	const int up = 100;
	matrix_t Mrand = GenRandomMatrix(mSize, mSize, lw, up);
	int k = 0;
	for (int i = 0; i < mSize; ++i) {
		const int jBg = (i >= iBg) ? ++k : 0;
		for (int j = jBg; j < mSize; ++j) {
			Mrand[i][j] = 0.0;
		}
	}	
	return Mrand;
}
matrix_t GetIdentity(std::size_t n)
{
    matrix_t I(n);
    for (std::size_t i = 0; i < n; ++i) {
        I[i].resize(n, 0);
        I[i][i] = 1.0;
    }
    return I;
}
fp_t RelTol(fp_t v1, fp_t v2)
{
    if (std::fabs(v1) < T_ZERO && std::fabs(v2) < T_ZERO) {
        return 0;
    }
    if (std::fabs(v2) < T_ZERO) {
        std::swap(v2, v1);
    }
    return std::fabs((v1 - v2) / v2);
}
void LdltTester::Set(const matrix_t& M, const std::vector<double>& S) {
    this->M = M;
    this->S = S;
    nR = M.size();
    nC = M.front().size();
    solver = std::make_unique<LDLT>(M, S);
}
void LdltTester::Add(unsigned int index) {
    rows.insert(index);
    ASSERT_LT(index, nR);
    solver->Add(index);
    Check();
}
void LdltTester::AddPvt(unsigned int index) {
    rows.insert(index);
    ASSERT_LT(index, nR);
    solver->AddPvt(index);
    Check();
}
void LdltTester::Delete(unsigned int index) {
    rows.erase(index);
    ASSERT_LT(index, nR);
    solver->Delete(index);
    Check();
}
void LdltTester::Check() {
    const auto& rows = solver->GetRows();
    const std::size_t nr = rows.size();
    const auto pivots = solver->GetPivots();
    if (!pivots.empty()) {
        ASSERT_EQ(pivots.size(), nr);
        double pvt = std::numeric_limits<double>::max();
        for (auto el : pivots){
            ASSERT_LE(el, pvt);
            pvt = el;
        }
    }
    ASSERT_EQ(nr, this->rows.size());
    matrix_t mmt(nr, std::vector<double>(nr));
    std::size_t i = 0;
    for (unsigned int iR : rows) {
        std::size_t j = 0;
        for (unsigned int jR : rows) {
            double sum = S[iR] * S[jR];
            for (std::size_t c = 0; c < nC; ++c) {
                sum += M[iR][c] * M[jR][c];
            }
            mmt[i][j++] = sum;
        }
        ++i;
    }
    const matrix_t& l = solver->GetL();
    const std::vector<double>& d = solver->GetD();
    matrix_t ld(nr, std::vector<double>(nr, 0.0));
    for (unsigned int i = 0; i < nr; ++i) {
        for (unsigned int j = 0; j < nr; ++j) {
            ld[i][j] = l[i][j] * d[j];
        }
    }
    matrix_t ldl(mmt);
    for (unsigned int i = 0; i < nr; ++i) {
        for (unsigned int j = 0; j < nr; ++j) {
            double sum = 0.0;
            for (std::size_t c = 0; c < nr; ++c) {
                sum += ld[i][c] * l[j][c];
            }
            EXPECT_NEAR(mmt[i][j], sum, 1.0e-7);
        }
    }
}
// set problem
// method SetProblem() preprocesses the input problem.
// Returns true if all the preprocessing procedures passed without numerical problems, false otherwise
// If SetProblem() returns "true" the properties of the problem can be seen using method GetInitStatus()
// QP_NNLS::InitStageStatus initStatus = solver.GetInitStatus()
// "SUCCESS" means that problem has positive definite H
// "D_Z" means positive semidefinite
// "D_ZN" or "D_N"  - indefinite or negative definite
// Solver can't solve indefinite or negative definite problems but can solve positive semidefinite problems,
// in this case H correction will take place
QPProblem::QPProblem(const Input& in, const Output& out, const Configuration& config)
{
    using namespace QP_NNLS;
    QPNNLS solver;    // create solver instance
    solver.Init(config);       // apply settings
    bool posDef = true;
    if (!solver.SetProblem(in)) {
        if (solver.GetInitStatus()) {
            posDef = false;
        }
    }
    solver.Solve();
    const QP_NNLS::Output& sOut = solver.GetOutput();
    EXPECT_EQ(sOut.dualExitStatus, out.dualExitStatus);
    EXPECT_EQ(sOut.primalExitStatus, out.primalExitStatus);
    EXPECT_EQ(sOut.isPositiveDefinite, out.isPositiveDefinite);
    EXPECT_EQ(posDef, out.isPositiveDefinite);
    EXPECT_EQ(sOut.nIterations, out.nIterations);
    EXPECT_EQ(sOut.x.size(), out.x.size());
    EXPECT_EQ(sOut.lambdaC.size(), out.lambdaC.size());
    EXPECT_EQ(sOut.lambdaLw.size(), out.lambdaLw.size());
    EXPECT_EQ(sOut.lambdaUp.size(), out.lambdaUp.size());
    EXPECT_NEAR(sOut.maxBViolation, out.maxBViolation, FP_TOL);
    EXPECT_NEAR(sOut.maxCViolation, out.maxCViolation, FP_TOL);
    if (sOut.x.size() == out.x.size()) {
        for (std::size_t i = 0; i < out.x.size(); ++i) {
            EXPECT_LT(RelTol(sOut.x[i], out.x[i]), FP_REL_TOL);
        }
    }
    if (sOut.lambdaC.size() == out.lambdaC.size()) {
        for (std::size_t i = 0; i < out.lambdaC.size(); ++i) {
            EXPECT_LT(RelTol(sOut.lambdaC[i], out.lambdaC[i]), FP_REL_TOL);
        }
    }
    if (sOut.lambdaLw.size() == out.lambdaLw.size()) {
        for (std::size_t i = 0; i < out.lambdaLw.size(); ++i) {
            EXPECT_LT(RelTol(sOut.lambdaLw[i], out.lambdaLw[i]), FP_REL_TOL);
        }
    }
    if (sOut.lambdaUp.size() == out.lambdaUp.size()) {
        for (std::size_t i = 0; i < out.lambdaUp.size(); ++i) {
            EXPECT_LT(RelTol(sOut.lambdaUp[i], out.lambdaUp[i]), FP_REL_TOL);
        }
    }
}











