#include "utils.h"
#include <cmath>
#include <iostream>

namespace QP_NNLS {
    fp_t g_GetMachineEps() {
        fp_t eps = 1.0;
        while (1.0 + eps > 1.0) {
            eps /= 2.0;
        }
        return 4.0 * eps;
    }
	matrix_t& operator-(matrix_t& M) {
		const std::size_t n = M.size();
		if (n == 0) {
			return M;
		}
		const std::size_t m = M.front().size();
		for (std::size_t i = 0; i < n; ++i) {
			for (std::size_t j = 0; j < m; ++j) {
				M[i][j] = -M[i][j];
			}
		}
		return M;
	}
    std::vector<fp_t> operator-(const std::vector<fp_t>& v) {
        const std::size_t n = v.size();
        std::vector<fp_t> vinv(n);
        for (std::size_t i = 0; i < n; ++i) {
            vinv[i] = - v[i];
        }
        return vinv;
    }

	void Mult(const matrix_t& M1, const matrix_t& M2, matrix_t& mult) { //M1*M2
		const std::size_t n1 = M1.size();
		const std::size_t m1 = M1.front().size();
		const std::size_t m2 = M2.front().size();
		for (std::size_t k = 0; k < n1; ++k) {
			for (std::size_t i = 0; i < m2; ++i) {
				mult[k][i] = 0.0;
				for (std::size_t j = 0; j < m1; ++j) {
					mult[k][i] += M1[k][j] * M2[j][i];
				}
			}
		}
	} 
    void MultTransp(const matrix_t& M, const std::vector<fp_t>& v, std::vector<fp_t>& res) { //MT*v
		const std::size_t nrows = M.size();
		const std::size_t ncols = M.front().size();
		for (int i = 0; i < ncols; ++i) {
			res[i] = 0.0;
			for (int j = 0; j < nrows; ++j) {
				res[i] += M[j][i] * v[j];
			}
		}
    }

    void MultTransp(const matrix_t& M, const std::vector<fp_t>& v, const std::set<unsg_t>& activesetIndices, std::vector<fp_t>& res) {
        //M_T * v on active set
        const std::size_t nrows = M.size();
        if (nrows == 0) {
            res.clear();
            return;
        }
        const std::size_t ncols = M.front().size();
        std::fill(res.begin(), res.end(), 0.0);
        if (!activesetIndices.empty()) {
            for (std::size_t i = 0; i < ncols; ++i) {
                res[i] = 0.0;
                for (auto iAct: activesetIndices) {
                    res[i] += M[iAct][i] * v[iAct];
                }
            }
        }
    }

	void M1M2T(const matrix_t& M1, const matrix_t& M2, matrix_t& MMT) {
		//MMT = M1 * M2T
		//m cols M1 = ncols M2
		const int nrM1 = M1.size();
		if (nrM1 == 0) {
			MMT.clear();
			return;
		}
		const int ncM1 = M1.front().size();
		if (ncM1 == 0) {
			MMT.clear();
			return;
		}
		const int nrM2 = M2.size();
		for (int i = 0; i < nrM1; ++i) {
			for (int j = 0; j < nrM2; ++j) {
				MMT[i][j] = 0.0;
				for (int k = 0; k < ncM1; ++k) {
					MMT[i][j] += M1[i][k] * M2[j][k];
				}
			}
		}
	}

    void InvertHermit(const matrix_t& Chol, matrix_t& Inv) {
        //Inv must be allocated in advance
        //M = Chol_T * Chol, Chol - low triangular matrix
        const std::size_t n = Chol.size();
        for (int r = 0; r < n; ++r) {
            const fp_t diagInv = 1.0 / Chol[r][r];
            for (int c = 0; c <= r; ++c) {
                Inv[r][c] = (c == r) ? diagInv : 0.0;
                for (int i = 0; i < r; ++i) {
                    Inv[r][c] -= Chol[r][i] * Inv[i][c];
                }
                Inv[r][c] *= diagInv;
                Inv[c][r] = Inv[r][c];
            }
        }
    }
    void InvertCholetsky(const matrix_t& Chol, matrix_t& Inv) {
        //Inv must be allocated with zeros in advance
        const std::size_t n = Chol.size();
        for (int r = 0; r < n; ++r) {
            Inv[r][r]  = 1.0 / Chol[r][r];
            for (int c = 0; c < r; ++c) {
                Inv[r][c] = 0.0;
                for (int i = c; i < r; ++i) {
                    Inv[r][c] -= Chol[r][i] * Inv[i][c];
                }
                Inv[r][c] *= Inv[r][r];
            }
        }
    }
    void InvertL(matrix_t& L) {
        // Having matrix LD: Low triangular part is matrix L, diagonal - matrix D
        // Compute L_-1 and save to upper triangular part
        const std::size_t n = L.size();
        for (std::size_t r = 0; r < n; ++r) {
            L[r][r] = 1.0;
            for (std::size_t c = 0; c < r; ++c) {
                L[c][r] = 0.0;
                for (int i = c; i < r; ++i) {
                    L[c][r] -= L[r][i] * L[c][i];
                }
                L[r][c] = 0.0;
            }
        }
    }

	void GetIdentityMatrix(int size, matrix_t& M) {
        M.resize(size, std::vector<fp_t>(size, 0.0));
		for (int i = 0; i < size; ++i) {
			M[i][i] = 1.0;
		}
	}

    unsigned char InPlaceLdlt(matrix_t& M, std::vector<int>& P, size_t& iZero, size_t& iNeg, bool pvt, bool posDefCor) {
        const std::size_t n = M.size();
        unsigned char rStatus = 0u;
        const fp_t mEps = g_GetMachineEps();
        const fp_t dCorrection = 2.0 * mEps;
        const fp_t machZero = 0.1 * mEps;
        for (std::size_t c = 0; c < n; ++c) { // by columns
            P[c] = -1.0;
            if (pvt) {
                std::size_t iMax = c;
                fp_t dmax = M[c][c];
                for (std::size_t i = c; i < n; ++i) {
                    if (M[i][i] > dmax) {
                        iMax = i;
                        dmax = M[i][i];
                    }
                }
                if (iMax != c) {
                    std::swap(M[c], M[iMax]); // swap rows
                    for (std::size_t i = 0; i < n; ++i) {
                        std::swap(M[i][c], M[i][iMax]); // swap cols
                    }
                    P[c] = iMax;
                }
            }
            fp_t d = 0.0;
            for (std::size_t r = 0; r < c; ++r) { // by rows
                d -= (M[c][r] * M[c][r] * M[r][r]); // d_rr
            }
            M[c][c] += d; //save d[c]
            if (isSame(M[c][c], 0.0, machZero)) {
                if (posDefCor) {
                    M[c][c] = dCorrection;
                }
                else {
                    rStatus |= 1u;
                }
            } else if (M[c][c] < 0.0) {
                rStatus |= (1u << 1);
            }
            // go by rows in low triangular part and rewrite M[r][c] with L[r][c], r < c
            // M[r][r] = 1.0
            // L_ij = (M_ij - Sum_k=1:j L_ik * L_jk * D_k) / D_j
            for (std::size_t r = c + 1; r < n; ++r) {
                fp_t sum = 0.0;
                for (std::size_t k = 0; k < c; ++k) {
                    sum -= (M[r][k] * M[c][k] * M[k][k]);
                }
                if (!isSame(M[c][c], 0.0, machZero)) {
                    M[r][c] = (M[r][c] + sum) / M[c][c];
                } else {
                    M[r][c] = 0.0; // minor [0:c]x[0:c] is zero, M is indefinite or semidefinite matrix
                }
            }
        }
        return rStatus;
    }

    void Mult(const matrix_t& M, const std::vector<fp_t>& v, std::vector<fp_t>& res) {
        //M_T * v
		const int n = M.size();
		const int m = M.front().size();
		for (int i = 0; i < n; ++i) {
			res[i] = 0.0;
			for (int j = 0; j < m; ++j) {
				res[i] += M[i][j] * v[j];
			}
		}
	}
    void VSum(const std::vector<fp_t>& v1, const std::vector<fp_t>& v2, std::vector<fp_t>& sum) { //v1+v2
		const int n = static_cast<int>(v1.size());
		for (int i = 0; i < n; ++i) {
			sum[i] = v1[i] + v2[i];
		}
	}
    void VAdd(std::vector<fp_t>& v1, const std::vector<fp_t>& v2) { // v1+=v2
		return;
	}
    fp_t DotProduct(const std::vector<fp_t>& v1, const std::vector<fp_t>& v2) {
        fp_t res = 0.0;
		const int sz = v1.size();
		for (int i = 0; i < sz; ++i) {
			res += v1[i] * v2[i];
		}
		return res;
	}

    fp_t DotProduct(const std::vector<fp_t>& v1, const std::vector<fp_t>& v2, const std::set<unsg_t>& activeSetIndices) {
        fp_t res = 0.0;
        for (auto iAct: activeSetIndices) {
            res += v1[iAct] * v2[iAct];
        }
        return res;
    }

    LDLT::LDLT(const matrix_t& M, const std::vector<fp_t>& S):
        d(0.0), maxSize(M.size()), nX(M.front().size()), curIndex(0), actSize(0),
        M(M), S(S),
        L(matrix_t(maxSize, std::vector<fp_t>(maxSize))),
        D(std::vector<fp_t>(maxSize)),
        norms2(std::vector<fp_t>(maxSize)),
        cache(matrix_t(maxSize, std::vector<fp_t>(maxSize, inf)))
    {
        for (std::size_t r = 0; r < maxSize; ++r) {
            fp_t sum = 0.0;
            for (std::size_t c = 0; c < nX ; ++c) {
                sum += M[r][c] * M[r][c];
            }
            norms2[r] = sum + S[r] * S[r];
        }
    }

    unsigned int LDLT::Compute(const std::set<unsigned int>& active) {
        //clear cache
        rows.clear();
        actSize = 0;
        ndzero = 0;
        for (std::size_t i = 0; i < maxSize; ++i) {
            for (std::size_t j = 0; j < maxSize; ++j) {
                cache[i][j] = inf;
            }
        }
        for (auto iAct : active) {
            Add(iAct);
        }
        return ndzero;
    }
    void LDLT::AddPvt(std::size_t rowNumber, bool delMode) {
        std::size_t pos = std::distance(pivots.begin(), pivots.insert(norms2[rowNumber]));
        auto itr = rows.begin();
        std::advance(itr, pos);
        std::vector<unsigned int> r(rows.begin(), rows.end());
        rows = std::list<unsigned int>(rows.begin(), itr);
        actSize = pos;
        Add(rowNumber, delMode);
        for (auto it = r.begin() + pos; it != r.end(); ++it) {
            Add(*it, delMode);
        }
    }
    void LDLT::Add(std::size_t rowNumber, bool delMode) {
        // add row with index rowNumber to last position and recompute L and D
        d = norms2[rowNumber];
        // check if row is already added. Note: in deleteMode it's OK
        bool exists = false;
        for (std::list<unsigned int>::iterator it = rows.begin(); it != rows.end(); ++it){
            if (*it == rowNumber) {
                exists = true;
                if (!delMode) {
                    std::cout << "LDLT Add() info: row " << rowNumber << " is already added" << std::endl;
                }
                break;
            }
        }
        if (!exists) {
            rows.push_back(rowNumber);
        }
        std::size_t i = 0;
        for (auto r : rows) {
            if (i >= actSize){
                break;
            }
            if (std::fabs(D[i]) < dTol) {
                L[actSize][i] = 0.0;
            } else {
                fp_t dot = 0.0;
                const fp_t cv = cache[rowNumber][r];
                if (isSame(cv, inf)) {
                    dot = S[r] * S[rowNumber];
                    for (std::size_t j = 0; j < nX; ++j) {
                        dot += M[r][j] * M[rowNumber][j];
                    }
                    cache[rowNumber][r] = dot;
                } else {
                    dot = cv;
                }
                for (std::size_t j = 0; j < i; ++j) {
                    dot -= L[i][j] * D[j] * L[actSize][j];
                }
                L[actSize][i] = dot / D[i];
                d -= L[actSize][i] * dot;
            }
            ++i;
        }
        // fix zero and negative D values
        if (d <= dTol ) {
            std::cout << "LDL warning: " << "d = " << d << " <= eps" << std::endl;
            d = dTol;
            ++ndzero;
        }
        D[actSize] = d;
        L[actSize][actSize] = 1.0;
        ++actSize;
    }

    void LDLT::Delete(std::size_t rowNumber) {
        bool deleted = false;
        actSize = 0;
        for (std::list<unsigned int>::iterator it = rows.begin(); it != rows.end(); ++it){
            if (!deleted && (*it == rowNumber)) {
                ndzero = 0;
                if (!pivots.empty()) {
                    auto itp = pivots.begin();
                    std::advance(itp, std::distance(rows.begin(), it));
                    pivots.erase(itp);
                }
                if((it = rows.erase(it)) == rows.end()) {
                    return;
                }
                std::cout << "LDLT: deleted " << rowNumber << std::endl;
                deleted = true;
            }           
            if (deleted) {           
                Add(*it, true);
            } else {
                ++actSize;
            }
        }
    }

    MmtLinSolver::MmtLinSolver(const matrix_t& M, const std::vector<fp_t>& S):
        nDZero(0), maxSize(S.size()), curSize(0), gamma(1.0), ldlt(M,S), S(S),
        forward(std::vector<fp_t>(maxSize)),
        backward(std::vector<fp_t>(maxSize))
    {}
    unsigned int MmtLinSolver::Solve(const std::set<unsigned int>& active, fp_t gamma) {
        nDZero = ldlt.Compute(active);
        curSize = active.size();
        this->gamma = gamma;
        Forward(active);
        Backward();
        return nDZero;
    }
    void MmtLinSolver::Forward(const std::set<unsigned int>& active) {
        const matrix_t& L = ldlt.GetL();
        std::size_t i = 0;
        for (auto iAct : active) {
            fp_t sum = 0.0;
            for (std::size_t j = 0; j < i; ++j) {
                sum += L[i][j] * forward[j];
            }
            forward[i] = gamma * S[iAct] - sum;
            ++i;
        }
    }
    void MmtLinSolver::Backward() {
        const matrix_t& L = ldlt.GetL();
        const std::vector<fp_t>& D = ldlt.GetD();
        for (int i = curSize - 1; i >= 0; --i) {
            fp_t sum = 0.0;
            for (int j = i + 1; j < curSize; ++j) {
                sum += L[j][i] * D[i] * backward[j];
            }
            backward[i] = (std::fabs(D[i]) < zeroTol) ? 0.0 : (forward[i] - sum) / D[i];
        }
    }

}
