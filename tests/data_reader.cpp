#include "data_reader.h"
namespace INPUT_DATA {

Matrix2D::Matrix2D(std::ifstream& fid):
fid(fid)
{
    try {
        Read();
    } catch(...) {
       fid.close();
       mat.clear();
    }
}

void Matrix2D::Read()
{
    std::string line;
    std::size_t iLine = 0;
    mat.clear();
    while(std::getline(fid, line)) {
        mat.emplace_back();
        std::size_t pos, prev = 0;
        while(pos = line.find(' ', prev) != std::string::npos) {
            if (pos != prev) {
                mat[iLine].push_back(stod(line.substr(prev, pos - prev)));
            }
            prev = pos + 1;
        }
        // if no ' ' after last value
        if (prev < line.size()) {
            mat[iLine].push_back(stod(line.substr(prev)));
        }
        ++iLine;
    }
}

const std::vector<std::vector<double>>& Matrix2D::Get() const
{
    return mat;
}

Vector::Vector(std::ifstream& fid)
{
    try {
        std::string line;
        if (std::getline(fid, line)) {
            vc.clear();
            std::size_t prev, pos = 0;
            while(pos = line.find(' ', prev) != std::string::npos) {
                if (pos != prev) {
                    vc.push_back(stod(line.substr(prev, pos - prev)));
                }
            }
        }
    } catch(...) {
        fid.close();
        vc.clear();
    }
}
const std::vector<double>& Vector::Get() const
{
    return vc;
}


ReadFromFolder::ReadFromFolder(const std::filesystem::path& folder)
{
    try {
        std::filesystem::path filepath = folder;
        if (!std::filesystem::is_directory(filepath)) {
            throw ;
        } else {
            const std::filesystem::path fh = folder / fH;
            const std::filesystem::path fa = folder / fA;
            const std::filesystem::path fc = folder / fC;
            const std::filesystem::path falw = folder / fAlw;
            const std::filesystem::path faup = folder / fAup;
            const std::filesystem::path faeq = folder / fAeq;
            const std::filesystem::path flw = folder / fLw;
            const std::filesystem::path fup = folder / fUp;


            auto hSize = ReadMatrix(fh, problem.H);
            if (hSize.first != hSize.second || !hSize.first) {
                throw std::runtime_error("H matrix size " + std::to_string(hSize.first) + "/" + std::to_string(hSize.second));
            }


            if (std::filesystem::exists(faup)) {
                Vector v(std::ifstream(faup));
                problem.upA = v.Get();
            }
            if ( std::filesystem::exists(faeq)) {
                Vector v(std::ifstream(faeq));
                problem.eqA = v.Get();
            }
            const std::size_t neq = problem.eqA.size();
            std::size_t nlw = problem.lwA.size();
            std::size_t nup = problem.upA.size();
            if (nlw && nup && nlw != nup) {
                throw;
            }
            if(std::filesystem::exists(fa)) {
                Matrix2D A(std::ifstream(fa));
                problem.A = A.Get();
            }
            if (problem.A.size() != std::max(nlw, nup) + neq) {
                throw;
            }
            if (std::filesystem::exists(flw)) {
                Vector v(std::ifstream(flw));
                problem.lw = v.Get();
            }
            if (std::filesystem::exists(fup)) {
                Vector v(std::ifstream(fup));
                problem.up = v.Get();
            }
            nlw = problem.lw.size();
            nup = problem.up.size();
            if (nlw && nup && nlw != nup) {
                throw;
            }
            std::size_t nB = std::max(nlw, nup);
            if (nB && (problem.H.size() != nB)) {
                throw;
            }
            if (std::filesystem::exists(fc)) {
                Vector v(std::ifstream(fup));
                problem.c = v.Get();
            }
            if (problem.c.size() && problem.c.size() != problem.H.size()) {
                throw;
            }
        }
    } catch(...) {

    }
}
std::pair<std::size_t, std::size_t> ReadFromFolder::ReadMatrix(const std::filesystem::path& file, std::vector<std::vector<double>>& mat)
{
    try {
        std::ifstream fid(file);
        if (!fid) {
            throw std::runtime_error("Failed to open file " + std::string(file));
        }
        Matrix2D mt(fid);
        mat = mt.Get();
        fid.close();
    } catch(...) {
        mat.clear();
    }
    std::size_t nCols = mat.size() ? mat.front().size() : 0;
    return std::make_pair<std::size_t, std::size_t>(mat.size(), std::forward(nCols));
}

std::size_t ReadFromFolder::ReadVector(const std::filesystem::path& file, std::vector<double>& vc)
{
    std::ifstream fid(file);
    if (!fid) {
        throw std::runtime_error("Failed to open file " + std::string(file));
    }
    Vector v(fid);
    vc = v.Get();
    fid.close();
    return vc.size();
}



} // INPUT_DATA
