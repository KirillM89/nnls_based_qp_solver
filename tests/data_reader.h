#ifndef DATA_READER_H
#define DATA_READER_H
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <filesystem>
#include "qp_types.h"
namespace INPUT_DATA {
class Matrix2D {
  public:
    Matrix2D() = delete;
    Matrix2D(std::ifstream& fid);
    ~Matrix2D() = default;
    const std::vector<std::vector<double>>& Get() const;
  private:
    std::ifstream& fid;
    std::vector<std::vector<double>> mat;
    void Read();
};

class Vector {
  public:
    Vector() = delete;
    Vector(std::ifstream& fid);
    ~Vector() = default;
    const std::vector<double>& Get() const;
  private:
    std::vector<double> vc;
};

class ReadFromFolder {
// read qp problem from folder
// folder must include files h.txt, c.txt
// aux. files a.txt, a_lw.txt, a_up.txt, a_eq.txt, lw.txt, up.txt
// if a.txt exist then must exist a_eq || a_lw || a_up
  public:
    ReadFromFolder() = delete;
    ReadFromFolder(const std::filesystem::path& folder);
    ~ReadFromFolder() = default;
  private:
    QP_SOLVERS::DenseProblem problem;
    const std::filesystem::path fH = "h.txt";
    const std::filesystem::path fA = "a.txt";
    const std::filesystem::path fC = "c.txt";
    const std::filesystem::path fAlw = "alw.txt";
    const std::filesystem::path fAup = "aup.txt";
    const std::filesystem::path fAeq = "aeq.txt";
    const std::filesystem::path fLw = "lw.txt";
    const std::filesystem::path fUp = "up.txt";
    std::pair<std::size_t, std::size_t> ReadMatrix(const std::filesystem::path& file, std::vector<std::vector<double>>& mat);
    std::size_t ReadVector(const std::filesystem::path& file, std::vector<double>& vc);
};
}

#endif // DATA_READER_H
