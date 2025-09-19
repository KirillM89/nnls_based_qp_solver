#include <iostream>
#include <fstream>
#include <sstream>
#include "tools.h"
std::vector<std::vector<double>> readMatrix(const std::string& fileName) {
    std::ifstream file(fileName);
    std::vector<std::vector<double>> matrix;
    if (file.is_open()) {
        std::string line;
        while (std::getline(file, line)) { // Read the file line by line
            std::vector<double> row;
            std::size_t prev = 0, pos;
            while ((pos = line.find(' ', prev)) != std::string::npos) {
                if (pos != prev) {
                    std::string token = line.substr(prev, pos - prev);
                    double value = std::stod(token);
                    row.push_back(value);
                }
                prev = pos + 1;
            }

            // Process the last token
            if (prev < line.length()) {
                std::string token = line.substr(prev);
                double value = std::stod(token);
                row.push_back(value);
            }
            matrix.push_back(row);
        }
    file.close();
    }
    return matrix;
}


