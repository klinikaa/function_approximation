#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include "approximation.h"

typedef std::pair<std::vector<float>, std::vector<float>> Points;

void labInfo() {
    std::cout << "==============================" << std::endl;
    std::cout << "CHESNOKOV ARKADY" << std::endl;
    std::cout << "P33111 GROUP" << std::endl;
    std::cout << "INTERPOLATION AND APPROXIMATION" << std::endl;
    std::cout << "==============================" << std::endl;
}

Points readFunctionPointsFromFile(std::string &fileName) {
    std::ifstream file(fileName);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open the file!");
    }

    std::string line;
    std::getline(file, line);

    std::istringstream xs_stream(line);
    std::vector<float> xs;

    float x;
    while (xs_stream >> x) {
        xs.push_back(x);
    }

    std::getline(file, line); // Считываем вторую строку с массивом ys
    std::istringstream ysStream(line);
    std::vector<float> ys;

    float y;
    while (ysStream >> y) {
        ys.push_back(y);
    }

    return {xs, ys};
}

int main() {
    std::string fileName = "test.txt";
    Points points = readFunctionPointsFromFile(fileName);

    if (points.first.size() != points.second.size()) {
        throw std::runtime_error("The number of points x and y don't match!");
    }

    std::pair<float, float> ab = approx_lineal(points.first, points.second);
    float a = ab.first;
    float b = ab.second;

    std::cout << "We got a = " + std::to_string(a) + " and b = " + std::to_string(b) << std::endl;

    return 0;
}
