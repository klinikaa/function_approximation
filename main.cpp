#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>

#include "approximation.h"
#include "process.h"
#include "deviation.h"
#include "util.h"
#include "graph.h"

void labInfo() {
    std::cout << "==============================" << std::endl;
    std::cout << "CHESNOKOV ARKADY" << std::endl;
    std::cout << "P33111 GROUP" << std::endl;
    std::cout << "INTERPOLATION AND APPROXIMATION" << std::endl;
    std::cout << "==============================" << std::endl;
}

std::pair<std::vector<float>, std::vector<float>> readFunctionPointsFromFile(std::string &fileName) {
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

    std::getline(file, line);
    std::istringstream ysStream(line);
    std::vector<float> ys;

    float y;
    while (ysStream >> y) {
        ys.push_back(y);
    }

    return {xs, ys};
}

int main() {
    labInfo();

    std::string fileName = "test.txt";
    std::pair<std::vector<float>, std::vector<float>> points = readFunctionPointsFromFile(fileName);

    if (points.first.size() != points.second.size()) {
        throw std::runtime_error("The number of points x and y don't match!");
    }

    std::vector<float> xs = points.first;
    std::vector<float> ys = points.second;

    bool isNegativeX = hasNegativeNumber(xs);
    bool isNegativeY = hasNegativeNumber(ys);

    std::vector<float> result;
    if (isNegativeX && isNegativeY) {
        result.push_back(process_lineal(xs, ys));
        result.push_back(process_quadratic(xs, ys));
        result.push_back(process_qube(xs, ys));

        float minValue = * std::min_element(result.begin(), result.end());

        std::cout << "Best approx: " << minValue << std::endl;

        if (minValue == result[0]) {
            std::cout << "The best approximation is lineal" << std::endl;
        } else if (minValue == result[1]) {
            std::cout << "The best approximation is quadratic" << std::endl;
        } else if (minValue == result[2]) {
            std::cout << "The best approximation is qube" << std::endl;
        }

        plotIfXAndYNeg(xs, ys);
    } else if (isNegativeX) {
        result.push_back(process_lineal(xs, ys));
        result.push_back(process_quadratic(xs, ys));
        result.push_back(process_qube(xs, ys));
        result.push_back(process_exp(xs, ys));

        float minValue = * std::min_element(result.begin(), result.end());

        std::cout << "Best approx: " << minValue << std::endl;

        if (minValue == result[0]) {
            std::cout << "The best approximation is lineal" << std::endl;
        } else if (minValue == result[1]) {
            std::cout << "The best approximation is quadratic" << std::endl;
        } else if (minValue == result[2]) {
            std::cout << "The best approximation is qube" << std::endl;
        } else if (minValue == result[3]) {
            std::cout << "The best approximation is exp" << std::endl;
        }

        plotIfXNeg(xs, ys);
    } else if (isNegativeY) {
        result.push_back(process_lineal(xs, ys));
        result.push_back(process_quadratic(xs, ys));
        result.push_back(process_qube(xs, ys));
        result.push_back(process_log(xs, ys));

        float minValue = *std::min_element(result.begin(), result.end());
        std::cout << "Best approx: " << minValue << std::endl;

        if (minValue == result[0]) {
            std::cout << "The best approximation is lineal" << std::endl;
        } else if (minValue == result[1]) {
            std::cout << "The best approximation is quadratic" << std::endl;
        } else if (minValue == result[2]) {
            std::cout << "The best approximation is qube" << std::endl;
        } else if (minValue == result[3]) {
            std::cout << "The best approximation is log" << std::endl;
        }

        plotIfYNeg(xs, ys);
    } else {
        result.push_back(process_lineal(xs, ys));
        result.push_back(process_quadratic(xs, ys));
        result.push_back(process_qube(xs, ys));
        result.push_back(process_power(xs, ys));
        result.push_back(process_exp(xs, ys));
        result.push_back(process_log(xs, ys));

        float minValue = *std::min_element(result.begin(), result.end());

        std::cout << "Best approx: " << minValue << std::endl;

        if (minValue == result[0]) {
            std::cout << "The best approximation is lineal" << std::endl;
        } else if (minValue == result[1]) {
            std::cout << "The best approximation is quadratic" << std::endl;
        } else if (minValue == result[2]) {
            std::cout << "The best approximation is qube" << std::endl;
        } else if (minValue == result[3]) {
            std::cout << "The best approximation is power" << std::endl;
        } else if (minValue == result[4]) {
            std::cout << "The best approximation is exp" << std::endl;
        } else if (minValue == result[5]) {
            std::cout << "The best approximation is log" << std::endl;
        }

        plotAllGraphs(xs, ys);
    }

    return 0;
}