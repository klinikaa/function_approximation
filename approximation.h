#ifndef FUNCTION_APPROXIMATION_APPROXIMATION_H
#define FUNCTION_APPROXIMATION_APPROXIMATION_H

#include <vector>
#include <numeric>
#include <valarray>
#include <stdexcept>
#include <iostream>

#include "/home/cleanyco/Downloads/eigen-3.4.0/Eigen/Core"
#include "/home/cleanyco/Downloads/eigen-3.4.0/Eigen/Dense"

std::pair<float, float> approx_lineal(std::vector<float> xs, std::vector<float> ys);

std::vector<float> quadratic_approximation(std::vector<float> xs, std::vector<float> ys);

std::vector<float> cube_approximation(std::vector<float> xs, std::vector<float> ys);

std::pair<float, float> approx_exponential(std::vector<float> &xs, std::vector<float> &ys);

std::pair<float, float> approx_power(std::vector<float> &xs, std::vector<float> &ys);

std::pair<float, float> approx_log(std::vector<float> &xs, std::vector<float> &ys);

#endif //FUNCTION_APPROXIMATION_APPROXIMATION_H
