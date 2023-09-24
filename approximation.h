#ifndef FUNCTION_APPROXIMATION_APPROXIMATION_H
#define FUNCTION_APPROXIMATION_APPROXIMATION_H

#include <vector>
#include <numeric>
#include <valarray>
#include <stdexcept>
<<<<<<< HEAD
#include <iostream>

#include "/home/cleanyco/Downloads/eigen-3.4.0/Eigen/Core"
#include "/home/cleanyco/Downloads/eigen-3.4.0/Eigen/Dense"

std::pair<float, float> approx_lineal(std::vector<float> xs, std::vector<float> ys);

std::vector<float> quadratic_approximation(std::vector<float> xs, std::vector<float> ys);

std::vector<float> cube_approximation(std::vector<float> xs, std::vector<float> ys);

std::pair<float, float> approx_exponential(std::vector<float> &xs, std::vector<float> &ys);

std::pair<float, float> approx_power(std::vector<float> &xs, std::vector<float> &ys);

std::pair<float, float> approx_log(std::vector<float> &xs, std::vector<float> &ys);

=======

std::pair<float, float> approx_lineal(std::vector<float> xs, std::vector<float> ys);

float deviation_lineal(float a, float b, std::vector<float> &xs, std::vector<float> &ys);


float standard_deviation(float S, size_t n);

float correlation_coefficient(std::vector<float> xs, std::vector<float> ys);


std::pair<float, float> approx_exponential(std::vector<float> &xs, std::vector<float> &ys);

float deviation_exponential(float a, float b, std::vector<float> &xs, std::vector<float> &ys);


std::pair<float, float> approx_power(std::vector<float> &xs, std::vector<float> &ys);

float deviation_power(float a, float b, std::vector<float> &xs, std::vector<float> &ys);


std::pair<float, float> approx_log(std::vector<float> &xs, std::vector<float> &ys);

float deviation_log(float a, float b, std::vector<float> &xs, std::vector<float> &ys);

>>>>>>> 14a64c6ca021cf574de1c7c89d39f8eca1b3870f
#endif //FUNCTION_APPROXIMATION_APPROXIMATION_H
