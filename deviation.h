#ifndef FUNCTION_APPROXIMATION_DEVIATION_H
#define FUNCTION_APPROXIMATION_DEVIATION_H

#include <vector>
#include <numeric>
#include <complex>

float average(std::vector<float> &v);

float correlation_coefficient(std::vector<float> xs, std::vector<float> ys);

float standard_deviation(float S, size_t n);

float deviation_lineal(float a, float b, std::vector<float> &xs, std::vector<float> &ys);

float deviation_exponential(float a, float b, std::vector<float> &xs, std::vector<float> &ys);

float deviation_power(float a, float b, std::vector<float> &xs, std::vector<float> &ys);

float deviation_log(float a, float b, std::vector<float> &xs, std::vector<float> &ys);

float deviation_quadratic(float a_0, float a_1, float a_2, std::vector<float> &xs, std::vector<float> &ys);

float deviation_qube(float a_0, float a_1, float a_2, float a_3, std::vector<float> &xs, std::vector<float> &ys);

#endif //FUNCTION_APPROXIMATION_DEVIATION_H
