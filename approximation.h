#ifndef FUNCTION_APPROXIMATION_APPROXIMATION_H
#define FUNCTION_APPROXIMATION_APPROXIMATION_H

#include <vector>
#include <numeric>
#include <valarray>

std::pair<float, float> approx_lineal(std::vector<float> xs, std::vector<float> ys);

float deviation_lineal(float a, float b, std::vector<float> &xs, std::vector<float> &ys);

float standard_deviation(float S, size_t n);

#endif //FUNCTION_APPROXIMATION_APPROXIMATION_H
