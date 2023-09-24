#ifndef FUNCTION_APPROXIMATION_PROCESS_H
#define FUNCTION_APPROXIMATION_PROCESS_H

#include "approximation.h"
#include "deviation.h"
#include "table.h"

typedef std::vector<float> Points;
typedef std::pair<float, float> Coefficients;

float process_lineal(Points &xs, Points &ys);

float process_quadratic(Points &xs, Points &ys);

float process_qube(Points &xs, Points &ys);

float process_power(Points &xs, Points &ys);

float process_exp(Points &xs, Points &ys);

float process_log(Points &xs, Points &ys);

#endif //FUNCTION_APPROXIMATION_PROCESS_H
