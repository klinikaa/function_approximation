#ifndef FUNCTION_APPROXIMATION_GRAPH_H
#define FUNCTION_APPROXIMATION_GRAPH_H

#include <vector>
#include "process.h"
#include <sciplot/sciplot.hpp>

void plotAllGraphs(std::vector<float> &xs, std::vector<float> &ys);

void plotIfYNeg(std::vector<float> &xs, std::vector<float> &ys);

void plotIfXNeg(std::vector<float> &xs, std::vector<float> &ys);

void plotIfXAndYNeg(std::vector<float> &xs, std::vector<float> &ys);

#endif //FUNCTION_APPROXIMATION_GRAPH_H
