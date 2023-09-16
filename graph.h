#ifndef FUNCTION_APPROXIMATION_GRAPH_H
#define FUNCTION_APPROXIMATION_GRAPH_H

#include <vector>
#include "../matplotlib-cpp/matplotlibcpp.h"

void setLinearFunction(std::vector<float> &xs, std::vector<float> &ys);

class Graph {
public:
    Graph();

    void setLinearFunction(std::vector<float> &xs, std::vector<float> &ys);

    void plot();
};

#endif //FUNCTION_APPROXIMATION_GRAPH_H
