#include "graph.h"

namespace plt = matplotlibcpp;

Graph::Graph() {
    plt::title("I NEED SOME SLEEP");
    plt::legend();
}

void Graph::setLinearFunction(std::vector<float> &xs, std::vector<float> &ys) {
    plt::named_plot("ax + b", xs, ys);
}

void Graph::plot() {
    plt::show();
}

