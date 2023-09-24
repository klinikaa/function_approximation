#include "graph.h"

void plotAllGraphs(std::vector<float> &xs, std::vector<float> &ys) {
    Coefficients cf = approx_lineal(xs, ys);
    float a = cf.first;
    float b = cf.second;

    auto phi_of_x_lin = [a, b](float x) -> float {
        return a * x + b;
    };

    std::vector<float> phi_lin;
    for (float x : xs) {
        phi_lin.push_back(phi_of_x_lin(x));
    }

    std::vector<float> acf = quadratic_approximation(xs, ys);
    float a_0 = acf[0];
    float a_1 = acf[1];
    float a_2 = acf[2];

    auto phi_of_x_quad = [a_0, a_1, a_2](float x) -> float {
        return a_2 * std::pow(x, 2) + a_1 * std::pow(x, 1) + a_0 * std::pow(x, 0);
    };

    std::vector<float> phi_quad;
    for (float x : xs) {
        phi_quad.push_back(phi_of_x_quad(x));
    }

    acf = cube_approximation(xs, ys);
    float a_00 = acf[0];
    float a_10 = acf[1];
    float a_20 = acf[2];
    float a_30 = acf[3];

    auto phi_of_x_cube = [a_00, a_10, a_20, a_30](float x) -> float {
        return a_30 * std::pow(x, 3) + a_20 * std::pow(x, 2) + a_10 * std::pow(x, 1) + a_00 * std::pow(x, 0);
    };

    std::vector<float> phi_cube;
    for (float x : xs) {
        phi_cube.push_back(phi_of_x_cube(x));
    }

    cf = approx_power(xs, ys);
    float ap = cf.first;
    float bp = cf.second;

    auto phi_of_x_power = [ap, bp](float x) -> float {
        return ap * std::pow(x, bp);
    };

    std::vector<float> phi_power;
    for (float x : xs) {
        phi_power.push_back(phi_of_x_power(x));
    }

    cf = approx_exponential(xs, ys);
    float ae = cf.first;
    float be = cf.second;

    auto phi_of_x_exp = [ae, be](float x) -> float {
        return ae * std::exp(be * x);
    };

    std::vector<float> phi_e;
    for (float x : xs) {
        phi_e.push_back(phi_of_x_exp(x));
    }

    cf = approx_log(xs, ys);
    float al = cf.first;
    float bl = cf.second;

    auto phi_of_x_log = [al, bl](float x) -> float {
        return al * std::log(x) + bl;
    };

    std::vector<float> phi_log;
    for (float x : xs) {
        phi_log.push_back(phi_of_x_log(x));
    }

    using namespace sciplot;

    Plot2D plot;

    plot.xlabel("x");
    plot.ylabel("y");
    plot.grid().show();

    plot.legend()
            .atOutsideBottom()
            .displayHorizontal()
            .displayExpandWidthBy(2);

    plot.drawCurve(xs, ys).label("y = sf").lineColor("#5B0888");
    plot.drawCurve(xs, phi_lin).label("y = " + std::to_string(a) + "x + " + std::to_string(b)).lineColor("#7D7C7C");
    plot.drawCurve(xs, phi_log).label("y = " + std::to_string(al) + "ln(x) + " + std::to_string(bl)).lineColor("#79AC78");
    plot.drawCurve(xs, phi_e).label("y = " + std::to_string(ae) + "e^(" + std::to_string(b) + "x)").lineColor("#E55604");
    plot.drawCurve(xs, phi_quad).label(std::to_string(a_0) + std::to_string(a_1) + "x + " + std::to_string(a_2) + "x^2").lineColor("#26577C");
    plot.drawCurve(xs, phi_cube).label(std::to_string(a_00) + std::to_string(a_10) + "x + " + std::to_string(a_20) + "x^2 + "
                                       + std::to_string(a_30) + "x^3").lineColor("#EDB7ED");
    plot.drawCurve(xs, phi_power).label(std::to_string(ap) + "x^(" + std::to_string(bp) + ")").lineColor("#EF9595");

    Figure fig = {{plot}};
    Canvas canvas = {{fig}};

    canvas.size(800, 600);

    canvas.show();
}

void plotIfXAndYNeg(std::vector<float> &xs, std::vector<float> &ys) {
    Coefficients cf = approx_lineal(xs, ys);
    float a = cf.first;
    float b = cf.second;

    auto phi_of_x_lin = [a, b](float x) -> float {
        return a * x + b;
    };

    std::vector<float> phi_lin;
    for (float x : xs) {
        phi_lin.push_back(phi_of_x_lin(x));
    }

    std::vector<float> acf = quadratic_approximation(xs, ys);
    float a_0 = acf[0];
    float a_1 = acf[1];
    float a_2 = acf[2];

    auto phi_of_x_quad = [a_0, a_1, a_2](float x) -> float {
        return a_2 * std::pow(x, 2) + a_1 * std::pow(x, 1) + a_0 * std::pow(x, 0);
    };

    std::vector<float> phi_quad;
    for (float x : xs) {
        phi_quad.push_back(phi_of_x_quad(x));
    }

    acf = cube_approximation(xs, ys);
    float a_00 = acf[0];
    float a_10 = acf[1];
    float a_20 = acf[2];
    float a_30 = acf[3];

    auto phi_of_x_cube = [a_00, a_10, a_20, a_30](float x) -> float {
        return a_30 * std::pow(x, 3) + a_20 * std::pow(x, 2) + a_10 * std::pow(x, 1) + a_00 * std::pow(x, 0);
    };

    std::vector<float> phi_cube;
    for (float x : xs) {
        phi_cube.push_back(phi_of_x_cube(x));
    }

    using namespace sciplot;

    Plot2D plot;

    plot.xlabel("x");
    plot.ylabel("y");
    plot.grid().show();

    plot.legend()
            .atOutsideBottom()
            .displayHorizontal()
            .displayExpandWidthBy(2);

    plot.drawCurve(xs, ys).label("y = sf").lineColor("#5B0888");
    plot.drawCurve(xs, phi_lin).label("y = " + std::to_string(a) + "x + " + std::to_string(b)).lineColor("#7D7C7C");
    plot.drawCurve(xs, phi_quad).label(std::to_string(a_0) + std::to_string(a_1) + "x + " + std::to_string(a_2) + "x^2").lineColor("#26577C");
    plot.drawCurve(xs, phi_cube).label(std::to_string(a_00) + std::to_string(a_10) + "x + " + std::to_string(a_20) + "x^2 + "
                                       + std::to_string(a_30) + "x^3").lineColor("#EDB7ED");

    Figure fig = {{plot}};
    Canvas canvas = {{fig}};

    canvas.size(800, 600);

    canvas.show();
}

void plotIfXNeg(std::vector<float> &xs, std::vector<float> &ys) {
    //lineal, quad, cube, exp
    Coefficients cf = approx_lineal(xs, ys);
    float a = cf.first;
    float b = cf.second;

    auto phi_of_x_lin = [a, b](float x) -> float {
        return a * x + b;
    };

    std::vector<float> phi_lin;
    for (float x : xs) {
        phi_lin.push_back(phi_of_x_lin(x));
    }

    std::vector<float> acf = quadratic_approximation(xs, ys);
    float a_0 = acf[0];
    float a_1 = acf[1];
    float a_2 = acf[2];

    auto phi_of_x_quad = [a_0, a_1, a_2](float x) -> float {
        return a_2 * std::pow(x, 2) + a_1 * std::pow(x, 1) + a_0 * std::pow(x, 0);
    };

    std::vector<float> phi_quad;
    for (float x : xs) {
        phi_quad.push_back(phi_of_x_quad(x));
    }

    acf = cube_approximation(xs, ys);
    float a_00 = acf[0];
    float a_10 = acf[1];
    float a_20 = acf[2];
    float a_30 = acf[3];

    auto phi_of_x_cube = [a_00, a_10, a_20, a_30](float x) -> float {
        return a_30 * std::pow(x, 3) + a_20 * std::pow(x, 2) + a_10 * std::pow(x, 1) + a_00 * std::pow(x, 0);
    };

    std::vector<float> phi_cube;
    for (float x : xs) {
        phi_cube.push_back(phi_of_x_cube(x));
    }

    cf = approx_exponential(xs, ys);
    float ae = cf.first;
    float be = cf.second;

    auto phi_of_x_exp = [ae, be](float x) -> float {
        return ae * std::exp(be * x);
    };

    std::vector<float> phi_e;
    for (float x : xs) {
        phi_e.push_back(phi_of_x_exp(x));
    }

    using namespace sciplot;

    Plot2D plot;

    plot.xlabel("x");
    plot.ylabel("y");
    plot.grid().show();

    plot.legend()
            .atOutsideBottom()
            .displayHorizontal()
            .displayExpandWidthBy(2);

    plot.drawCurve(xs, ys).label("y = sf").lineColor("#5B0888");
    plot.drawCurve(xs, phi_lin).label("y = " + std::to_string(a) + "x + " + std::to_string(b)).lineColor("#7D7C7C");
    plot.drawCurve(xs, phi_e).label("y = " + std::to_string(ae) + "e^(" + std::to_string(b) + "x)").lineColor("#E55604");
    plot.drawCurve(xs, phi_quad).label(std::to_string(a_0) + std::to_string(a_1) + "x + " + std::to_string(a_2) + "x^2").lineColor("#26577C");
    plot.drawCurve(xs, phi_cube).label(std::to_string(a_00) + std::to_string(a_10) + "x + " + std::to_string(a_20) + "x^2 + "
                                       + std::to_string(a_30) + "x^3").lineColor("#EDB7ED");

    Figure fig = {{plot}};
    Canvas canvas = {{fig}};

    canvas.size(800, 600);

    canvas.show();
}

void plotIfYNeg(std::vector<float> &xs, std::vector<float> &ys) {
    //lin, quad, cube, log
    Coefficients cf = approx_lineal(xs, ys);
    float a = cf.first;
    float b = cf.second;

    auto phi_of_x_lin = [a, b](float x) -> float {
        return a * x + b;
    };

    std::vector<float> phi_lin;
    for (float x : xs) {
        phi_lin.push_back(phi_of_x_lin(x));
    }

    std::vector<float> acf = quadratic_approximation(xs, ys);
    float a_0 = acf[0];
    float a_1 = acf[1];
    float a_2 = acf[2];

    auto phi_of_x_quad = [a_0, a_1, a_2](float x) -> float {
        return a_2 * std::pow(x, 2) + a_1 * std::pow(x, 1) + a_0 * std::pow(x, 0);
    };

    std::vector<float> phi_quad;
    for (float x : xs) {
        phi_quad.push_back(phi_of_x_quad(x));
    }

    acf = cube_approximation(xs, ys);
    float a_00 = acf[0];
    float a_10 = acf[1];
    float a_20 = acf[2];
    float a_30 = acf[3];

    auto phi_of_x_cube = [a_00, a_10, a_20, a_30](float x) -> float {
        return a_30 * std::pow(x, 3) + a_20 * std::pow(x, 2) + a_10 * std::pow(x, 1) + a_00 * std::pow(x, 0);
    };

    std::vector<float> phi_cube;
    for (float x : xs) {
        phi_cube.push_back(phi_of_x_cube(x));
    }

    using namespace sciplot;

    Plot2D plot;

    plot.xlabel("x");
    plot.ylabel("y");
    plot.grid().show();

    plot.legend()
            .atOutsideBottom()
            .displayHorizontal()
            .displayExpandWidthBy(2);

    plot.drawCurve(xs, ys).label("y = sf").lineColor("#5B0888");
    plot.drawCurve(xs, phi_lin).label("y = " + std::to_string(a) + "x + " + std::to_string(b)).lineColor("#7D7C7C");
    plot.drawCurve(xs, phi_quad).label(std::to_string(a_0) + std::to_string(a_1) + "x + " + std::to_string(a_2) + "x^2").lineColor("#26577C");
    plot.drawCurve(xs, phi_cube).label(std::to_string(a_00) + std::to_string(a_10) + "x + " + std::to_string(a_20) + "x^2 + "
                                       + std::to_string(a_30) + "x^3").lineColor("#EDB7ED");

    Figure fig = {{plot}};
    Canvas canvas = {{fig}};

    canvas.size(800, 600);

    canvas.show();
}