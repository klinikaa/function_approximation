#include "process.h"

float process_lineal(Points &xs, Points &ys) {
    std::cout << "<lineal approximation>" << std::endl;

    Coefficients cf = approx_lineal(xs, ys);
    float a = cf.first;
    float b = cf.second;

    auto phi_of_x = [a, b](float x) -> float {
        return a * x + b;
    };

    std::vector<float> phi;
    for (float x : xs) {
        phi.push_back(phi_of_x(x));
    }

    std::string eq = std::to_string(a) + "x + " + std::to_string(b);
    const std::vector<std::string> HEADERS = {"i", "x", "y", eq, "epsilon"};
    std::vector<std::vector<std::string>> LINES;

    for (size_t i = 0; i < xs.size(); i++) {
        std::vector<std::string> LINE;

        LINE.push_back(std::to_string(i + 1));
        LINE.push_back(std::to_string(xs[i]));
        LINE.push_back(std::to_string(ys[i]));
        LINE.push_back(std::to_string(phi[i]));
        LINE.push_back(std::to_string(std::abs(phi[i] - ys[i])));
        LINES.push_back(LINE);
    }

    std::cout << "<TABLE>" << std::endl;
    printTable(HEADERS, LINES);

    std::cout << "We got a = " << a << " and b = " << b << std::endl;

    float linealDeviation = deviation_lineal(a, b, xs, ys);
    std::cout << "Deviation measure for linear approximation = " << linealDeviation << std::endl;

    size_t n = xs.size();
    float linealStandardDeviation = standard_deviation(linealDeviation, n);
    std::cout << "Standard deviation for lineal approximation (δ)= " << linealStandardDeviation << std::endl;

    float pearsonCoefficient = correlation_coefficient(xs, ys);
    std::cout << "Pearson coefficient for linear approximation (r) = " << pearsonCoefficient << std::endl;

    std::cout << "<lineal approximation> [END]" << std::endl;

    return linealStandardDeviation;
}

float process_quadratic(Points &xs, Points &ys) {
    std::cout << "<quadratic approximation>" << std::endl;

    std::vector<float> cf = quadratic_approximation(xs, ys);
    float a_0 = cf[0];
    float a_1 = cf[1];
    float a_2 = cf[2];

    auto phi_of_x = [a_0, a_1, a_2](float x) -> float {
        return a_2 * std::pow(x, 2) + a_1 * std::pow(x, 1) + a_0 * std::pow(x, 0);
    };

    std::vector<float> phi;
    for (float x : xs) {
        phi.push_back(phi_of_x(x));
    }

    std::string eq = std::to_string(a_0) + std::to_string(a_1) + "x + " + std::to_string(a_2) + "x^2";
    const std::vector<std::string> HEADERS = {"i", "x", "y", eq, "epsilon"};
    std::vector<std::vector<std::string>> LINES;

    for (size_t i = 0; i < xs.size(); i++) {
        std::vector<std::string> LINE;

        LINE.push_back(std::to_string(i + 1));
        LINE.push_back(std::to_string(xs[i]));
        LINE.push_back(std::to_string(ys[i]));
        LINE.push_back(std::to_string(phi[i]));
        LINE.push_back(std::to_string(std::abs(phi[i] - ys[i])));
        LINES.push_back(LINE);
    }

    std::cout << "<TABLE>" << std::endl;
    printTable(HEADERS, LINES);

    std::cout << "We got a_0 = " << a_0 << " and a_1 = " << a_1 << " and a_2 = " << a_2 << std::endl;

    float quadraticDeviation = deviation_quadratic(a_0, a_1, a_2 ,xs, ys);
    std::cout << "Deviation measure for quadratic approximation = " << quadraticDeviation << std::endl;

    size_t n = xs.size();
    float quadraticStandardDeviation = standard_deviation(quadraticDeviation, n);
    std::cout << "Standard deviation for quadratic approximation (δ)= " << quadraticStandardDeviation << std::endl;

    std::cout << "<quadratic approximation> [END]" << std::endl;

    return quadraticStandardDeviation;
}

float process_qube(Points &xs, Points &ys) {
    std::cout << "<qube approximation>" << std::endl;

    std::vector<float> cf = cube_approximation(xs, ys);
    float a_0 = cf[0];
    float a_1 = cf[1];
    float a_2 = cf[2];
    float a_3 = cf[3];

    auto phi_of_x = [a_0, a_1, a_2, a_3](float x) -> float {
        return a_3 * std::pow(x, 3) + a_2 * std::pow(x, 2) + a_1 * std::pow(x, 1) + a_0 * std::pow(x, 0);
    };

    std::vector<float> phi;
    for (float x : xs) {
        phi.push_back(phi_of_x(x));
    }

    std::string eq = std::to_string(a_0) + std::to_string(a_1) + "x + " + std::to_string(a_2) + "x^2"
            + std::to_string(a_3) + "x^3";
    const std::vector<std::string> HEADERS = {"i", "x", "y", eq, "epsilon"};
    std::vector<std::vector<std::string>> LINES;

    for (size_t i = 0; i < xs.size(); i++) {
        std::vector<std::string> LINE;

        LINE.push_back(std::to_string(i + 1));
        LINE.push_back(std::to_string(xs[i]));
        LINE.push_back(std::to_string(ys[i]));
        LINE.push_back(std::to_string(phi[i]));
        LINE.push_back(std::to_string(std::abs(phi[i] - ys[i])));
        LINES.push_back(LINE);
    }

    std::cout << "<TABLE>" << std::endl;
    printTable(HEADERS, LINES);

    std::cout << "We got a_0 = " << a_0 << " and a_1 = " << a_1 << " and a_2 = " << a_2 << " and a_3 = " << a_3 <<std::endl;

    float cubeDeviation = deviation_qube(a_0, a_1, a_2, a_3, xs, ys);
    std::cout << "Deviation measure for cube approximation = " << cubeDeviation << std::endl;

    size_t n = xs.size();
    float cubeStandardDeviation = standard_deviation(cubeDeviation, n);
    std::cout << "Standard deviation for cube approximation (δ)= " << cubeStandardDeviation << std::endl;

    std::cout << "<cube approximation> [END]" << std::endl;

    return cubeStandardDeviation;
}

float process_power(Points &xs, Points &ys) {
    std::cout << "<power approximation>" << std::endl;

    Coefficients cf = approx_power(xs, ys);
    float a = cf.first;
    float b = cf.second;

    auto phi_of_x = [a, b](float x) -> float {
        return a * std::pow(x, b);
    };

    std::vector<float> phi;
    for (float x : xs) {
        phi.push_back(phi_of_x(x));
    }

    std::string eq = std::to_string(a) + "x^" + std::to_string(b);
    const std::vector<std::string> HEADERS = {"i", "x", "y", eq, "epsilon"};
    std::vector<std::vector<std::string>> LINES;

    for (size_t i = 0; i < xs.size(); i++) {
        std::vector<std::string> LINE;

        LINE.push_back(std::to_string(i + 1));
        LINE.push_back(std::to_string(xs[i]));
        LINE.push_back(std::to_string(ys[i]));
        LINE.push_back(std::to_string(phi[i]));
        LINE.push_back(std::to_string(std::abs(phi[i] - ys[i])));
        LINES.push_back(LINE);
    }

    std::cout << "<TABLE>" << std::endl;
    printTable(HEADERS, LINES);

    std::cout << "We got a = " << a << " and b = " << b << std::endl;

    float powerDeviation = deviation_power(a, b, xs, ys);
    std::cout << "Deviation measure for power approximation = " << powerDeviation << std::endl;

    size_t n = xs.size();
    float powerStandardDeviation = standard_deviation(powerDeviation, n);
    std::cout << "Standard deviation for power approximation (δ)= " << powerStandardDeviation << std::endl;

    return powerStandardDeviation;
}

float process_exp(Points &xs, Points &ys) {
    std::cout << "<exp approximation>" << std::endl;

    Coefficients cf = approx_exponential(xs, ys);
    float a = cf.first;
    float b = cf.second;

    auto phi_of_x = [a, b](float x) -> float {
        return a * std::exp(b * x);
    };

    std::vector<float> phi;
    for (float x : xs) {
        phi.push_back(phi_of_x(x));
    }

    std::string eq = std::to_string(a) + "e^("  + std::to_string(b) + "x)";
    const std::vector<std::string> HEADERS = {"i", "x", "y", eq, "epsilon"};
    std::vector<std::vector<std::string>> LINES;

    for (size_t i = 0; i < xs.size(); i++) {
        std::vector<std::string> LINE;

        LINE.push_back(std::to_string(i + 1));
        LINE.push_back(std::to_string(xs[i]));
        LINE.push_back(std::to_string(ys[i]));
        LINE.push_back(std::to_string(phi[i]));
        LINE.push_back(std::to_string(std::abs(phi[i] - ys[i])));
        LINES.push_back(LINE);
    }

    std::cout << "<TABLE>" << std::endl;
    printTable(HEADERS, LINES);

    std::cout << "We got a = " << a << " and b = " << b << std::endl;

    float exponentialDeviation = deviation_exponential(a, b, xs, ys);
    std::cout << "Deviation measure for exponential approximation = " << exponentialDeviation << std::endl;

    size_t e_n = xs.size();
    float exponentialStandardDeviation = standard_deviation(exponentialDeviation, e_n);
    std::cout << "Standard deviation for exponential approximation (δ)= " << exponentialStandardDeviation << std::endl;

    return exponentialStandardDeviation;
}

float process_log(Points &xs, Points &ys) {
    std::cout << "<log approximation>" << std::endl;

    Coefficients cf = approx_log(xs, ys);
    float a = cf.first;
    float b = cf.second;

    auto phi_of_x = [a, b](float x) -> float {
        return a * std::log(x) + b;
    };

    std::vector<float> phi;
    for (float x : xs) {
        phi.push_back(phi_of_x(x));
    }

    std::string eq = std::to_string(a) + "ln(x) + " + std::to_string(b);
    const std::vector<std::string> HEADERS = {"i", "x", "y", eq, "epsilon"};
    std::vector<std::vector<std::string>> LINES;

    for (size_t i = 0; i < xs.size(); i++) {
        std::vector<std::string> LINE;

        LINE.push_back(std::to_string(i + 1));
        LINE.push_back(std::to_string(xs[i]));
        LINE.push_back(std::to_string(ys[i]));
        LINE.push_back(std::to_string(phi[i]));
        LINE.push_back(std::to_string(std::abs(phi[i] - ys[i])));
        LINES.push_back(LINE);
    }

    std::cout << "<TABLE>" << std::endl;
    printTable(HEADERS, LINES);

    std::cout << "We got a = " << a << " and b = " << b << std::endl;
    float logDeviation = deviation_log(a, b, xs, ys);
    std::cout << "Deviation measure for log approximation = " << logDeviation << std::endl;

    size_t l_n = xs.size();
    float logStandardDeviation = standard_deviation(logDeviation, l_n);
    std::cout << "Standard deviation for log approximation (δ)= " << logStandardDeviation << std::endl;
    std::cout << "log approximation [END]" << std::endl;

    return logStandardDeviation;
}