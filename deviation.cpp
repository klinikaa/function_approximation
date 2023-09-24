#include "deviation.h"



float average(std::vector<float> &v) {
    if (v.empty()) {
        return 0;
    }

    auto const size = static_cast<float>(v.size());
    return std::accumulate(v.begin(), v.end(), 0.0f) / size;
}

/* also known as Pearson coefficient
 * used in this program only for linear approximation
 * formula: r = (Σ((x_i - x_avg) * (y_i - y_avg))) / (sqrt(Σ(x_i - x_avg)^2) * sqrt(Σ(y_i - y_avg)^2))
 *
 * The correlation coefficient is used to answer the question:
 * is there a linear relationship between the variables and how strong is it?
 *
 * r = +- 1 => strong relationship
 * r = 0 => no relationship
 * etc.
 */
float correlation_coefficient(std::vector<float> xs, std::vector<float> ys) {
    float xAverage = average(xs);
    float yAverage = average(ys);

    float numerator = 0;
    for (size_t i = 0; i < xs.size(); i++) {
        numerator += (xs[i] - xAverage) * (ys[i] - yAverage);
    }

    float denominator;

    float leftSum = 0;
    for (size_t j = 0; j < xs.size(); j++) {
        leftSum += std::pow(xs[j] - xAverage, 2);
    }

    float rightSum = 0;
    for (size_t k = 0; k < ys.size(); k++) {
        rightSum += std::pow(ys[k] - yAverage, 2);
    }

    denominator = std::sqrt(leftSum * rightSum);

    if (denominator == 0) {
        throw std::runtime_error("While processing Pearson's coefficient dominator become 0!");
    } else {
        return numerator / denominator;
    }
}

float standard_deviation(float S, size_t n) { /* S =  ∑[1, n](a*x_i + b - y_i)^2 */
    return std::sqrt(S / n);
}

float deviation_lineal(float a, float b, std::vector<float> &xs, std::vector<float> &ys) {
    /* φ(x) = ax + b
     * this is approximating function
     *
     * we substitute the values from the vector xs into it,
     * and then compare the resulting y using the least squares method
     *
     * least squares function: S = S(a, b) = ∑[1, n](ε_i^2) =
     * ∑[1, n](φ(x_i) - y_i)^2 = ∑[1, n](a*x_i + b - y_i)^2 -> min
     * */

    auto phi_of_x = [a, b](float x) -> float {
        return a * x + b;
    };

    float S = 0;
    std::vector<float> y_phi; y_phi.reserve(xs.size());

    for (float x : xs) {
        y_phi.push_back(phi_of_x(x));
    }

    for (size_t i = 0; i < y_phi.size(); i++) {
        S += std::pow(y_phi[i] - ys[i], 2);
    }

    return S;
}

float deviation_exponential(float a, float b, std::vector<float> &xs, std::vector<float> &ys) {
    /* φ(x) = a * exp(b * x)
     * this is approximating function
     *
     * we substitute the values from the vector xs into it,
     * and then compare the resulting y using the least squares method
     *
     * least squares function: S = S(a, b) = ∑[1, n](ε_i^2) =
     * ∑[1, n](φ(x_i) - y_i)^2 = ∑[1, n](a*x_i + b - y_i)^2 -> min
     * */
    auto phi_of_x = [a, b](float x) -> float {
        return a * std::exp(b * x);
    };

    float S = 0;
    std::vector<float> y_phi; y_phi.reserve(xs.size());

    for (float x : xs) {
        y_phi.push_back(phi_of_x(x));
    }

    for (size_t i = 0; i < y_phi.size(); i++) {
        S += std::pow(y_phi[i] - ys[i], 2);
    }

    return S;
}

//TODO: add comments
float deviation_power(float a, float b, std::vector<float> &xs, std::vector<float> &ys) {
    auto phi_of_x = [a, b](float x) -> float {
        return a * std::pow(x, b);
    };

    float S = 0;
    std::vector<float> y_phi; y_phi.reserve(xs.size());

    for (float x : xs) {
        y_phi.push_back(phi_of_x(x));
    }

    for (size_t i = 0; i < y_phi.size(); i++) {
        S += std::pow(y_phi[i] - ys[i], 2);
    }

    return S;
}

float deviation_log(float a, float b, std::vector<float> &xs, std::vector<float> &ys) {
    auto phi_of_x = [a, b](float x) -> float {
        return a * std::log(x) + b;
    };

    float S = 0;
    std::vector<float> y_phi; y_phi.reserve(xs.size());

    for (float x : xs) {
        y_phi.push_back(phi_of_x(x));
    }

    for (size_t i = 0; i < y_phi.size(); i++) {
        S += std::pow(y_phi[i] - ys[i], 2);
    }

    return S;
}

float deviation_quadratic(float a_0, float a_1, float a_2, std::vector<float> &xs, std::vector<float> &ys) {
    auto phi_of_x = [a_0, a_1, a_2](float x) -> float {
        return a_2 * std::pow(x, 2) + a_1 * std::pow(x, 1) + a_0 * std::pow(x, 0);
    };

    float S = 0;
    std::vector<float> y_phi; y_phi.reserve(xs.size());

    for (float x : xs) {
        y_phi.push_back(phi_of_x(x));
    }

    for (size_t i = 0; i < y_phi.size(); i++) {
        S += std::pow(y_phi[i] - ys[i], 2);
    }

    return S;
}

float deviation_qube(float a_0, float a_1, float a_2, float a_3, std::vector<float> &xs, std::vector<float> &ys) {
    auto phi_of_x = [a_0, a_1, a_2, a_3](float x) -> float {
        return a_3 * std::pow(x, 3) + a_2 * std::pow(x, 2) + a_1 * std::pow(x, 1) + a_0 * std::pow(x, 0);
    };

    float S = 0;
    std::vector<float> y_phi; y_phi.reserve(xs.size());

    for (float x : xs) {
        y_phi.push_back(phi_of_x(x));
    }

    for (size_t i = 0; i < y_phi.size(); i++) {
        S += std::pow(y_phi[i] - ys[i], 2);
    }

    return S;
}