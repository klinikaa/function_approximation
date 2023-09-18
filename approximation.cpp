#include "approximation.h"

class LinearApproximationException : public std::exception {
    public:
        [[nodiscard]] const char* what() const noexcept override {
            return "There is no minimum for a linear approximation function!";
        }
    };

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

std::pair<float, float> approx_lineal(std::vector<float> xs, std::vector<float> ys) {
    auto squareSum =
            [](float acc, float num)-> float {
                return acc + (num * num);
            };

    float sx = std::accumulate(xs.begin(), xs.end(), 0.0f);
    float sxx = std::accumulate(xs.begin(), xs.end(), 0.0f, squareSum);

    float sy = std::accumulate(ys.begin(), ys.end(), 0.0f);
    float sxy = std::inner_product(xs.begin(), xs.end(), ys.begin(), 0.0f);

    unsigned int n = xs.size();
    float delta_0 = sxx * n - sx * sx;
    float delta_1 = sxy * n - sx * sy;
    float delta_2 = sxx * sy - sx * sxy;

    float a = delta_1 / delta_0;
    float b = delta_2 / delta_0;

    /* checking the necessary condition the existence of a minimum for the function S */
    float sum_1 = 0;
    for (size_t i = 0; i < n; i++) {
        sum_1 += (a * xs[i] + b - ys[i]) * xs[i];
    }

    float sum_2 = 0;
    for (size_t i = 0; i < n; i++) {
        sum_2 += (a * xs[i] + b - ys[i]);
    }

    sum_1 *= 2;
    sum_2 *= 2;

    /*
     * A function is unbounded: If a function has no lower bound and continues to decrease indefinitely,
     * then it will have no minimum.
     *
     * A function is not differentiable: If a function is not differentiable in a given domain,
     * then it may not have extremum points, including a minimum.
     *
     * A function has discontinuities: If a function has discontinuities or breaking points in a given region,
     * then it may not have a minimum.
     */

    /*
     * Due to errors arising from the representation of floating point numbers in a computer,
     * we will compare not with 0, but with epsilon.
     *
     * Here, epsilon is the acceptable tolerance or margin of error when comparing floating point numbers.
     * Instead of using exact equality (==) to compare sum_1 and sum_2 with 0, we use a comparison with epsilon.
     *
     * By comparing the absolute values of sum_1 and sum_2 with epsilon,
     * we can determine if they are close enough to zero.
     */

    float epsilon = 1e-3;

    if (std::abs(sum_1) < epsilon && std::abs(sum_2) < epsilon) {
        return {a, b};
    } else {
        throw LinearApproximationException();
    }
}

/* φ(x) = a * exp(b * x)
 * to apply the least squares method, the function is linearized
 * this is necessary because the least squares method assumes
 * that the relationship between variables can be described by a linear model
 *
 * we make a change of variables: A = ln(a), B = b
 * then we look for these coefficients by solving the system of equations
 * (you can reuse the linear approximation function, since we have a linear dependence)
 * at the end we do the reverse change of variables => a = exp(A)
 */
std::pair<float, float> approx_exponential(std::vector<float> &xs, std::vector<float> &ys) {
    std::vector<double> ln_ys;
    ln_ys.reserve(ys.size());

    /* data linearization */
    for (const auto& y : ys) {
        ln_ys.push_back(std::log(y));
    }

    double sum_xs = 0.0;
    double sum_ln_ys = 0.0;
    double sum_xs_ln_ys = 0.0;
    double sum_xs_squared = 0.0;

    for (size_t i = 0; i < xs.size(); i++) {
        sum_xs += xs[i];
        sum_ln_ys += ln_ys[i];
        sum_xs_ln_ys += xs[i] * ln_ys[i];
        sum_xs_squared += xs[i] * xs[i];
    }

    double n = xs.size();

    double B = (n * sum_xs_ln_ys - sum_xs * sum_ln_ys) / (n * sum_xs_squared - sum_xs * sum_xs);
    double A = std::exp((sum_ln_ys - B * sum_xs) / n);

    return std::make_pair(A, B);
}