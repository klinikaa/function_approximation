#include "approximation.h"

class LinearApproximationException : public std::exception {
    public:
        [[nodiscard]] const char* what() const noexcept override {
            return "There is no minimum for a linear approximation function!";
        }
    };

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
