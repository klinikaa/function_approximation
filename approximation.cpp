#include "approximation.h"

class LinearApproximationException : public std::exception {
    public:
        [[nodiscard]] const char* what() const noexcept override {
            return "There is no minimum for a linear approximation function!";
        }
    };

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
    std::vector<float> ln_ys;
    ln_ys.reserve(ys.size());

    /* data linearization */
    for (const auto& y : ys) {
        ln_ys.push_back(std::log(y));
    }

    float sum_xs = 0.0;
    float sum_ln_ys = 0.0;
    float sum_xs_ln_ys = 0.0;
    float sum_xs_squared = 0.0;

    for (size_t i = 0; i < xs.size(); i++) {
        sum_xs += xs[i];
        sum_ln_ys += ln_ys[i];
        sum_xs_ln_ys += xs[i] * ln_ys[i];
        sum_xs_squared += xs[i] * xs[i];
    }

    float n = xs.size();

    float B = (n * sum_xs_ln_ys - sum_xs * sum_ln_ys) / (n * sum_xs_squared - sum_xs * sum_xs);
    float A = std::exp((sum_ln_ys - B * sum_xs) / n);

    return {A, B};
}

std::pair<float, float> approx_power(std::vector<float> &xs, std::vector<float> &ys) {
    std::vector<float> ln_ys;
    ln_ys.reserve(ys.size());

    std::vector<float> ln_xs;
    ln_xs.reserve(xs.size());

    size_t n = xs.size();
    for (size_t i = 0; i < n; i++) {
        ln_ys.push_back(std::log(ys[i]));
        ln_xs.push_back(std::log(xs[i]));
    }

    float sum_xs = 0.0;
    float sum_ln_ys = 0.0;
    float sum_xs_ln_ys = 0.0;
    float sum_xs_squared = 0.0;

    for (size_t i = 0; i < n; i++) {
        sum_xs += ln_xs[i];
        sum_ln_ys += ln_ys[i];
        sum_xs_ln_ys += ln_xs[i] * ln_ys[i];
        sum_xs_squared += ln_xs[i] * ln_xs[i];
    }

    float B = (n * sum_xs_ln_ys - sum_xs * sum_ln_ys) / (n * sum_xs_squared - sum_xs * sum_xs);
    float A = std::exp((sum_ln_ys - B * sum_xs) / n);

    return {A, B};
}

//FIXME
std::pair<float, float> approx_log(std::vector<float> &xs, std::vector<float> &ys) {
    std::vector<float> ln_xs;
    ln_xs.reserve(xs.size());

    size_t n = xs.size();
    for (size_t i = 0; i < n; i++) {
        ln_xs.push_back(std::log(xs[i]));
    }

    float sum_xs = 0.0;
    float sum_ln_ys = 0.0;
    float sum_xs_ln_ys = 0.0;
    float sum_xs_squared = 0.0;

    for (size_t i = 0; i < n; i++) {
        sum_xs += ln_xs[i];
        sum_ln_ys += ys[i];
        sum_xs_ln_ys += ln_xs[i] * ys[i];
        sum_xs_squared += ln_xs[i] * ln_xs[i];
    }

    float B = (n * sum_xs_ln_ys - sum_xs * sum_ln_ys) / (n * sum_xs_squared - sum_xs * sum_xs);
    float A = (sum_ln_ys - B * sum_xs) / n;

    return {B, A};
}

/*
 *
 *
 * return value contains 3 float coefficients
 */
//TODO: checking the necessary condition the existence of a minimum for the function S
std::vector<float> quadratic_approximation(std::vector<float> xs, std::vector<float> ys) {
    auto squareSum =
            [](float acc, float num) -> float {
                return acc + (num * num);
            };

    auto cubeSum = [](float acc, float num) -> float {
        return acc + (num * num * num);
    };

    auto fourSum = [](float acc, float num) -> float {
        return acc + (num * num * num * num);
    };

    float sx = std::accumulate(xs.begin(), xs.end(), 0.0f);
    float sxx = std::accumulate(xs.begin(), xs.end(), 0.0f, squareSum);
    float sxxx = std::accumulate(xs.begin(), xs.end(), 0.0f, cubeSum);
    float sxxxx = std::accumulate(xs.begin(), xs.end(), 0.0f, fourSum);

    float sy = std::accumulate(ys.begin(), ys.end(), 0.0f);
    float sxy = std::inner_product(xs.begin(), xs.end(), ys.begin(), 0.0f);
    float sxxy = 0.0f;
    for (size_t i = 0; i < xs.size(); i++) {
        sxxy += xs[i] * xs[i] * ys[i];
    }

    size_t n = xs.size();

    Eigen::Matrix3f A;
    A << n, sx, sxx,
            sx, sxx, sxxx,
            sxx, sxxx, sxxxx;

    Eigen::Vector3f B;
    B << sy, sxy, sxxy;

    float detA = A.determinant();

    if (detA == 0) {
        throw std::runtime_error("The system of equations has no unique solution!");
    } else {
        std::cout << "DET = " << detA << std::endl;
    }

    Eigen::Vector3f a;
    a = A.colPivHouseholderQr().solve(B);

    return {
            a(0), a(1), a(2)
    };
}

std::vector<float> cube_approximation(std::vector<float> xs, std::vector<float> ys) {
    auto squareSum =
            [](float acc, float num) -> float {
                return acc + (num * num);
            };

    auto cubeSum = [](float acc, float num) -> float {
        return acc + (num * num * num);
    };

    auto fourSum = [](float acc, float num) -> float {
        return acc + (num * num * num * num);
    };

    auto fiveSum = [](float acc, float num) -> float {
        return acc + (num * num * num * num * num);
    };

    float sx = std::accumulate(xs.begin(), xs.end(), 0.0f);
    float sxx = std::accumulate(xs.begin(), xs.end(), 0.0f, squareSum);
    float sxxx = std::accumulate(xs.begin(), xs.end(), 0.0f, cubeSum);
    float sxxxx = std::accumulate(xs.begin(), xs.end(), 0.0f, fourSum);
    float sxxxxx = std::accumulate(xs.begin(), xs.end(), 0.0f, fiveSum);

    float sy = std::accumulate(ys.begin(), ys.end(), 0.0f);
    float sxy = std::inner_product(xs.begin(), xs.end(), ys.begin(), 0.0f);
    float sxxy = 0.0f;
    float sxxxy = 0.0f;
    for (size_t i = 0; i < xs.size(); i++) {
        sxxy += xs[i] * xs[i] * ys[i];
        sxxxy += xs[i] * xs[i] * xs[i] * ys[i];
    }

    size_t n = xs.size();

    Eigen::Matrix4f A;
    A << n, sx, sxx, sxxx,
            sx, sxx, sxxx, sxxxx,
            sxx, sxxx, sxxxx, sxxxxx,
            sxxx, sxxxx, sxxxxx, sxxxy;

    Eigen::Vector4f B;
    B << sy, sxy, sxxy, sxxxy;

    float detA = A.determinant();

    if (detA == 0) {
        throw std::runtime_error("The system of equations has no unique solution!");
    }

    Eigen::Vector4f a;
    a = A.colPivHouseholderQr().solve(B);

    return {
            a(0), a(1), a(2), a(3)
    };
}