#include <vector>
#include <numeric>
#include <valarray>

class LinearApproximationException : public std::exception {
    public:
        [[nodiscard]] const char* what() const noexcept override {
            return "There is no minimum for a linear approximation function!";
        }
    };

struct result {
    std::vector<float> x;
    std::vector<float> y;
};

std::pair<float, float> approx_lineal(std::vector<float> xs, std::vector<float> ys) {
    auto squareSum {
            [](int acc, int num)-> float {
                return acc + (num * num);
            }
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

    // let's check the necessary condition the existence of a minimum for the function S
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

    if (sum_1 == 0 && sum_2 == 0) {
        return {a, b};
    } else {
        throw LinearApproximationException();
    }
}
