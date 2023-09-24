#include <vector>
#include "util.h"

bool hasNegativeNumber(std::vector<float>& numbers) {
    for (int number : numbers) {
        if (number < 0) {
            return true;
        }
    }

    return false;
}
