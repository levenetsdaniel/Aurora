#ifndef AURORA_RANDOM_HPP
#define AURORA_RANDOM_HPP

#include <cstdint>
#include <random>

class Random {
private:
    std::mt19937_64 engine;

public:
    explicit Random(std::uint64_t seed) : engine(seed) {}

    double uniform01() {
        // The top 53 random bits map reproducibly to [0, 1).
        return static_cast<double>(engine() >> 11) * 0x1.0p-53;
    }
};

#endif // AURORA_RANDOM_HPP
