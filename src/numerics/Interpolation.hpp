#ifndef AURORA_INTERPOLATION_HPP
#define AURORA_INTERPOLATION_HPP

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>
#include <vector>

struct TabulatedPoint {
    double energy_eV{};
    double sigma_m2{};
};

// Linear interpolation in eV/m²; values outside the table are never extrapolated.
class TabulatedFunction {
private:
    std::vector<TabulatedPoint> points;

public:
    explicit TabulatedFunction(std::vector<TabulatedPoint> points)
        : points(std::move(points)) {
        if (this->points.size() < 2) {
            throw std::invalid_argument("TabulatedFunction: at least two points are required");
        }
        double previousEnergy_eV = -1.0;
        for (const TabulatedPoint& point : this->points) {
            if (!std::isfinite(point.energy_eV) || point.energy_eV < 0.0 ||
                point.energy_eV <= previousEnergy_eV ||
                !std::isfinite(point.sigma_m2) || point.sigma_m2 < 0.0) {
                throw std::invalid_argument("TabulatedFunction: invalid or unordered points");
            }
            previousEnergy_eV = point.energy_eV;
        }
    }

    double evaluate(double energy_eV) const {
        if (!std::isfinite(energy_eV) || energy_eV < 0.0) {
            throw std::invalid_argument("TabulatedFunction: energy must be finite and nonnegative");
        }
        if (energy_eV < minEnergy_eV() || energy_eV > maxEnergy_eV()) {
            throw std::out_of_range("TabulatedFunction: energy is outside the table");
        }

        const auto upper = std::lower_bound(
            points.begin(), points.end(), energy_eV,
            [](const TabulatedPoint& point, double energy) {
                return point.energy_eV < energy;
            });
        if (upper->energy_eV == energy_eV) return upper->sigma_m2;

        const TabulatedPoint& lower = *(upper - 1);
        const double fraction = (energy_eV - lower.energy_eV) /
                                (upper->energy_eV - lower.energy_eV);
        return lower.sigma_m2 + fraction * (upper->sigma_m2 - lower.sigma_m2);
    }

    double minEnergy_eV() const {
        return points.front().energy_eV;
    }

    double maxEnergy_eV() const {
        return points.back().energy_eV;
    }

    double firstValue_m2() const {
        return points.front().sigma_m2;
    }

    double lastValue_m2() const {
        return points.back().sigma_m2;
    }
};

#endif // AURORA_INTERPOLATION_HPP
