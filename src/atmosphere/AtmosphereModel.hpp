#ifndef AURORA_ATMOSPHEREMODEL_HPP
#define AURORA_ATMOSPHEREMODEL_HPP

#include <algorithm>
#include <cmath>
#include <numbers>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

struct AtmosphereSample {
    double altitude_m;
    double nO_m3;
    double nO2_m3;
    double nN2_m3;
    double temperature_K;
};

class AtmosphereModel {
private:
    std::vector<AtmosphereSample> samples;
    std::string sourceMetadata;

    static double interpolateDensity(double lower, double upper, double fraction) {
        if (lower == 0.0 || upper == 0.0) {
            return lower + fraction * (upper - lower);
        }
        return std::exp(std::log(lower) + fraction * (std::log(upper) - std::log(lower)));
    }

public:
    AtmosphereModel(std::vector<AtmosphereSample> samples, std::string sourceMetadata)
        : samples(std::move(samples)), sourceMetadata(std::move(sourceMetadata)) {
        if (this->samples.size() < 2 || this->sourceMetadata.empty()) {
            throw std::invalid_argument("AtmosphereModel: at least two samples and source metadata are required");
        }

        double previousAltitude_m = -1.0;
        for (const AtmosphereSample& sample : this->samples) {
            if (!std::isfinite(sample.altitude_m) || sample.altitude_m < 0.0 ||
                sample.altitude_m <= previousAltitude_m ||
                !std::isfinite(sample.nO_m3) || sample.nO_m3 < 0.0 ||
                !std::isfinite(sample.nO2_m3) || sample.nO2_m3 < 0.0 ||
                !std::isfinite(sample.nN2_m3) || sample.nN2_m3 < 0.0 ||
                !std::isfinite(sample.temperature_K) || sample.temperature_K <= 0.0) {
                throw std::invalid_argument("AtmosphereModel: invalid or unordered atmosphere table");
            }
            previousAltitude_m = sample.altitude_m;
        }
    }

    AtmosphereSample sample(double altitude_m, double latitude_rad,
                            double longitude_rad, double time_s) const {
        if (!std::isfinite(altitude_m) || !std::isfinite(latitude_rad) ||
            !std::isfinite(longitude_rad) || !std::isfinite(time_s) ||
            std::abs(latitude_rad) > std::numbers::pi / 2.0 || time_s < 0.0) {
            throw std::invalid_argument("AtmosphereModel: invalid coordinates or time");
        }
        if (altitude_m < minAltitude_m() || altitude_m > maxAltitude_m()) {
            throw std::out_of_range("AtmosphereModel: altitude is outside the table");
        }

        const auto upper = std::lower_bound(
            samples.begin(), samples.end(), altitude_m,
            [](const AtmosphereSample& entry, double altitude) {
                return entry.altitude_m < altitude;
            });
        if (upper->altitude_m == altitude_m) {
            return *upper;
        }

        const AtmosphereSample& lower = *(upper - 1);
        const double fraction = (altitude_m - lower.altitude_m) /
                                (upper->altitude_m - lower.altitude_m);
        return {
            altitude_m,
            interpolateDensity(lower.nO_m3, upper->nO_m3, fraction),
            interpolateDensity(lower.nO2_m3, upper->nO2_m3, fraction),
            interpolateDensity(lower.nN2_m3, upper->nN2_m3, fraction),
            lower.temperature_K + fraction * (upper->temperature_K - lower.temperature_K)
        };
    }

    double minAltitude_m() const {
        return samples.front().altitude_m;
    }

    double maxAltitude_m() const {
        return samples.back().altitude_m;
    }

    const std::string& source() const {
        return sourceMetadata;
    }
};

#endif // AURORA_ATMOSPHEREMODEL_HPP
