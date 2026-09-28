#ifndef AURORA_REFERENCECROSSSECTIONS_HPP
#define AURORA_REFERENCECROSSSECTIONS_HPP

#include "CrossSectionDatabase.hpp"

#include <initializer_list>
#include <utility>
#include <vector>

// Laher & Gilmore, J. Phys. Chem. Ref. Data 19, 277 (1990), Tables 1 and 3.
// https://doi.org/10.1063/1.555872
// Table 3 is in 10^-18 cm^2 = 10^-22 m^2. No values above 200 eV are
// tabulated here; the database deliberately does not extrapolate them.
inline std::vector<CollisionChannel> atomicOxygenExcitationChannels() {
    constexpr double tableUnit_m2 = 1e-22;
    const char* source = "Laher & Gilmore (1990), JPCRD 19, 277, Tables 1 and 3; doi:10.1063/1.555872";

    const auto table = [](std::initializer_list<std::pair<double, double>> values) {
        std::vector<TabulatedPoint> points;
        points.reserve(values.size());
        for (const auto& [energy_eV, value] : values) {
            points.push_back({energy_eV, value * tableUnit_m2});
        }
        return TabulatedFunction(std::move(points));
    };

    return {
        {"O-3P-to-1D", "O", CollisionType::Excitation, "O(1D)", 1.96,
         table({
             {1.96, 0.0}, {2.1, 0.82}, {2.4, 5.41}, {2.7, 10.7},
             {3.0, 15.4}, {4.0, 24.9}, {4.4, 26.7}, {4.8, 27.7},
             {5.0, 28.0}, {6.0, 28.3}, {7.0, 27.4}, {8.0, 26.0},
             {9.0, 24.4}, {10.0, 22.8}, {12.0, 19.8}, {14.0, 17.1},
             {16.0, 14.8}, {18.0, 12.9}, {20.0, 11.2}, {25.0, 8.07},
             {30.0, 5.90}, {40.0, 3.34}, {45.0, 2.58}, {50.0, 2.02},
             {55.0, 1.61}, {70.0, 0.87}, {100.0, 0.33},
             {150.0, 0.11}, {200.0, 0.05}
         }), source},
        {"O-3P-to-1S", "O", CollisionType::Excitation, "O(1S)", 4.18,
         table({
             {4.18, 0.0}, {4.4, 0.54}, {4.8, 1.30}, {5.0, 1.60},
             {6.0, 2.57}, {7.0, 3.02}, {8.0, 3.22}, {9.0, 3.28},
             {10.0, 3.27}, {12.0, 3.15}, {14.0, 2.97}, {16.0, 2.79},
             {18.0, 2.61}, {20.0, 2.44}, {25.0, 2.09}, {30.0, 1.80},
             {40.0, 1.37}, {45.0, 1.21}, {50.0, 1.07}, {55.0, 0.95},
             {70.0, 0.69}, {100.0, 0.38}, {150.0, 0.16},
             {200.0, 0.08}
         }), source}
    };
}

#endif // AURORA_REFERENCECROSSSECTIONS_HPP
