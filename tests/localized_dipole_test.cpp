#include "../src/physics/LocalizedDipoleField.hpp"

#include <cmath>
#include <iostream>
#include <stdexcept>

namespace {

void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

double divergence(const MagneticField& field, const Vector3& position, double step_m) {
    const Vector3 x(step_m, 0, 0);
    const Vector3 y(0, step_m, 0);
    const Vector3 z(0, 0, step_m);
    return ((field.evaluate(position + x).x - field.evaluate(position - x).x) +
            (field.evaluate(position + y).y - field.evaluate(position - y).y) +
            (field.evaluate(position + z).z - field.evaluate(position - z).z)) /
           (2.0 * step_m);
}

} // namespace

int main() {
    try {
        const DipoleField background;
        const LocalizedDipoleField zero(0.32, 500000.0, 0.0);
        const LocalizedDipoleField anomaly(0.32, 500000.0, 1e-6);
        const Vector3 radial(std::sin(0.32), 0, std::cos(0.32));
        const Vector3 surface = radial * C::earthRadius_m;
        const Vector3 source = radial * (C::earthRadius_m + 450000.0);
        require((zero.evaluate(source) - background.evaluate(source)).norm() == 0.0,
                "Zero anomaly must recover the background field exactly");
        require(std::abs(anomaly.anomaly(surface).dot(radial) - 1e-6) < 1e-20,
                "Surface axial field does not match the requested amplitude");
        const double axialExpected = 1e-6 * std::pow(500000.0 / 950000.0, 3);
        require(std::abs(anomaly.anomaly(source).dot(radial) - axialExpected) < 1e-20,
                "Anomaly does not decay as an axial dipole");
        for (const Vector3 point : {source, Vector3(1.5e6, 3.0e5, 6.6e6),
                                    Vector3(-1.2e6, 2.1e6, 6.3e6)}) {
            const double numerical = divergence(anomaly, point, 10.0);
            const double scale = anomaly.evaluate(point).norm() / point.norm();
            require(std::abs(numerical) / scale < 1e-7,
                    "Numerical div B is not approximately zero outside the source");
        }
        bool rejected = false;
        try { LocalizedDipoleField invalid(0.32, 0.0, 1e-6); }
        catch (const std::invalid_argument&) { rejected = true; }
        require(rejected, "Invalid source depth was accepted");
        std::cout << "Localized dipole field passed\n";
    } catch (const std::exception& error) {
        std::cerr << "Localized dipole field failed: " << error.what() << '\n';
        return 1;
    }
}
