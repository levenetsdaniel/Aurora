#ifndef AURORA_DIPOLEFIELD_HPP
#define AURORA_DIPOLEFIELD_HPP


#include "MagneticField.hpp"
#include "Constants.hpp"
#include "../math/Vector3.hpp"
#include <cmath>
#include <stdexcept>


class DipoleField final : public MagneticField {
private:
    double b0_T_;
    double earthRadius_m_;

public:
    DipoleField(double b0_T = C::b0_T, double earthRadius_m = C::earthRadius_m) :
    b0_T_(b0_T), earthRadius_m_(earthRadius_m) {
        if (b0_T_ <= 0 || earthRadius_m_ <= 0)
            throw std::invalid_argument("Invalid dipole parameters");
    }

    Vector3 evaluate(const Vector3& r) const override {
        const double r2 = r.norm2();
        if (r2 == 0)
            throw std::invalid_argument("Dipole field is undefined at origin");

        const double radius = std::sqrt(r2);
        const Vector3 zAxis(0, 0, 1);
        const double z = r.dot(zAxis);
        const double scale = b0_T_ * std::pow(earthRadius_m_ / radius, 3);

        return (r * (3.0 * z / r2) - zAxis) * scale;
    }
};


#endif //AURORA_DIPOLEFIELD_HPP