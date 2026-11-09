#pragma once

#include "aurora/math/Vec3.hpp"

namespace aurora {
class MagneticField {
public:
    virtual ~MagneticField() = default;
    virtual Vec3 at(Vec3 position_m, double time_s) const = 0;
};

class UniformField final : public MagneticField {
public:
    explicit UniformField(Vec3 field_T) : field_T_(field_T) {
        if (!finite(field_T) || !std::isfinite(norm(field_T)))
            throw std::invalid_argument("Magnetic field must be finite");
    }
    Vec3 at(Vec3, double) const override { return field_T_; }
    Vec3 value_T() const { return field_T_; }
private:
    Vec3 field_T_;
};

class ZeroField final : public MagneticField {
public:
    Vec3 at(Vec3, double) const override { return {}; }
};
} // namespace aurora
