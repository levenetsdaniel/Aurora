#ifndef AURORA_UNIFORMFIELD_HPP
#define AURORA_UNIFORMFIELD_HPP


#include "MagneticField.hpp"


class UniformField final : public MagneticField {
private:
    Vector3 field_T;

public:
    UniformField(Vector3 field_T) : field_T(field_T) {};
    Vector3 evaluate(const Vector3&) const override {
        return field_T;
    }
};


#endif //AURORA_UNIFORMFIELD_HPP
