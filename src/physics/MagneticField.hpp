#ifndef AURORA_MAGNETICFIELD_HPP
#define AURORA_MAGNETICFIELD_HPP


#include "../math/Vector3.hpp"


class MagneticField {
public:
    virtual ~MagneticField() = default;
    virtual Vector3 evaluate(const Vector3& r) const = 0;
};


#endif //AURORA_MAGNETICFIELD_HPP