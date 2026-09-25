#ifndef AURORA_VECTOR3_HPP
#define AURORA_VECTOR3_HPP


#include <cmath>
#include <ostream>


class Vector3 {
public:
    double x, y, z;

    Vector3(double x = 0, double y = 0, double z = 0) : x(x), y(y), z(z) {};

    Vector3 operator+(const Vector3 &v) const {
        return Vector3(x + v.x, y + v.y, z + v.z);
    }

    Vector3 operator-(const Vector3 &v) const {
        return Vector3(x - v.x, y - v.y, z - v.z);
    }

    Vector3 operator*(double a) const {
        return Vector3(x * a, y * a, z * a);
    }

    Vector3 operator/(double a) const {
        return Vector3(x / a, y / a, z / a);
    }

    double dot(const Vector3 &v) const {
        return x * v.x + y * v.y + z * v.z;
    }

    Vector3 cross(const Vector3 &v) const {
        return Vector3(y * v.z - z * v.y, z * v.x - x * v.z, x * v.y - y * v.x);
    }

    double norm2() const {
        return x * x + y * y + z * z;
    }

    double norm() const {
        return std::sqrt(norm2());
    }

    bool isFinite() const {
        return std::isfinite(x) && std::isfinite(y) && std::isfinite(z);
    }

    Vector3 normalize() {
        return *this / norm();
    }

    friend std::ostream& operator<<(std::ostream& os, const Vector3& v) {
        return os <<'(' << v.x << ", " << v.y << ", " << v.z << ')';
    }
};


#endif //AURORA_VECTOR3_HPP
