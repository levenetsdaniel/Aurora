#pragma once

#include <cmath>
#include <stdexcept>

namespace aurora {
struct Vec3 {
    double x{}, y{}, z{};

    constexpr Vec3 operator+(Vec3 b) const {
        return {x + b.x, y + b.y, z + b.z};
    }

    constexpr Vec3 operator-(Vec3 b) const {
        return {x - b.x, y - b.y, z - b.z};
    }

    constexpr Vec3 operator-() const {
        return {-x, -y, -z};
    }

    constexpr Vec3 operator*(double s) const {
        return {x * s, y * s, z * s};
    }

    constexpr Vec3 operator/(double s) const {
        return {x / s, y / s, z / s};
    }

    constexpr Vec3& operator+=(Vec3 b) {
        *this = *this + b;
        return *this;
    }
};

constexpr Vec3 operator*(double s, Vec3 a) {
    return a * s;
}

constexpr double dot(Vec3 a, Vec3 b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

constexpr Vec3 cross(Vec3 a, Vec3 b) {
    return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}

inline double norm(Vec3 a) {
    return std::hypot(a.x, a.y, a.z);
}

inline bool finite(Vec3 a) {
    return std::isfinite(a.x) && std::isfinite(a.y) && std::isfinite(a.z);
}

inline Vec3 normalized(Vec3 a) {
    const double length = norm(a);
    if (!finite(a) || !(length > 0.0) || !std::isfinite(length))
        throw std::invalid_argument("Cannot normalize a zero or non-finite vector");
    return a / length;
}
} // namespace aurora
