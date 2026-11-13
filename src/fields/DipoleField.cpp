#include "aurora/fields/DipoleField.hpp"

namespace aurora {
DipoleField::DipoleField(double equatorial_T, Vec3 axis, double radius, double minimum)
    : equatorial_T_(equatorial_T), radius_m_(radius), minimum_m_(minimum), axis_(normalized(axis)) {
    if (!std::isfinite(equatorial_T) || equatorial_T<=0 || !std::isfinite(radius) || radius<=0
        || !std::isfinite(minimum) || minimum<=0 || minimum>radius)
        throw std::invalid_argument("Dipole: require B_eq>0 and 0<minimum_radius<=reference_radius");
}
double DipoleField::checked_radius(Vec3 p) const {
    const double r=norm(p);
    if (!finite(p) || !std::isfinite(r) || r<minimum_m_)
        throw std::domain_error("Position is outside the dipole field domain");
    return r;
}
Vec3 DipoleField::at(Vec3 p, double) const {
    const double r=checked_radius(p);
    const Vec3 n=p/r;
    const double ratio=radius_m_/r;
    const Vec3 b=(n*(3.0*dot(axis_,n))-axis_)*(equatorial_T_*ratio*ratio*ratio);
    if (!finite(b) || !(norm(b)>0)) throw std::overflow_error("Dipole field is outside double range");
    return b;
}
Vec3 DipoleField::magnitude_gradient(Vec3 p) const {
    const double r=checked_radius(p);
    const Vec3 n=p/r;
    const double u=dot(axis_,n);
    return (-3.0*n+(axis_-n*u)*(3.0*u/(1.0+3.0*u*u)))*(norm(at(p,0))/r);
}
} // namespace aurora
