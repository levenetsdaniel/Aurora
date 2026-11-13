#include "aurora/simulation/EventLocator.hpp"
#include "aurora/physics/ParticleProperties.hpp"
#include "aurora/numerics/RelativisticBorisPusher.hpp"
#include <algorithm>
#include <functional>

namespace aurora {
namespace {
RelativisticBorisPusher::State advance_raw(
    const ParticleState& start, const MagneticField& field, double dt) {
    RelativisticBorisPusher pusher(dt);
    auto s = pusher.initialize(start, field);
    pusher.advance(s, field);
    return s;
}

LocatedEvent locate(const ParticleState& start, const MagneticField& field, double dt, double lo,
    double hi, double tolerance, const std::function<double(double)>& residual) {
    if (!(dt > 0) || !std::isfinite(dt) || !(tolerance > 0) || !std::isfinite(tolerance) ||
        lo < 0 || !(lo < hi) || hi > 1 || !std::isfinite(lo) || !std::isfinite(hi))
        throw std::invalid_argument("Invalid event bracket or tolerance");

    double a = residual(lo), b = residual(hi);
    if (!std::isfinite(a) || !std::isfinite(b) ||
        (a != 0 && b != 0 && std::signbit(a) == std::signbit(b)))
        throw std::invalid_argument("Event root is not bracketed");

    if (a == 0)
        return {fractional_step(start, field, lo * dt), 0};
    if (b == 0)
        return {fractional_step(start, field, hi * dt), 0};

    for (int iteration = 0; (hi - lo) * dt > tolerance; ++iteration) {
        const double mid = lo + (hi - lo) * 0.5;
        if (iteration == 100 || mid == lo || mid == hi)
            throw std::runtime_error("Event tolerance is below numeric resolution");

        const double value = residual(mid);
        if (!std::isfinite(value))
            throw std::runtime_error("Non-finite event residual");
        if (value == 0)
            return {fractional_step(start, field, mid * dt), 0};

        if (std::signbit(value) == std::signbit(a)) {
            lo = mid;
            a = value;
        } else
            hi = mid;
    }

    return {fractional_step(start, field, (lo + (hi - lo) * 0.5) * dt), (hi - lo) * dt};
}
}

ParticleState fractional_step(
    const ParticleState& start, const MagneticField& field, double interval) {
    if (interval == 0)
        return start;

    RelativisticBorisPusher pusher(interval);
    return pusher.synchronize(advance_raw(start, field, interval), field);
}

LocatedEvent locate_radius_event(const ParticleState& start, const MagneticField& field, double dt,
    double radius, double lo, double hi, double tolerance) {
    if (!(radius > 0) || !std::isfinite(radius))
        throw std::invalid_argument("Invalid event radius");

    return locate(start, field, dt, lo, hi, tolerance, [&](double f) {
        const Vec3 position =
            f == 0 ? start.position_m : advance_raw(start, field, f * dt).position_m();
        // F(f) = |r(t_start + f*dt)| - R_boundary; на сфере F = 0.
        return norm(position) - radius;
    });
}

LocatedEvent locate_mirror_event(
    const ParticleState& start, const MagneticField& field, double dt, double tolerance) {
    return locate(start, field, dt, 0, 1, tolerance, [&](double f) {
        // F(f) = p_parallel/|p|; при отражении продольный импульс проходит через ноль.
        return particle_properties(fractional_step(start, field, f * dt), field).parallel_fraction;
    });
}

double inner_crossing_bracket(Vec3 start, Vec3 end, double radius) {
    if (norm(end) < radius)
        return 1;

    // d = r_end - r_start; прямой отрезок шага: r(f) = r_start + f*d.
    const Vec3 delta = end - start;
    const double length = norm(delta);
    if (!(length > 0))
        return 0;

    const double closest = std::clamp(-dot(start, delta / length) / length, 0.0, 1.0);
    return norm(start + delta * closest) < radius ? closest : 0;
}
} // namespace aurora
