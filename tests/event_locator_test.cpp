#include "../src/simulation/EventDetector.hpp"
#include "../src/physics/UniformField.hpp"
#include "../src/physics/Constants.hpp"
#include "../src/physics/Relativity.hpp"

#include <cmath>
#include <iostream>
#include <stdexcept>

namespace {

void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

template <typename Exception, typename Action>
void requireThrows(Action action, const char* message) {
    try {
        action();
    } catch (const Exception&) {
        return;
    }
    throw std::runtime_error(message);
}

double testMomentum() {
    return Relativity::momentumMagnitudeFromKineticEnergy_J(
        1000.0 * C::elementaryCharge_C, C::electronMass_kg);
}

void checkRadialCrossings() {
    const double radius_m = C::earthRadius_m + 150000.0;
    const double distance_m = 200.0;
    const double dt_s = 5e-5;
    const UniformField field(Vector3(1e-9, 0.0, 0.0));
    const BorisIntegrator integrator(field, -C::elementaryCharge_C, C::electronMass_kg);
    const double p = testMomentum();
    const double speed_m_s = Relativity::velocityFromMomentum(
        Vector3(p, 0.0, 0.0), C::electronMass_kg).norm();
    const double exactTime_s = distance_m / speed_m_s;

    const ParticleState inwardBefore{Vector3(radius_m + distance_m, 0.0, 0.0),
                                     Vector3(-p, 0.0, 0.0), 0.0};
    const ParticleState inwardAfter = integrator.step(inwardBefore, dt_s);
    const SphereCrossingLocator loose(radius_m, 1e-7, 100);
    const SphereCrossingLocator tight(radius_m, 1e-12, 100);
    const EventLocation coarse = loose.locate(inwardBefore, inwardAfter, integrator);
    const EventLocation fine = tight.locate(inwardBefore, inwardAfter, integrator);

    require(std::abs(fine.state.time_s - exactTime_s) <= 1e-12,
            "inward event time differs from analytic crossing");
    require(std::abs(fine.radiusResidual_m) < 1e-5,
            "inward event is not on the sphere");
    require(std::abs(fine.radiusResidual_m - (fine.state.r.norm() - radius_m)) < 1e-12,
            "reported radius residual disagrees with event state");
    require(fine.iterations > coarse.iterations,
            "tightening tolerance did not refine the event");
    require(std::abs(fine.state.time_s - exactTime_s) <=
            std::abs(coarse.state.time_s - exactTime_s) + 1e-12,
            "tightening tolerance worsened event time");

    const ParticleState outwardBefore{Vector3(radius_m - distance_m, 0.0, 0.0),
                                      Vector3(p, 0.0, 0.0), 0.0};
    const ParticleState outwardAfter = integrator.step(outwardBefore, dt_s);
    const EventLocation outward = tight.locate(outwardBefore, outwardAfter, integrator);
    require(std::abs(outward.state.time_s - exactTime_s) <= 1e-12,
            "outward event time differs from analytic crossing");
    require(std::abs(outward.radiusResidual_m) < 1e-5,
            "outward event is not on the sphere");

    requireThrows<std::invalid_argument>([&] {
        (void)tight.locate(inwardBefore, inwardBefore, integrator);
    }, "non-crossing states were accepted");
    requireThrows<std::runtime_error>([&] {
        (void)SphereCrossingLocator(radius_m, 1e-15, 1).locate(
            inwardBefore, inwardAfter, integrator);
    }, "iteration limit did not report failure to reach tolerance");
}

void checkPhysicalSubstep() {
    const double radius_m = C::earthRadius_m + 150000.0;
    const UniformField field(Vector3(0.0, 0.0, 5e-5));
    const BorisIntegrator integrator(field, -C::elementaryCharge_C, C::electronMass_kg);
    const ParticleState before{Vector3(radius_m + 0.2, 0.0, 0.0),
                               Vector3(-testMomentum(), 0.0, 0.0), 0.0};
    const ParticleState after = integrator.step(before, 5e-8);
    require(before.r.norm() > radius_m && after.r.norm() < radius_m,
            "curved reference step does not cross the sphere");

    const SphereCrossingLocator locator(radius_m, 1e-13, 100);
    const EventLocation event = locator.locate(before, after, integrator);
    const ParticleState physical = integrator.step(before, event.state.time_s - before.time_s);
    require((event.state.r - physical.r).norm() < 1e-6,
            "event position does not match a Boris substep");
    require((event.state.p - physical.p).norm() / before.p.norm() < 1e-11,
            "event momentum was interpolated instead of reintegrated");
    require(std::abs(event.radiusResidual_m) < 1e-5,
            "curved event is not close to the sphere");
}

} // namespace

int main() {
    try {
        checkRadialCrossings();
        checkPhysicalSubstep();
        requireThrows<std::invalid_argument>([] {
            (void)SphereCrossingLocator(-1.0, 1e-12, 100);
        }, "negative radius accepted");
        requireThrows<std::invalid_argument>([] {
            (void)SphereCrossingLocator(1.0, 0.0, 100);
        }, "zero time tolerance accepted");
        requireThrows<std::invalid_argument>([] {
            (void)SphereCrossingLocator(1.0, 1e-12, 0);
        }, "zero iteration limit accepted");
        std::cout << "SphereCrossingLocator checks passed\n";
    } catch (const std::exception& error) {
        std::cerr << "SphereCrossingLocator check failed: " << error.what() << '\n';
        return 1;
    }
}
