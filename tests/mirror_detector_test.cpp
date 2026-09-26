#include "../src/simulation/MirrorDetector.hpp"
#include "../src/physics/Constants.hpp"
#include "../src/physics/Relativity.hpp"
#include "../src/physics/UniformField.hpp"

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

class TiltingField final : public MagneticField {
private:
    double centerX_m;
    double scale_m;
    double strength_T;

public:
    TiltingField(double centerX_m, double scale_m, double strength_T)
        : centerX_m(centerX_m), scale_m(scale_m), strength_T(strength_T) {}

    Vector3 evaluate(const Vector3& r) const override {
        return Vector3((r.x - centerX_m) / scale_m, 0.0, 1.0) * strength_T;
    }
};

double electronMomentum() {
    return Relativity::momentumMagnitudeFromKineticEnergy_J(
        1000.0 * C::elementaryCharge_C, C::electronMass_kg);
}

void checkMirrorLocation() {
    const double centerX_m = C::earthRadius_m;
    const double height_m = 100000.0;
    const TiltingField field(centerX_m, 100.0, 1e-9);
    const BorisIntegrator integrator(field, 0.0, C::electronMass_kg);
    const double p = electronMomentum();
    const double speed_m_s = Relativity::velocityFromMomentum(
        Vector3(p, 0.0, 0.0), C::electronMass_kg).norm();
    const ParticleState before{Vector3(centerX_m - 100.0, 0.0, height_m),
                               Vector3(p, 0.0, 0.0), 0.0};
    const ParticleState after = integrator.step(before, 1e-5);
    const double expectedTime_s = 100.0 / speed_m_s;

    const MirrorDetector loose(p * 0.01, 1e-7, 100);
    const MirrorDetector tight(p * 0.01, 1e-12, 100);
    require(tight.bracketed(before, after, field), "parallel momentum reversal was not bracketed");

    const MirrorEvent coarse = loose.locate(before, after, field, integrator);
    const MirrorEvent event = tight.locate(before, after, field, integrator);
    const double expectedLatitude_rad = std::asin(height_m /
        std::sqrt(centerX_m * centerX_m + height_m * height_m));
    require(std::abs(event.state.time_s - expectedTime_s) < 1e-12,
            "mirror time differs from known straight-line root");
    require(std::abs(event.latitude_rad - expectedLatitude_rad) < 1e-12,
            "mirror latitude does not match event position");
    require(std::abs(event.state.p.dot(field.evaluate(event.state.r)) /
                     field.evaluate(event.state.r).norm()) < p * 1e-6,
            "parallel momentum at located event is not close to zero");
    require(std::abs(event.state.time_s - expectedTime_s) <=
            std::abs(coarse.state.time_s - expectedTime_s) + 1e-12,
            "tightening time tolerance worsened mirror time");

    const ParticleState physical = integrator.step(before, event.state.time_s - before.time_s);
    require((event.state.r - physical.r).norm() < 1e-6 &&
            (event.state.p - physical.p).norm() / p < 1e-12,
            "located mirror state does not match a Boris substep");

    requireThrows<std::runtime_error>([&] {
        (void)MirrorDetector(p * 0.01, 1e-15, 1).locate(before, after, field, integrator);
    }, "iteration limit did not report failure");
}

void checkFalseCandidates() {
    const double p = electronMomentum();
    const UniformField field(Vector3(0.0, 0.0, 1e-9));
    const MirrorDetector detector(p * 0.01);
    const ParticleState before{Vector3(C::earthRadius_m, 0.0, -100.0),
                               Vector3(0.0, 0.0, p), 0.0};
    const ParticleState after{Vector3(C::earthRadius_m, 0.0, 100.0),
                              Vector3(0.0, 0.0, p), 1e-5};
    require(!detector.bracketed(before, after, field),
            "changing position was mistaken for momentum reversal");

    const ParticleState smallPositive{before.r, Vector3(0.0, 0.0, p * 0.001), 0.0};
    const ParticleState smallNegative{after.r, Vector3(0.0, 0.0, -p * 0.001), 1e-5};
    require(!detector.bracketed(smallPositive, smallNegative, field),
            "momentum tolerance did not suppress a tiny sign change");

    const BorisIntegrator integrator(field, 0.0, C::electronMass_kg);
    requireThrows<std::invalid_argument>([&] {
        (void)detector.locate(before, after, field, integrator);
    }, "locate accepted states without a bracketed reversal");
}

} // namespace

int main() {
    try {
        checkMirrorLocation();
        checkFalseCandidates();
        requireThrows<std::invalid_argument>([] {
            (void)MirrorDetector(-1.0);
        }, "negative momentum tolerance accepted");
        requireThrows<std::invalid_argument>([] {
            (void)MirrorDetector(0.0, 0.0);
        }, "zero time tolerance accepted");
        std::cout << "MirrorDetector checks passed\n";
    } catch (const std::exception& error) {
        std::cerr << "MirrorDetector check failed: " << error.what() << '\n';
        return 1;
    }
}
