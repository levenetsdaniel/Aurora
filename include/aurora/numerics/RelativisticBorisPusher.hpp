#pragma once

#include "aurora/fields/MagneticField.hpp"
#include "aurora/simulation/ParticleState.hpp"
#include <cstddef>

namespace aurora {
class RelativisticBorisPusher {
  public:
    class State {
        friend class RelativisticBorisPusher;
        Vec3 position_m_, momentum_before_;
        double initial_time_s_{}, dt_s_{};
        std::size_t steps_{};

      public:
        double time_s() const {
            return initial_time_s_ + static_cast<double>(steps_) * dt_s_;
        }

        Vec3 position_m() const {
            return position_m_;
        }
    };

    explicit RelativisticBorisPusher(double dt_s);
    State initialize(const ParticleState& initial, const MagneticField& field) const;
    void advance(State& state, const MagneticField& field) const;
    ParticleState synchronize(const State& state, const MagneticField& field) const;

    double dt_s() const {
        return dt_s_;
    }

  private:
    double dt_s_;
    void check_step(const State& state) const;
    static Vec3 rotate(Vec3 momentum, Vec3 field_T, double interval_s);
};
} // namespace aurora
