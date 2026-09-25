#ifndef AURORA_PARTICLE_HPP
#define AURORA_PARTICLE_HPP


#include "../math/Vector3.hpp"


struct ParticleState {
    Vector3 r;
    Vector3 p;
    double time_s{};
};

enum class ParticleOutcome { Precipitated, Mirrored, Escaped, Unresolved };


#endif //AURORA_PARTICLE_HPP