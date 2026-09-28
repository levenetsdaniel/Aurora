#ifndef AURORA_EXCITEDSTATEKINETICS_HPP
#define AURORA_EXCITEDSTATEKINETICS_HPP

#include "ExcitationModel.hpp"
#include "../atmosphere/AtmosphereModel.hpp"
#include "../physics/Constants.hpp"

#include <cmath>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

struct RadiativeState {
    std::string id;
    std::string species;
    double totalDecay_s1{};
    std::string source;
};

struct QuenchingTerm {
    std::string stateId;
    std::string quencher;
    double coefficient_m3_s{};
    double referenceTemperature_K{};
    double temperatureExponent{};
    double activationTemperature_K{};
    std::string source;

    // k(T) = coefficient * (T/Tref)^exponent * exp(-activation/T).
    double rateCoefficient_m3_s(double temperature_K) const {
        if (!std::isfinite(temperature_K) || temperature_K <= 0.0) {
            throw std::invalid_argument("QuenchingTerm: invalid temperature");
        }
        const double value = coefficient_m3_s *
            std::pow(temperature_K / referenceTemperature_K, temperatureExponent) *
            std::exp(-activationTemperature_K / temperature_K);
        if (!std::isfinite(value) || value < 0.0) {
            throw std::overflow_error("QuenchingTerm: rate coefficient overflow");
        }
        return value;
    }
};

struct StateAlias {
    ExcitedStateId eventState;
    std::string stateId;
};

class KineticsCatalog {
private:
    std::map<std::string, RadiativeState> statesById;
    std::multimap<std::string, QuenchingTerm> quenchingByState;
    std::map<ExcitedStateId, std::string> aliases;

public:
    KineticsCatalog(std::vector<RadiativeState> states,
                    std::vector<QuenchingTerm> terms,
                    std::vector<StateAlias> aliases) {
        if (states.empty()) {
            throw std::invalid_argument("KineticsCatalog: no states");
        }
        for (RadiativeState& state : states) {
            const std::string id = state.id;
            if (state.id.empty() || state.species.empty() || state.source.empty() ||
                !std::isfinite(state.totalDecay_s1) || state.totalDecay_s1 <= 0.0 ||
                !statesById.emplace(id, std::move(state)).second) {
                throw std::invalid_argument("KineticsCatalog: invalid or duplicate state");
            }
        }
        for (QuenchingTerm& term : terms) {
            const std::string id = term.stateId;
            if (!statesById.contains(term.stateId) ||
                (term.quencher != "O" && term.quencher != "O2" &&
                 term.quencher != "N2") ||
                !std::isfinite(term.coefficient_m3_s) || term.coefficient_m3_s < 0.0 ||
                !std::isfinite(term.referenceTemperature_K) ||
                term.referenceTemperature_K <= 0.0 ||
                !std::isfinite(term.temperatureExponent) ||
                !std::isfinite(term.activationTemperature_K) ||
                term.source.empty()) {
                throw std::invalid_argument("KineticsCatalog: invalid quenching term");
            }
            quenchingByState.emplace(id, std::move(term));
        }
        for (StateAlias& alias : aliases) {
            const auto state = statesById.find(alias.stateId);
            if (alias.eventState.species.empty() || alias.eventState.state.empty() ||
                state == statesById.end() ||
                state->second.species != alias.eventState.species ||
                !this->aliases.emplace(alias.eventState, alias.stateId).second) {
                throw std::invalid_argument("KineticsCatalog: invalid or duplicate alias");
            }
        }
        for (const auto& [id, state] : statesById) {
            if (quenchingByState.count(id) == 0) {
                throw std::invalid_argument("KineticsCatalog: state has no quenching data");
            }
            for (const std::string& quencher : {"O", "O2", "N2"}) {
                bool found = false;
                const auto [first, last] = quenchingByState.equal_range(id);
                for (auto term = first; term != last; ++term) {
                    found |= term->second.quencher == quencher;
                }
                if (!found) {
                    throw std::invalid_argument("KineticsCatalog: missing quencher " +
                                                quencher + " for " + id);
                }
            }
        }
    }

    const RadiativeState& resolve(const ExcitedStateId& eventState) const {
        const auto alias = aliases.find(eventState);
        if (alias == aliases.end()) {
            throw std::out_of_range("KineticsCatalog: no coefficients for " +
                                    eventState.species + ":" + eventState.state);
        }
        return statesById.at(alias->second);
    }

    const RadiativeState& stateById(const std::string& id) const {
        return statesById.at(id);
    }

    double radiativeLifetime_s(const ExcitedStateId& eventState) const {
        return 1.0 / resolve(eventState).totalDecay_s1;
    }

    std::vector<QuenchingTerm> termsFor(const ExcitedStateId& eventState) const {
        const auto [first, last] = quenchingByState.equal_range(resolve(eventState).id);
        std::vector<QuenchingTerm> result;
        for (auto term = first; term != last; ++term) result.push_back(term->second);
        return result;
    }

    double quenchingFrequency_s1(const std::string& stateId,
                                  const AtmosphereSample& atmosphere) const {
        if (!statesById.contains(stateId)) {
            throw std::out_of_range("KineticsCatalog: unknown state");
        }
        if (!std::isfinite(atmosphere.temperature_K) ||
            atmosphere.temperature_K <= 0.0 ||
            !std::isfinite(atmosphere.nO_m3) || atmosphere.nO_m3 < 0.0 ||
            !std::isfinite(atmosphere.nO2_m3) || atmosphere.nO2_m3 < 0.0 ||
            !std::isfinite(atmosphere.nN2_m3) || atmosphere.nN2_m3 < 0.0) {
            throw std::invalid_argument("KineticsCatalog: invalid atmosphere sample");
        }
        double total_s1 = 0.0;
        const auto [first, last] = quenchingByState.equal_range(stateId);
        for (auto term = first; term != last; ++term) {
            const QuenchingTerm& coefficient = term->second;
            const double density_m3 = coefficient.quencher == "O" ? atmosphere.nO_m3 :
                coefficient.quencher == "O2" ? atmosphere.nO2_m3 : atmosphere.nN2_m3;
            total_s1 += coefficient.rateCoefficient_m3_s(atmosphere.temperature_K) *
                        density_m3;
        }
        if (!std::isfinite(total_s1) || total_s1 < 0.0) {
            throw std::overflow_error("KineticsCatalog: quenching frequency overflow");
        }
        return total_s1;
    }
};

struct StatePopulation {
    double production_m3_s1{};
    double radiativeDecay_s1{};
    double quenching_s1{};
    double effectiveLifetime_s{};
    double density_m3{};
};

class ExcitedStateKinetics {
private:
    const KineticsCatalog& catalog;
    const AtmosphereModel& atmosphere;

public:
    ExcitedStateKinetics(const KineticsCatalog& catalog,
                         const AtmosphereModel& atmosphere)
        : catalog(catalog), atmosphere(atmosphere) {}

    StatePopulation solve(const ExcitationModel& excitation,
                          const ExcitedStateId& state,
                          std::size_t x, std::size_t y, std::size_t z) const {
        const RadiativeState& coefficients = catalog.resolve(state);
        const Vector3 center_m = excitation.cellCenter_m(x, y, z);
        const double altitude_m = center_m.norm() - C::earthRadius_m;
        const AtmosphereSample sample = atmosphere.sample(altitude_m, 0.0, 0.0, 0.0);
        const double production_m3_s1 = excitation.productionRate_m3_s1(state, x, y, z);
        const double quenching_s1 = catalog.quenchingFrequency_s1(coefficients.id, sample);
        const double loss_s1 = coefficients.totalDecay_s1 + quenching_s1;
        if (!std::isfinite(loss_s1)) {
            throw std::overflow_error("ExcitedStateKinetics: total loss rate overflow");
        }
        const double lifetime_s = 1.0 / loss_s1;
        const double density_m3 = production_m3_s1 * lifetime_s;
        if (!std::isfinite(lifetime_s) || !std::isfinite(density_m3) ||
            density_m3 < 0.0) {
            throw std::overflow_error("ExcitedStateKinetics: population overflow");
        }
        return {production_m3_s1, coefficients.totalDecay_s1,
                quenching_s1, lifetime_s, density_m3};
    }

    std::vector<StatePopulation> solveGrid(const ExcitationModel& excitation,
                                            const ExcitedStateId& state) const {
        const auto& grid = excitation.configuration();
        std::vector<StatePopulation> cells;
        cells.reserve(grid.cellsX * grid.cellsY * grid.cellsZ);
        for (std::size_t z = 0; z < grid.cellsZ; ++z) {
            for (std::size_t y = 0; y < grid.cellsY; ++y) {
                for (std::size_t x = 0; x < grid.cellsX; ++x) {
                    cells.push_back(solve(excitation, state, x, y, z));
                }
            }
        }
        return cells;
    }
};

#endif // AURORA_EXCITEDSTATEKINETICS_HPP
