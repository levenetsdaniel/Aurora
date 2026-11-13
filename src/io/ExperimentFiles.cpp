#include "aurora/io/ExperimentFiles.hpp"
#include "aurora/io/TrajectoryCsv.hpp"
#include "aurora/io/EventCsv.hpp"
#include "aurora/physics/Relativity.hpp"
#include "aurora/BuildInfo.hpp"
#include <charconv>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <map>
#include <string>

namespace aurora {
namespace {
long long now_ms() {
    return std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::system_clock::now().time_since_epoch())
        .count();
}

std::ofstream file(const std::filesystem::path& path) {
    std::ofstream out(path);
    out.exceptions(std::ios::failbit | std::ios::badbit);
    out << std::setprecision(17);
    return out;
}

std::filesystem::path new_directory(const std::filesystem::path& root) {
    std::filesystem::create_directories(root);

    const std::string prefix = "run-" + std::to_string(now_ms());

    for (int i = 0; i < 10000; ++i) {
        const auto path = std::filesystem::absolute(root / (prefix + "-" + std::to_string(i)));
        if (std::filesystem::create_directory(path))
            return path;
    }

    throw std::runtime_error("Cannot allocate a unique run directory");
}

void manifest(const std::filesystem::path& directory, const char* status, long long started) {
    auto out = file(directory / "manifest.txt");
    out << "schema_version=1\nexecution_status=" << status << "\nstarted_unix_ms=" << started
        << "\nupdated_unix_ms=" << now_ms() << "\nsource_sha256=" << build_info::source_sha256
        << "\ngit_commit_at_configuration=" << build_info::git_commit
        << "\ngit_state_at_configuration=" << build_info::git_state
        << "\ncompiler=" << build_info::compiler << "\nplatform=" << build_info::platform
        << "\nbuild_type=" << AURORA_BUILD_TYPE << "\ncompiler_flags=" << build_info::compiler_flags
        << "\ncxx_standard=20\nmodel=static_centered_dipole_E_zero\ncoordinates=geocentric_cartesian_SI"
        << "\nrandom_generator=none\nsource_manifest=source/source.sha256"
        << "\nelectron_mass_kg=" << constants::electron_mass_kg
        << "\nelectron_charge_C=" << constants::electron_charge_C
        << "\nspeed_of_light_m_s=" << constants::speed_of_light_m_s
        << "\nearth_radius_m=" << constants::earth_radius_m << '\n';
    out.close();
}

std::string quote_csv(const std::string& value) {
    std::string result = "\"";
    for (char c : value) {
        if (c == '"')
            result += '"';
        result += c;
    }
    return result + '"';
}

void summary(const std::filesystem::path& directory, const PreparedExperiment& e,
    const ExperimentResult& r) {
    const auto& s = r.trajectory.final_state;
    auto out = file(directory / "summary.csv");
    out << "outcome,prediction,agreement,stop_reason,resolved,completed_steps,events,final_time_s,"
           "final_energy_J,max_relative_energy_error,max_relative_mu_error,x_m,y_m,z_m,"
           "px_kg_m_s,py_kg_m_s,pz_kg_m_s,error_message\n"
        << outcome_name(r.outcome) << ','
        << (e.prediction ? prediction_name(*e.prediction) : "not_applicable") << ','
        << agreement_name(r.agreement) << ',' << stop_reason_name(r.trajectory.stop_reason) << ','
        << is_resolved(r.outcome) << ',' << r.trajectory.completed_steps << ','
        << r.trajectory.events.size() << ',' << s.time_s << ','
        << kinetic_energy_J(s.momentum_kg_m_s) << ',' << r.trajectory.max_relative_energy_error
        << ',' << r.trajectory.max_relative_magnetic_moment_error << ',' << s.position_m.x << ','
        << s.position_m.y << ',' << s.position_m.z << ',' << s.momentum_kg_m_s.x << ','
        << s.momentum_kg_m_s.y << ',' << s.momentum_kg_m_s.z << ','
        << quote_csv(r.trajectory.error_message) << '\n';
    out.close();
}

std::string trim(const std::string& text) {
    const auto begin = text.find_first_not_of(" \t\r");
    return begin == std::string::npos
        ? ""
        : text.substr(begin, text.find_last_not_of(" \t\r") - begin + 1);
}
}

void write_experiment_parameters(std::ostream& out, const ExperimentParameters& p) {
    out << std::setprecision(17) << "schema_version=1\nenergy_keV=" << p.energy_keV
        << "\npitch_deg=" << p.pitch_deg << "\ngyro_phase_deg=" << p.gyro_phase_deg
        << "\nequatorial_distance_RE=" << p.equatorial_distance_RE << "\ndt_s=" << p.dt_s
        << "\nmax_time_s=" << p.max_time_s << "\natmosphere_height_m=" << p.atmosphere_height_m
        << "\nouter_distance_RE=" << p.outer_distance_RE
        << "\ndipole_equatorial_T=" << p.dipole_equatorial_T
        << "\ndipole_axis_x=" << p.dipole_axis.x << "\ndipole_axis_y=" << p.dipole_axis.y
        << "\ndipole_axis_z=" << p.dipole_axis.z << "\nmax_steps=" << p.max_steps
        << "\nsample_every=" << p.sample_every << "\nevent_tolerance_s=" << p.event_tolerance_s
        << "\nmirror_hysteresis=" << p.mirror_hysteresis
        << "\nstop_at_first_mirror=" << (p.stop_at_first_mirror ? 1 : 0) << '\n';
    if (!out)
        throw std::runtime_error("Cannot write experiment parameters");
}

ExperimentParameters read_experiment_parameters(std::istream& input) {
    std::map<std::string, std::string> values;
    std::string line;
    while (std::getline(input, line)) {
        line = trim(line);
        if (line.empty() || line.starts_with('#'))
            continue;
        const auto separator = line.find('=');
        if (separator == std::string::npos)
            throw std::invalid_argument("Expected key=value in parameters");
        const auto key = trim(line.substr(0, separator));
        const auto value = trim(line.substr(separator + 1));
        if (key.empty() || value.empty() || !values.emplace(key, value).second)
            throw std::invalid_argument("Empty or duplicate parameter: " + key);
    }
    if (input.bad())
        throw std::runtime_error("Cannot read experiment parameters");

    auto take = [&](const std::string& key) {
        const auto it = values.find(key);
        if (it == values.end())
            throw std::invalid_argument("Missing parameter: " + key);
        const auto value = it->second;
        values.erase(it);
        return value;
    };

    auto number = [&](const std::string& key) {
        const auto value = take(key);
        std::size_t used = 0;
        const double result = std::stod(value, &used);
        if (used != value.size() || !std::isfinite(result))
            throw std::invalid_argument("Invalid number: " + key);
        return result;
    };

    auto count = [&](const std::string& key) {
        const auto value = take(key);
        std::size_t result = 0;
        const auto parsed = std::from_chars(value.data(), value.data() + value.size(), result);
        if (parsed.ec != std::errc{} || parsed.ptr != value.data() + value.size())
            throw std::invalid_argument("Invalid integer: " + key);
        return result;
    };
    if (take("schema_version") != "1")
        throw std::invalid_argument("Unsupported parameter schema");

    ExperimentParameters p;
    p.energy_keV = number("energy_keV");
    p.pitch_deg = number("pitch_deg");
    p.gyro_phase_deg = number("gyro_phase_deg");
    p.equatorial_distance_RE = number("equatorial_distance_RE");
    p.dt_s = number("dt_s");
    p.max_time_s = number("max_time_s");
    p.atmosphere_height_m = number("atmosphere_height_m");
    p.outer_distance_RE = number("outer_distance_RE");
    p.dipole_equatorial_T = number("dipole_equatorial_T");
    p.dipole_axis = {number("dipole_axis_x"), number("dipole_axis_y"), number("dipole_axis_z")};
    p.max_steps = count("max_steps");
    p.sample_every = count("sample_every");
    p.event_tolerance_s = number("event_tolerance_s");
    p.mirror_hysteresis = number("mirror_hysteresis");
    const auto stop = count("stop_at_first_mirror");
    if (stop > 1)
        throw std::invalid_argument("stop_at_first_mirror must be 0 or 1");
    p.stop_at_first_mirror = stop != 0;

    if (!values.empty())
        throw std::invalid_argument("Unknown parameter: " + values.begin()->first);

    (void)prepare_experiment(p);
    return p;
}

SavedExperiment save_experiment(
    const ExperimentParameters& parameters, const std::filesystem::path& root) {
    const auto e = prepare_experiment(parameters);
    const auto directory = new_directory(root);
    const auto started_ms = now_ms();
    manifest(directory, "in_progress", started_ms);
    try {
        auto config = file(directory / "parameters.txt");
        write_experiment_parameters(config, parameters);
        config.close();
        std::filesystem::copy(build_info::snapshot_directory, directory / "source",
            std::filesystem::copy_options::recursive);
        auto csv = file(directory / "trajectory.csv");
        write_magnetic_trajectory_header(csv);
        const auto start = std::chrono::steady_clock::now();

        auto result = run_atmosphere_experiment(e, [&](std::size_t step, const ParticleState& s) {
            write_magnetic_trajectory_row(csv, step, s, e.field);
        });
        const double elapsed =
            std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();

        csv.close();

        auto events = file(directory / "events.csv");
        write_event_header(events);
        for (const auto& event : result.trajectory.events)
            write_event_row(events, event, e.field, "atmosphere_entry");
        events.close();

        summary(directory, e, result);
        manifest(directory,
            result.outcome == ExperimentOutcome::NumericalError ? "numerical_error" : "finished",
            started_ms);

        return {std::move(result), directory, elapsed};
    } catch (...) {
        try {
            manifest(directory, "failed", started_ms);
        } catch (...) {
        }
        throw;
    }
}
} // namespace aurora
