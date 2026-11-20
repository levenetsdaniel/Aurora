#include "aurora/visualization/DemoData.hpp"
#include "aurora/io/TrajectoryCsv.hpp"
#include "aurora/physics/Gyromotion.hpp"
#include "aurora/physics/ParticleProperties.hpp"
#include "aurora/physics/Relativity.hpp"
#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>
#include <numbers>
#include <sstream>

namespace aurora::visualization {
namespace {
std::vector<std::string> csv_columns(const std::string& line) {
    std::vector<std::string> columns;
    std::string value;
    bool quoted = false;
    for (std::size_t i = 0; i < line.size(); ++i) {
        const char c = line[i];
        if (c == '"') {
            if (quoted && i + 1 < line.size() && line[i + 1] == '"') {
                value += '"';
                ++i;
            } else
                quoted = !quoted;
        } else if (c == ',' && !quoted) {
            columns.push_back(value);
            value.clear();
        } else if (c != '\r')
            value += c;
    }
    if (quoted)
        throw std::invalid_argument("Unterminated CSV field");
    columns.push_back(value);
    return columns;
}

struct Csv {
    std::ifstream input;
    std::map<std::string, std::size_t> header;
    std::vector<std::string> row;

    explicit Csv(const std::filesystem::path& path) : input(path) {
        if (!input)
            throw std::runtime_error("Cannot open " + path.string());
        std::string line;
        if (!std::getline(input, line))
            throw std::invalid_argument("Missing CSV header");
        const auto names = csv_columns(line);
        for (std::size_t i = 0; i < names.size(); ++i)
            if (!header.emplace(names[i], i).second)
                throw std::invalid_argument("Duplicate CSV column");
    }

    bool next() {
        std::string line;
        if (!std::getline(input, line)) {
            if (input.bad())
                throw std::runtime_error("CSV read error");
            return false;
        }
        row = csv_columns(line);
        if (row.size() != header.size())
            throw std::invalid_argument("CSV column count mismatch");
        return true;
    }

    const std::string& text(const std::string& key) const {
        return row.at(header.at(key));
    }

    double number(const std::string& key) const {
        const auto& value = text(key);
        std::size_t used = 0;
        const double result = std::stod(value, &used);
        if (used != value.size() || !std::isfinite(result))
            throw std::invalid_argument("Invalid CSV value: " + key);
        return result;
    }

    ParticleState state() const {
        return {{number("x_m"), number("y_m"), number("z_m")},
            {number("px_kg_m_s"), number("py_kg_m_s"), number("pz_kg_m_s")}, number("time_s")};
    }
};

VisualSample visual_sample(const ParticleState& s, const MagneticField& field) {
    return {s.time_s, s.position_m, (norm(s.position_m) - constants::earth_radius_m) / 1000,
        particle_properties(s, field).parallel_fraction};
}

std::vector<VisualSample> read_trajectory(
    const std::filesystem::path& path, const MagneticField& field) {
    Csv csv(path);
    std::vector<VisualSample> result;
    while (csv.next()) {
        const auto state = csv.state();
        if (!result.empty() && state.time_s <= result.back().time_s)
            throw std::invalid_argument("Trajectory times must be strictly increasing");
        result.push_back(visual_sample(state, field));
    }
    if (result.empty())
        throw std::invalid_argument("Empty trajectory");
    return result;
}

void decimate(std::vector<VisualSample>& samples) {
    constexpr std::size_t maximum = 8000;
    if (samples.size() <= maximum)
        return;
    const auto stride = (samples.size() + maximum - 1) / maximum;
    std::vector<VisualSample> reduced;
    for (std::size_t i = 0; i < samples.size(); i += stride)
        reduced.push_back(samples[i]);
    if (reduced.back().time_s != samples.back().time_s)
        reduced.push_back(samples.back());
    samples = std::move(reduced);
}
}

DemoCase load_demo_case(const std::filesystem::path& directory) {
    std::ifstream manifest(directory / "manifest.txt");
    if (!manifest)
        throw std::runtime_error("Missing run manifest");
    const std::string metadata((std::istreambuf_iterator<char>(manifest)), {});
    if (metadata.find("execution_status=finished\n") == std::string::npos)
        throw std::invalid_argument("Viewer requires a finished run");
    std::ifstream config(directory / "parameters.txt");
    if (!config)
        throw std::runtime_error("Missing run parameters");
    const auto parameters = read_experiment_parameters(config);
    const auto prepared = prepare_experiment(parameters);
    DemoCase result;
    result.parameters = parameters;
    result.directory = directory;
    result.trajectory = read_trajectory(directory / "trajectory.csv", prepared.field);
    decimate(result.trajectory);
    if (std::filesystem::exists(directory / "gyro_segment.csv"))
        result.gyro_segment = read_trajectory(directory / "gyro_segment.csv", prepared.field);

    Csv summary(directory / "summary.csv");
    if (!summary.next())
        throw std::invalid_argument("Missing run summary");
    bool recognized = false;
    for (auto outcome : {ExperimentOutcome::AtmosphereReached, ExperimentOutcome::Reflected,
             ExperimentOutcome::Escaped, ExperimentOutcome::TimeLimit,
             ExperimentOutcome::StepLimit}) {
        if (summary.text("outcome") == outcome_name(outcome)) {
            result.outcome = outcome;
            recognized = true;
        }
    }
    if (!recognized)
        throw std::invalid_argument("Cannot display an invalid or failed run");
    result.energy_error = summary.number("max_relative_energy_error");
    result.mu_error = summary.number("max_relative_mu_error");
    if (summary.next())
        throw std::invalid_argument("Unexpected additional summary rows");

    Csv events(directory / "events.csv");
    const double first_time = result.trajectory.front().time_s;
    const double last_time = result.trajectory.back().time_s;
    while (events.next()) {
        const auto& name = events.text("type");
        EventType type;
        if (name == "mirror")
            type = EventType::Mirror;
        else if (name == "atmosphere_entry" || name == "inner_boundary")
            type = EventType::InnerBoundary;
        else if (name == "outer_boundary")
            type = EventType::OuterBoundary;
        else
            throw std::invalid_argument("Unknown event type");
        const auto state = events.state();
        if (state.time_s < first_time || state.time_s > last_time)
            throw std::invalid_argument("Event is outside the saved trajectory");
        result.events.push_back({type, state, events.number("bracket_width_s")});
        result.trajectory.push_back(visual_sample(state, prepared.field));
    }
    std::stable_sort(
        result.trajectory.begin(), result.trajectory.end(), [](const auto& a, const auto& b) {
            return a.time_s < b.time_s;
        });
    result.trajectory.erase(std::unique(result.trajectory.begin(), result.trajectory.end(),
                                [](const auto& a, const auto& b) {
                                    return a.time_s == b.time_s;
                                }),
        result.trajectory.end());
    return result;
}

std::vector<DemoCase> create_demo_cases(
    const ExperimentParameters& base, const std::filesystem::path& root) {
    if (!(base.energy_keV > 0))
        throw std::invalid_argument("The visual demonstration requires positive electron energy");
    const auto initial = prepare_experiment(base);
    const double speed = norm(velocity_from_momentum(initial.initial.momentum_kg_m_s));
    const double reference_speed =
        norm(velocity_from_momentum(momentum_from_energy(constants::eV_to_J(100000), {1, 0, 0})));
    const double duration = std::max(base.max_time_s, 0.12 * reference_speed / speed);
    const double needed_steps = std::ceil(duration / base.dt_s) + 1;
    if (!std::isfinite(needed_steps) || needed_steps > 100000000)
        throw std::invalid_argument(
            "Visual presets exceed 100 million steps; choose a suitable energy and step");
    std::vector<DemoCase> cases;
    for (std::size_t i = 0; i < 3; ++i) {
        auto p = base;
        p.pitch_deg = i == 0 ? base.pitch_deg : (i == 1 ? 5.0 : 175.0);
        p.stop_at_first_mirror = i != 0;
        p.max_time_s = duration;
        p.max_steps = std::max(base.max_steps, static_cast<std::size_t>(needed_steps));
        const auto prepared = prepare_experiment(p);
        const double segment_end = 3 *
            gyroperiod_s(prepared.initial.momentum_kg_m_s,
                prepared.field.at(prepared.initial.position_m, 0));
        std::vector<std::pair<std::size_t, ParticleState>> dense;
        const auto saved =
            save_experiment(p, root, [&](std::size_t step, const ParticleState& state) {
                if (state.time_s <= segment_end && dense.size() < 20000)
                    dense.emplace_back(step, state);
            });
        if (saved.result.outcome == ExperimentOutcome::NumericalError)
            throw std::runtime_error(saved.result.trajectory.error_message);
        std::ofstream segment(saved.directory / "gyro_segment.csv");
        segment.exceptions(std::ios::failbit | std::ios::badbit);
        write_magnetic_trajectory_header(segment);
        for (const auto& [step, state] : dense)
            write_magnetic_trajectory_row(segment, step, state, prepared.field);
        segment.close();
        cases.push_back(load_demo_case(saved.directory));
        std::cout << "VTK example " << i + 1 << ": " << outcome_name(saved.result.outcome)
                  << ", events=" << saved.result.trajectory.events.size() << ", "
                  << saved.directory.string() << '\n'
                  << std::flush;
    }
    return cases;
}

VisualSample sample_at_time(const std::vector<VisualSample>& samples, double time) {
    if (samples.empty() || !std::isfinite(time))
        throw std::invalid_argument("Invalid animation time or samples");
    if (time <= samples.front().time_s)
        return samples.front();
    if (time >= samples.back().time_s)
        return samples.back();
    const auto end =
        std::upper_bound(samples.begin(), samples.end(), time, [](double t, const VisualSample& s) {
            return t < s.time_s;
        });
    const auto& b = *end;
    const auto& a = *(end - 1);
    const double f = (time - a.time_s) / (b.time_s - a.time_s);
    return {time, a.position_m + (b.position_m - a.position_m) * f,
        a.altitude_km + (b.altitude_km - a.altitude_km) * f,
        a.parallel_fraction + (b.parallel_fraction - a.parallel_fraction) * f};
}

std::vector<Vec3> dipole_field_line(
    const DipoleField& field, Vec3 start, double atmosphere_radius, std::size_t count) {
    if (count < 2 || count > 100000)
        throw std::invalid_argument("Invalid field line point count");
    (void)equatorial_loss_cone(field, start, atmosphere_radius);
    const double r_eq = norm(start);
    const Vec3 radial = normalized(start);
    const double latitude_limit = std::acos(std::sqrt(atmosphere_radius / r_eq));
    std::vector<Vec3> points;
    for (std::size_t i = 0; i < count; ++i) {
        const double latitude = -latitude_limit +
            2 * latitude_limit * static_cast<double>(i) / static_cast<double>(count - 1);
        const double cosine = std::cos(latitude);
        // Dipole field line: r(lambda) = r_eq * cos^2(lambda).
        const double r = r_eq * cosine * cosine;
        points.push_back((radial * cosine + field.axis() * std::sin(latitude)) * r);
    }
    return points;
}
} // namespace aurora::visualization
