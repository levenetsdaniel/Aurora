#include "aurora/visualization/DemoData.hpp"
#include "aurora/io/EventCsv.hpp"
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>

namespace {
using namespace aurora;
using namespace aurora::visualization;

void require(bool condition, const char* message) {
    if (!condition) {
        throw std::runtime_error(message);
    }
}

template <class Action> void require_rejected(Action action) {
    bool rejected = false;
    try {
        action();
    } catch (const std::exception&) {
        rejected = true;
    }
    require(rejected, "Invalid visualization data was accepted");
}

std::string contents(const std::filesystem::path& path) {
    std::ifstream input(path);
    require(static_cast<bool>(input), "Missing experiment file");
    return {std::istreambuf_iterator<char>(input), {}};
}

void interpolation() {
    const std::vector<VisualSample> points{
        {0, {0, 0, 0}, 100, 1}, {1, {2, 4, 6}, 200, 0}, {4, {8, 16, 24}, 500, -1}};
    const auto sample = sample_at_time(points, 2.5);
    require(norm(sample.position_m - Vec3{5, 10, 15}) < 1e-12 && sample.altitude_km == 350 &&
            sample.parallel_fraction == -0.5,
        "Interpolation must use physical time on a nonuniform grid");
    require(sample_at_time(points, -1).time_s == 0 && sample_at_time(points, 10).time_s == 4,
        "Animation must clamp to the saved time interval");
    require_rejected([&] {
        sample_at_time({}, 0);
    });
    require_rejected([&] {
        sample_at_time(points, std::numeric_limits<double>::quiet_NaN());
    });
}

void field_geometry() {
    const DipoleField field;
    const double radius = constants::earth_radius_m;
    const auto points = dipole_field_line(field, {2 * radius, 0, 0}, radius + 100000);
    require(std::abs(norm(points.front()) - radius - 100000) < 1e-6 &&
            std::abs(norm(points.back()) - radius - 100000) < 1e-6,
        "Field line must end on the atmospheric sphere");
    require(norm(points[points.size() / 2] - Vec3{2 * radius, 0, 0}) < 1e-6,
        "Field line must pass through the equatorial starting point");
    for (const auto& point : points) {
        const double sine = dot(normalized(point), field.axis());
        require(std::abs(norm(point) / (1 - sine * sine) - 2 * radius) < 1e-6,
            "Field line violates the dipole shell invariant");
    }
    require_rejected([&] {
        dipole_field_line(field, {2 * radius, 0, 0}, radius, 1);
    });
    require_rejected([&] {
        dipole_field_line(field, {0, 0, 2 * radius}, radius);
    });
}

void saved_data() {
    ExperimentParameters parameters;
    parameters.max_steps = 10;
    parameters.sample_every = 3;
    const auto ordinary = save_experiment(parameters, "test-results/visualization");
    std::size_t observed = 0;
    const auto dense = save_experiment(
        parameters, "test-results/visualization", [&](std::size_t step, const ParticleState&) {
            require(step == observed, "Observer missed or repeated an integration step");
            ++observed;
        });
    require(observed == 11, "Dense observation must include initial and final states");
    require(contents(ordinary.directory / "trajectory.csv") ==
            contents(dense.directory / "trajectory.csv"),
        "Dense observation changed the sampled trajectory CSV");
    require(
        contents(ordinary.directory / "summary.csv") == contents(dense.directory / "summary.csv"),
        "Dense observation changed the physical result");
    const auto loaded = load_demo_case(dense.directory);
    require(loaded.trajectory.size() == 5 && loaded.gyro_segment.empty(),
        "Loader must preserve the terminal sample and allow runs without a dense segment");
    require(loaded.outcome == ExperimentOutcome::StepLimit && loaded.events.empty(),
        "Loader must preserve an unresolved physical outcome");

    const auto prepared = prepare_experiment(parameters);
    std::ofstream events(dense.directory / "events.csv");
    write_event_header(events);
    for (double time : {2 * parameters.dt_s, 7 * parameters.dt_s}) {
        auto state = prepared.initial;
        state.time_s = time;
        write_event_row(events, {EventType::Mirror, state, 1e-12}, prepared.field);
    }
    events.close();
    const auto with_events = load_demo_case(dense.directory);
    require(with_events.events.size() == 2 && with_events.trajectory.size() == 7,
        "Loader must insert multiple exact event samples into the saved time interval");

    std::ofstream manifest(dense.directory / "manifest.txt");
    manifest << "execution_status=in_progress\n";
    manifest.close();
    require_rejected([&] {
        load_demo_case(dense.directory);
    });
    std::ofstream trajectory(ordinary.directory / "trajectory.csv", std::ios::app);
    trajectory << "broken,row\n";
    trajectory.close();
    require_rejected([&] {
        load_demo_case(ordinary.directory);
    });
}
}

int main() {
    try {
        interpolation();
        field_geometry();
        saved_data();
        std::cout << "Visualization data checks passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "FAIL: " << error.what() << '\n';
        return 1;
    }
}
