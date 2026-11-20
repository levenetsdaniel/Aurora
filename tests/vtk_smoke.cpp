#include "aurora/visualization/VtkViewer.hpp"
#include <fstream>
#include <iostream>

int main(int argc, char** argv) {
    try {
        using namespace aurora::visualization;
        const std::filesystem::path output = argc > 1 ? argv[1] : "results/visualization-preview";
        std::filesystem::create_directories(output);
        std::vector<DemoCase> cases;
        std::ifstream previous(output / "runs.txt");
        std::string path;
        while (std::getline(previous, path)) {
            cases.push_back(load_demo_case(path));
        }

        if (cases.empty()) {
            cases = create_demo_cases(aurora::ExperimentParameters{}, output);
            std::ofstream runs(output / "runs.txt");
            runs.exceptions(std::ios::failbit | std::ios::badbit);
            for (const auto& demo : cases) {
                runs << std::filesystem::absolute(demo.directory).string() << '\n';
            }
        }

        if (cases.size() != 3 || cases[0].events.empty() || cases[0].gyro_segment.size() < 100 ||
            cases[1].outcome != aurora::ExperimentOutcome::AtmosphereReached ||
            cases[2].outcome != aurora::ExperimentOutcome::AtmosphereReached ||
            cases[1].events.front().state.position_m.z *
                    cases[2].events.front().state.position_m.z >=
                0) {
            throw std::runtime_error(
                "Visual demonstration does not contain the expected physical cases");
        }

        for (std::size_t i = 0; i < cases.size(); ++i) {
            render_demo_snapshot(cases, i, 1, output / ("case-" + std::to_string(i + 1) + ".png"));
        }
        render_demo_snapshot(cases, 0, 0.35, output / "case-1-moving.png");
        std::cout << "Saved four VTK frames in " << output << '\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "FAIL: " << error.what() << '\n';
        return 1;
    }
}
