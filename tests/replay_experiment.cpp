#include "aurora/io/ExperimentFiles.hpp"
#include <fstream>
#include <iostream>

int main(int argc, char** argv) {
    try {
        if (argc < 2 || argc > 3)
            throw std::invalid_argument("Usage: aurora_replay parameters.txt [output-root]");
        std::ifstream input(argv[1]);
        if (!input)
            throw std::runtime_error("Cannot open parameters file");
        const auto p = aurora::read_experiment_parameters(input);
        const auto saved = aurora::save_experiment(p, argc == 3 ? argv[2] : "results/replay");
        std::cout << aurora::outcome_name(saved.result.outcome) << ": " << saved.directory.string()
                  << '\n';
        return saved.result.outcome == aurora::ExperimentOutcome::NumericalError ? 1 : 0;
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
