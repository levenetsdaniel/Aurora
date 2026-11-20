#pragma once

#include "aurora/visualization/DemoData.hpp"

namespace aurora::visualization {
void show_demo(const std::vector<DemoCase>& cases);
// A deterministic noninteractive render for screenshot review and graphics smoke checks.
void render_demo_snapshot(const std::vector<DemoCase>& cases, std::size_t selected_case,
    double time_fraction, const std::filesystem::path& png_path);
} // namespace aurora::visualization
