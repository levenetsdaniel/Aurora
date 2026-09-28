#include "../src/emission/OpticalRenderer.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <numbers>
#include <sstream>
#include <stdexcept>

namespace {

void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

template <typename Exception, typename Action>
void requireThrows(Action action, const char* message) {
    try { action(); } catch (const Exception&) { return; }
    throw std::runtime_error(message);
}

bool near(double actual, double expected) {
    return std::abs(actual - expected) <= 1e-12 * std::max(1.0, std::abs(expected));
}

VolumeEmissionGrid cube() {
    VolumeEmissionGrid volume({{0, 0, 0}, {1, 2, 2}, 2, 1, 1},
                              {{"red", 630, "red", 1, "test"},
                               {"green", 557.7, "green", 1, "test"}});
    volume.setEmissivity("red", 0, 0, 0, 2.0);
    volume.setEmissivity("red", 1, 0, 0, 4.0);
    volume.setEmissivity("green", 0, 0, 0, 8.0);
    volume.setEmissivity("green", 1, 0, 0, 8.0);
    return volume;
}

void checkAnalyticRays() {
    const auto volume = cube();
    const OpticalRenderer camera({{-1, 1, 1}, {1, 0, 0}, {0, 0, 1}, 60, 1, 1});
    const auto images = camera.render(volume);
    require(images.size() == 2 &&
            near(images[0].at(0, 0), 6.0 / (4.0 * std::numbers::pi)) &&
            near(images[1].at(0, 0), 16.0 / (4.0 * std::numbers::pi)),
            "uniform cube or multi-voxel path integral is wrong");
    const OpticalRenderer inside({{0.5, 1, 1}, {1, 0, 0}, {0, 0, 1}, 60, 1, 1});
    require(near(inside.render(volume)[0].at(0, 0),
                 5.0 / (4.0 * std::numbers::pi)),
            "camera inside volume has wrong path length");
    const OpticalRenderer away({{-1, 1, 1}, {-1, 0, 0}, {0, 0, 1}, 60, 1, 1});
    require(away.render(volume)[0].at(0, 0) == 0.0,
            "ray looking away intersects volume");
}

void checkPreviewAndValidation() {
    const auto volume = cube();
    const OpticalRenderer camera({{-1, 1, 1}, {2, 0, 0}, {0, 0, 3}, 60, 1, 1});
    const auto images = camera.render(volume);
    const auto preview = OpticalRenderer::makeRgbPreview(images[0], images[1], 1.0);
    require(preview.rgb.size() == 3 && preview.rgb[0] > 0 &&
            preview.rgb[1] > 0 && preview.rgb[2] == 0,
            "RGB preview did not map spectral images to display channels");
    std::ostringstream spectralCsv;
    images[0].writeCsv(spectralCsv);
    require(spectralCsv.str().find("x,y,I_ph_m2_s1_sr1\n") != std::string::npos &&
            spectralCsv.str().find("0,0,") != std::string::npos,
            "spectral CSV is missing radiance or units");
    std::ostringstream rgbPpm;
    preview.writePpm(rgbPpm);
    require(rgbPpm.str().substr(0, 11) == "P6\n1 1\n255\n" &&
            rgbPpm.str().size() == 14,
            "RGB PPM encoding is wrong");
    requireThrows<std::invalid_argument>([&] {
        auto unsafe = images;
        unsafe[0].lineId = "../outside";
        OpticalRenderer::saveSpectralCsv(unsafe, "unused-output");
    }, "unsafe line ID was accepted as a filename");
    requireThrows<std::invalid_argument>([] {
        (void)OpticalRenderer({{0, 0, 0}, {1, 0, 0}, {1, 0, 0}, 60, 1, 1});
    }, "parallel camera axes were accepted");
    requireThrows<std::invalid_argument>([] {
        (void)OpticalRenderer({{0, 0, 0}, {1, 0, 0}, {0, 0, 1}, 180, 1, 1});
    }, "invalid FOV was accepted");
    requireThrows<std::invalid_argument>([&] {
        auto copy = volume;
        copy.setEmissivity("red", 0, 0, 0, -1.0);
    }, "negative emissivity was accepted");
}

} // namespace

int main() {
    try {
        checkAnalyticRays();
        checkPreviewAndValidation();
        std::cout << "OpticalRenderer checks passed\n";
    } catch (const std::exception& error) {
        std::cerr << "OpticalRenderer check failed: " << error.what() << '\n';
        return 1;
    }
}
