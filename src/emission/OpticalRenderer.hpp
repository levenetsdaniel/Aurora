#ifndef AURORA_OPTICALRENDERER_HPP
#define AURORA_OPTICALRENDERER_HPP

#include "VolumeEmissionGrid.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <numbers>
#include <ostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

struct OpticalCamera {
    Vector3 position_m;
    Vector3 forward;
    Vector3 up;
    double verticalFov_deg{};
    std::size_t width{};
    std::size_t height{};
};

struct SpectralImage {
    std::string lineId;
    double wavelength_nm{};
    std::size_t width{};
    std::size_t height{};
    // Pixel-center radiance: photons m^-2 s^-1 sr^-1.
    std::vector<double> radiance_ph_m2_s1_sr1;

    double at(std::size_t x, std::size_t y) const {
        if (x >= width || y >= height) throw std::out_of_range("SpectralImage: pixel outside image");
        return radiance_ph_m2_s1_sr1.at(y * width + x);
    }

    void writeCsv(std::ostream& output) const {
        if (width == 0 || height == 0 ||
            width > std::numeric_limits<std::size_t>::max() / height ||
            radiance_ph_m2_s1_sr1.size() != width * height) {
            throw std::invalid_argument("SpectralImage: invalid image dimensions");
        }
        output << "# line_id=" << lineId << ",wavelength_nm=" << wavelength_nm << '\n';
        output << "x,y,I_ph_m2_s1_sr1\n" << std::setprecision(17);
        for (std::size_t y = 0; y < height; ++y) {
            for (std::size_t x = 0; x < width; ++x) {
                output << x << ',' << y << ',' << at(x, y) << '\n';
            }
        }
        if (!output) throw std::runtime_error("SpectralImage: write failed");
    }

    void saveCsv(const std::filesystem::path& path) const {
        std::ofstream output(path);
        if (!output) throw std::runtime_error("SpectralImage: cannot open " + path.string());
        writeCsv(output);
    }
};

struct RgbPreview {
    std::size_t width{};
    std::size_t height{};
    std::vector<unsigned char> rgb;

    void writePpm(std::ostream& output) const {
        if (width == 0 || height == 0 ||
            width > std::numeric_limits<std::size_t>::max() / height ||
            width * height > std::numeric_limits<std::size_t>::max() / 3 ||
            rgb.size() != width * height * 3) {
            throw std::invalid_argument("RgbPreview: invalid image dimensions");
        }
        output << "P6\n" << width << ' ' << height << "\n255\n";
        output.write(reinterpret_cast<const char*>(rgb.data()),
                     static_cast<std::streamsize>(rgb.size()));
        if (!output) throw std::runtime_error("RgbPreview: write failed");
    }

    void savePpm(const std::filesystem::path& path) const {
        std::ofstream output(path, std::ios::binary);
        if (!output) throw std::runtime_error("RgbPreview: cannot open " + path.string());
        writePpm(output);
    }
};

class OpticalRenderer {
private:
    OpticalCamera camera;
    Vector3 forward;
    Vector3 right;
    Vector3 up;
    double halfHeight{};
    double halfWidth{};

    static Vector3 unit(const Vector3& value) {
        const double length = value.norm();
        if (!value.isFinite() || !std::isfinite(length) || length <= 0.0) {
            throw std::invalid_argument("OpticalRenderer: invalid camera basis");
        }
        return value / length;
    }

    static bool intersectBox(const Vector3& origin, const Vector3& direction,
                             const Vector3& lower, const Vector3& upper,
                             double& enter_m, double& exit_m) {
        double near_m = -std::numeric_limits<double>::infinity();
        double far_m = std::numeric_limits<double>::infinity();
        for (int axis = 0; axis < 3; ++axis) {
            const double position = axis == 0 ? origin.x : axis == 1 ? origin.y : origin.z;
            const double speed = axis == 0 ? direction.x : axis == 1 ? direction.y : direction.z;
            const double low = axis == 0 ? lower.x : axis == 1 ? lower.y : lower.z;
            const double high = axis == 0 ? upper.x : axis == 1 ? upper.y : upper.z;
            if (speed == 0.0) {
                if (position < low || position > high) return false;
                continue;
            }
            const double first = (low - position) / speed;
            const double second = (high - position) / speed;
            near_m = std::max(near_m, std::min(first, second));
            far_m = std::min(far_m, std::max(first, second));
        }
        enter_m = std::max(0.0, near_m);
        exit_m = far_m;
        return std::isfinite(enter_m) && std::isfinite(exit_m) && exit_m > enter_m;
    }

    static std::size_t cellIndex(double coordinate, double lower, double cellSize,
                                 std::size_t count) {
        const double raw = std::floor((coordinate - lower) / cellSize);
        if (raw <= 0.0) return 0;
        if (raw >= static_cast<double>(count)) return count - 1;
        return static_cast<std::size_t>(raw);
    }

    static void appendPlanes(std::vector<double>& points, double origin,
                             double direction, double lower, double cellSize,
                             std::size_t count, double enter, double exit) {
        if (direction == 0.0) return;
        for (std::size_t boundary = 1; boundary < count; ++boundary) {
            const double plane = lower + static_cast<double>(boundary) * cellSize;
            const double distance = (plane - origin) / direction;
            if (distance > enter && distance < exit) points.push_back(distance);
        }
    }

public:
    explicit OpticalRenderer(OpticalCamera config) : camera(config) {
        if (!camera.position_m.isFinite() || !camera.forward.isFinite() ||
            !camera.up.isFinite() || !std::isfinite(camera.verticalFov_deg) ||
            camera.verticalFov_deg <= 0.0 || camera.verticalFov_deg >= 180.0 ||
            camera.width == 0 || camera.height == 0 ||
            camera.width > std::numeric_limits<std::size_t>::max() / camera.height) {
            throw std::invalid_argument("OpticalRenderer: invalid camera");
        }
        forward = unit(camera.forward);
        right = unit(camera.up.cross(forward));
        up = unit(forward.cross(right));
        halfHeight = std::tan(camera.verticalFov_deg * std::numbers::pi / 360.0);
        halfWidth = halfHeight * static_cast<double>(camera.width) / camera.height;
        if (!std::isfinite(halfWidth)) {
            throw std::invalid_argument("OpticalRenderer: invalid field of view");
        }
    }

    Vector3 pixelDirection(std::size_t x, std::size_t y) const {
        if (x >= camera.width || y >= camera.height) {
            throw std::out_of_range("OpticalRenderer: pixel outside image");
        }
        const double horizontal = (2.0 * (static_cast<double>(x) + 0.5) /
                                   static_cast<double>(camera.width) - 1.0) * halfWidth;
        const double vertical = (1.0 - 2.0 * (static_cast<double>(y) + 0.5) /
                                 static_cast<double>(camera.height)) * halfHeight;
        return unit(forward + right * horizontal + up * vertical);
    }

    std::vector<SpectralImage> render(const VolumeEmissionGrid& volume) const {
        std::vector<SpectralImage> images;
        for (const SpectralLine& line : volume.lines()) {
            images.push_back({line.id, line.wavelength_nm, camera.width, camera.height,
                              std::vector<double>(camera.width * camera.height, 0.0)});
        }
        const auto& grid = volume.configuration();
        for (std::size_t y = 0; y < camera.height; ++y) {
            for (std::size_t x = 0; x < camera.width; ++x) {
                const Vector3 direction = pixelDirection(x, y);
                double enter_m{}, exit_m{};
                if (!intersectBox(camera.position_m, direction, grid.origin_m,
                                  volume.upperCorner_m(), enter_m, exit_m)) continue;
                std::vector<double> points{enter_m, exit_m};
                appendPlanes(points, camera.position_m.x, direction.x, grid.origin_m.x,
                             grid.cellSize_m.x, grid.cellsX, enter_m, exit_m);
                appendPlanes(points, camera.position_m.y, direction.y, grid.origin_m.y,
                             grid.cellSize_m.y, grid.cellsY, enter_m, exit_m);
                appendPlanes(points, camera.position_m.z, direction.z, grid.origin_m.z,
                             grid.cellSize_m.z, grid.cellsZ, enter_m, exit_m);
                std::sort(points.begin(), points.end());
                for (std::size_t segment = 1; segment < points.size(); ++segment) {
                    const double length_m = points[segment] - points[segment - 1];
                    if (length_m <= 0.0) continue;
                    const double midpoint_m = points[segment - 1] + length_m / 2.0;
                    const Vector3 center = camera.position_m + direction * midpoint_m;
                    const auto cx = cellIndex(center.x, grid.origin_m.x,
                                              grid.cellSize_m.x, grid.cellsX);
                    const auto cy = cellIndex(center.y, grid.origin_m.y,
                                              grid.cellSize_m.y, grid.cellsY);
                    const auto cz = cellIndex(center.z, grid.origin_m.z,
                                              grid.cellSize_m.z, grid.cellsZ);
                    for (SpectralImage& image : images) {
                        double& pixel = image.radiance_ph_m2_s1_sr1[y * camera.width + x];
                        pixel += volume.emissivity(image.lineId, cx, cy, cz) *
                                 length_m / (4.0 * std::numbers::pi);
                        if (!std::isfinite(pixel)) {
                            throw std::overflow_error("OpticalRenderer: radiance overflow");
                        }
                    }
                }
            }
        }
        return images;
    }

    // Display-only red/green mapping; never use the PPM values as radiance.
    static RgbPreview makeRgbPreview(const SpectralImage& red,
                                     const SpectralImage& green,
                                     double referenceRadiance_ph_m2_s1_sr1) {
        if (red.width != green.width || red.height != green.height ||
            red.width == 0 || red.height == 0 ||
            red.width > std::numeric_limits<std::size_t>::max() / red.height ||
            red.width * red.height > std::numeric_limits<std::size_t>::max() / 3 ||
            red.radiance_ph_m2_s1_sr1.size() != red.width * red.height ||
            green.radiance_ph_m2_s1_sr1.size() != green.width * green.height ||
            !std::isfinite(referenceRadiance_ph_m2_s1_sr1) ||
            referenceRadiance_ph_m2_s1_sr1 <= 0.0) {
            throw std::invalid_argument("OpticalRenderer: invalid RGB preview inputs");
        }
        RgbPreview result{red.width, red.height,
                          std::vector<unsigned char>(red.width * red.height * 3, 0)};
        for (std::size_t pixel = 0; pixel < red.width * red.height; ++pixel) {
            const double redValue = red.radiance_ph_m2_s1_sr1[pixel];
            const double greenValue = green.radiance_ph_m2_s1_sr1[pixel];
            if (!std::isfinite(redValue) || redValue < 0.0 ||
                !std::isfinite(greenValue) || greenValue < 0.0) {
                throw std::invalid_argument("OpticalRenderer: invalid RGB radiance");
            }
            result.rgb[3 * pixel] = static_cast<unsigned char>(
                std::lround(255.0 * (1.0 - std::exp(-redValue /
                                                referenceRadiance_ph_m2_s1_sr1))));
            result.rgb[3 * pixel + 1] = static_cast<unsigned char>(
                std::lround(255.0 * (1.0 - std::exp(-greenValue /
                                                referenceRadiance_ph_m2_s1_sr1))));
        }
        return result;
    }

    static void saveSpectralCsv(const std::vector<SpectralImage>& images,
                                const std::filesystem::path& directory) {
        if (images.empty()) {
            throw std::invalid_argument("OpticalRenderer: no spectral images to save");
        }
        for (const SpectralImage& image : images) {
            if (image.lineId.empty() ||
                !std::all_of(image.lineId.begin(), image.lineId.end(),
                    [](unsigned char character) {
                        return std::isalnum(character) || character == '_' || character == '-';
                    })) {
                throw std::invalid_argument("OpticalRenderer: invalid line ID for filename");
            }
        }
        std::filesystem::create_directories(directory);
        for (const SpectralImage& image : images) {
            image.saveCsv(directory / (image.lineId + ".csv"));
        }
    }
};

#endif // AURORA_OPTICALRENDERER_HPP
