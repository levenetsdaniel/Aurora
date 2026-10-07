#ifndef AURORA_VTKELECTRONSCENE_HPP
#define AURORA_VTKELECTRONSCENE_HPP

#include "ElectronSceneData.hpp"
#include "ElectronAnimation.hpp"
#include "../emission/VolumeEmissionGrid.hpp"

#include <vtkActor.h>
#include <vtkCamera.h>
#include <vtkCellArray.h>
#include <vtkColorTransferFunction.h>
#include <vtkGlyph3D.h>
#include <vtkImageData.h>
#include <vtkNew.h>
#include <vtkPNGWriter.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkPolyLine.h>
#include <vtkPiecewiseFunction.h>
#include <vtkProperty.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkSmartPointer.h>
#include <vtkSphereSource.h>
#include <vtkSmartVolumeMapper.h>
#include <vtkTextActor.h>
#include <vtkTextProperty.h>
#include <vtkTubeFilter.h>
#include <vtkVolume.h>
#include <vtkVolumeProperty.h>
#include <vtkWindowToImageFilter.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <filesystem>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

// VTK coordinates are Earth radii. The source data remain in SI metres.
class VtkElectronScene {
private:
    vtkSmartPointer<vtkRenderer> renderer;
    vtkSmartPointer<vtkRenderWindow> window;
    bool animated{};
    std::array<vtkSmartPointer<vtkPolyData>, 4> animationTrails;
    std::array<vtkSmartPointer<vtkPolyData>, 4> animationMarkers;
    vtkSmartPointer<vtkTextActor> animationClock;
    struct EmissionOpacity {
        vtkSmartPointer<vtkPiecewiseFunction> function;
        double maximum{};
        double relativeOpacity{};
    };
    std::vector<EmissionOpacity> emissionOpacities;

    static constexpr std::array<std::array<double, 3>, 4> pathColors{{
        {{1.0, 0.78, 0.20}}, {{0.85, 0.36, 0.95}},
        {{0.60, 0.80, 1.0}}, {{0.55, 0.58, 0.65}}
    }};

    static std::size_t pathCategory(ParticleOutcome outcome) {
        switch (outcome) {
            case ParticleOutcome::Precipitated: return 0;
            case ParticleOutcome::Mirrored: return 1;
            case ParticleOutcome::Escaped: return 2;
            case ParticleOutcome::Unresolved: return 3;
        }
        throw std::invalid_argument("VtkElectronScene: invalid outcome");
    }

    static void setOpacity(EmissionOpacity& entry, double reveal) {
        entry.function->RemoveAllPoints();
        entry.function->AddPoint(0.0, 0.0);
        entry.function->AddPoint(entry.maximum * 0.001, 0.0);
        entry.function->AddPoint(entry.maximum * 0.05,
                                 0.09 * entry.relativeOpacity * reveal);
        entry.function->AddPoint(entry.maximum,
                                 0.55 * entry.relativeOpacity * reveal);
        entry.function->Modified();
    }

    static void addLine(vtkPoints* points, vtkCellArray* cells,
                        const std::vector<Vector3>& positions_m) {
        if (positions_m.size() < 2) return;
        vtkNew<vtkPolyLine> polyline;
        polyline->GetPointIds()->SetNumberOfIds(
            static_cast<vtkIdType>(positions_m.size()));
        for (std::size_t index = 0; index < positions_m.size(); ++index) {
            const Vector3& position = positions_m[index];
            if (!position.isFinite()) {
                throw std::invalid_argument("VtkElectronScene: nonfinite path point");
            }
            const vtkIdType id = points->InsertNextPoint(
                position.x / C::earthRadius_m,
                position.y / C::earthRadius_m,
                position.z / C::earthRadius_m);
            polyline->GetPointIds()->SetId(static_cast<vtkIdType>(index), id);
        }
        cells->InsertNextCell(polyline);
    }

    void addSphere(double radius, double red, double green, double blue,
                   double opacity) {
        vtkNew<vtkSphereSource> sphere;
        sphere->SetRadius(radius);
        sphere->SetThetaResolution(128);
        sphere->SetPhiResolution(96);
        vtkNew<vtkPolyDataMapper> mapper;
        mapper->SetInputConnection(sphere->GetOutputPort());
        vtkNew<vtkActor> actor;
        actor->SetMapper(mapper);
        actor->GetProperty()->SetColor(red, green, blue);
        actor->GetProperty()->SetOpacity(opacity);
        actor->GetProperty()->SetAmbient(0.25);
        actor->GetProperty()->SetDiffuse(0.75);
        actor->GetProperty()->SetSpecular(0.35);
        actor->GetProperty()->SetSpecularPower(25.0);
        renderer->AddActor(actor);
    }

    void addMarkers(const std::vector<Vector3>& positions_m,
                    double radius_Re, double red, double green, double blue,
                    double opacity = 1.0) {
        if (positions_m.empty()) return;
        vtkNew<vtkPoints> points;
        vtkNew<vtkCellArray> vertices;
        for (const Vector3& position : positions_m) {
            if (!position.isFinite()) {
                throw std::invalid_argument("VtkElectronScene: nonfinite marker point");
            }
            const vtkIdType id = points->InsertNextPoint(
                position.x / C::earthRadius_m,
                position.y / C::earthRadius_m,
                position.z / C::earthRadius_m);
            vertices->InsertNextCell(1, &id);
        }
        vtkNew<vtkPolyData> data;
        data->SetPoints(points);
        data->SetVerts(vertices);
        vtkNew<vtkSphereSource> sphere;
        sphere->SetRadius(radius_Re);
        sphere->SetThetaResolution(10);
        sphere->SetPhiResolution(10);
        vtkNew<vtkGlyph3D> glyph;
        glyph->SetInputData(data);
        glyph->SetSourceConnection(sphere->GetOutputPort());
        glyph->ScalingOff();
        vtkNew<vtkPolyDataMapper> mapper;
        mapper->SetInputConnection(glyph->GetOutputPort());
        vtkNew<vtkActor> actor;
        actor->SetMapper(mapper);
        actor->GetProperty()->SetColor(red, green, blue);
        actor->GetProperty()->SetOpacity(opacity);
        actor->GetProperty()->SetAmbient(0.6);
        renderer->AddActor(actor);
    }

    void addMagneticLines(const ElectronSceneData& scene) {
        vtkNew<vtkPoints> points;
        vtkNew<vtkCellArray> lines;
        for (const auto& path : scene.magneticLines_m) addLine(points, lines, path);
        if (lines->GetNumberOfCells() == 0) return;
        vtkNew<vtkPolyData> data;
        data->SetPoints(points);
        data->SetLines(lines);
        vtkNew<vtkTubeFilter> tubes;
        tubes->SetInputData(data);
        tubes->SetRadius(0.00065);
        tubes->SetNumberOfSides(5);
        vtkNew<vtkPolyDataMapper> mapper;
        mapper->SetInputConnection(tubes->GetOutputPort());
        vtkNew<vtkActor> actor;
        actor->SetMapper(mapper);
        actor->GetProperty()->SetColor(0.25, 0.48, 0.75);
        actor->GetProperty()->SetOpacity(0.55);
        actor->GetProperty()->LightingOff();
        renderer->AddActor(actor);
    }

    void addElectronPaths(const ElectronSceneData& scene) {
        std::vector<Vector3> sources;
        std::vector<Vector3> impacts;
        std::vector<Vector3> mirrors;
        const std::size_t markerStride =
            std::max<std::size_t>(1, (scene.paths.size() + 199) / 200);
        for (std::size_t pathIndex = 0; pathIndex < scene.paths.size(); ++pathIndex) {
            const ElectronPath& path = scene.paths[pathIndex];
            if (path.positions_m.empty() || !std::isfinite(path.energy_eV) ||
                path.energy_eV <= 0.0 || !std::isfinite(path.weight_s1) ||
                path.weight_s1 < 0.0) {
                throw std::invalid_argument("VtkElectronScene: invalid electron path");
            }
            if (pathIndex % markerStride == 0) {
                sources.push_back(path.positions_m.front());
                if (path.outcome == ParticleOutcome::Precipitated) {
                    impacts.push_back(path.positions_m.back());
                } else if (path.outcome == ParticleOutcome::Mirrored) {
                    mirrors.push_back(path.positions_m.back());
                }
            }
        }

        const std::array<ParticleOutcome, 4> outcomes{
            ParticleOutcome::Precipitated, ParticleOutcome::Mirrored,
            ParticleOutcome::Escaped, ParticleOutcome::Unresolved};
        for (std::size_t category = 0; category < outcomes.size(); ++category) {
            vtkNew<vtkPoints> points;
            vtkNew<vtkCellArray> lines;
            for (const ElectronPath& path : scene.paths) {
                if (path.outcome == outcomes[category]) {
                    addLine(points, lines, path.positions_m);
                }
            }
            if (lines->GetNumberOfCells() == 0) continue;
            vtkNew<vtkPolyData> data;
            data->SetPoints(points);
            data->SetLines(lines);
            vtkNew<vtkTubeFilter> tubes;
            tubes->SetInputData(data);
            tubes->SetRadius(scene.paths.size() > 300 ? 0.00055 : 0.0012);
            tubes->SetNumberOfSides(6);
            vtkNew<vtkPolyDataMapper> mapper;
            mapper->SetInputConnection(tubes->GetOutputPort());
            vtkNew<vtkActor> actor;
            actor->SetMapper(mapper);
            actor->GetProperty()->SetAmbient(0.5);
            actor->GetProperty()->SetOpacity(scene.paths.size() > 300 ? 0.42 : 1.0);
            actor->GetProperty()->SetColor(pathColors[category][0],
                                           pathColors[category][1],
                                           pathColors[category][2]);
            renderer->AddActor(actor);
        }
        const bool dense = scene.paths.size() > 300;
        addMarkers(sources, dense ? 0.0012 : 0.0022,
                   0.25, 0.95, 1.0, dense ? 0.65 : 1.0);
        addMarkers(impacts, dense ? 0.0016 : 0.0032,
                   1.0, 0.5, 0.16, dense ? 0.65 : 1.0);
        addMarkers(mirrors, dense ? 0.0016 : 0.0032,
                   0.95, 0.45, 1.0, dense ? 0.65 : 1.0);
    }

    void addAnimatedActors(const ElectronSceneData& scene) {
        for (const auto& path : scene.paths) ElectronAnimation::validate(path);
        for (std::size_t category = 0; category < 4; ++category) {
            animationTrails[category] = vtkSmartPointer<vtkPolyData>::New();
            vtkNew<vtkTubeFilter> tubes;
            tubes->SetInputData(animationTrails[category]);
            tubes->SetRadius(scene.paths.size() > 300 ? 0.00055 : 0.0012);
            tubes->SetNumberOfSides(6);
            vtkNew<vtkPolyDataMapper> mapper;
            mapper->SetInputConnection(tubes->GetOutputPort());
            vtkNew<vtkActor> actor;
            actor->SetMapper(mapper);
            actor->GetProperty()->SetColor(pathColors[category][0],
                                           pathColors[category][1],
                                           pathColors[category][2]);
            actor->GetProperty()->SetAmbient(0.5);
            actor->GetProperty()->SetOpacity(scene.paths.size() > 300 ? 0.55 : 1.0);
            renderer->AddActor(actor);

            animationMarkers[category] = vtkSmartPointer<vtkPolyData>::New();
            vtkNew<vtkSphereSource> sphere;
            sphere->SetRadius(scene.paths.size() > 300 ? 0.0022 : 0.0035);
            sphere->SetThetaResolution(8);
            sphere->SetPhiResolution(8);
            vtkNew<vtkGlyph3D> glyph;
            glyph->SetInputData(animationMarkers[category]);
            glyph->SetSourceConnection(sphere->GetOutputPort());
            glyph->ScalingOff();
            vtkNew<vtkPolyDataMapper> pointMapper;
            pointMapper->SetInputConnection(glyph->GetOutputPort());
            vtkNew<vtkActor> pointActor;
            pointActor->SetMapper(pointMapper);
            const std::array<std::array<double, 3>, 4> markerColors{{
                {{0.25, 0.95, 1.0}}, {{1.0, 0.5, 0.16}},
                {{0.95, 0.45, 1.0}}, {{0.60, 0.80, 1.0}}
            }};
            pointActor->GetProperty()->SetColor(markerColors[category][0],
                                                markerColors[category][1],
                                                markerColors[category][2]);
            pointActor->GetProperty()->SetAmbient(0.6);
            renderer->AddActor(pointActor);
        }
        animationClock = vtkSmartPointer<vtkTextActor>::New();
        animationClock->SetDisplayPosition(24, 55);
        animationClock->GetTextProperty()->SetFontSize(
            std::min(18, window->GetSize()[0] / 45));
        animationClock->GetTextProperty()->SetColor(0.9, 0.96, 1.0);
        renderer->AddViewProp(animationClock);
    }

public:
    VtkElectronScene(const ElectronSceneData& scene, int width = 1600,
                     int height = 1000, bool animate = false) : animated(animate) {
        if (!std::isfinite(scene.atmosphereHeight_m) ||
            scene.atmosphereHeight_m < 0.0 || scene.paths.empty() ||
            width <= 0 || height <= 0) {
            throw std::invalid_argument("VtkElectronScene: invalid scene configuration");
        }
        renderer = vtkSmartPointer<vtkRenderer>::New();
        window = vtkSmartPointer<vtkRenderWindow>::New();
        window->AddRenderer(renderer);
        window->SetSize(width, height);
        window->SetWindowName("Aurora: electron precipitation");
        renderer->SetBackground(0.015, 0.025, 0.06);
        renderer->SetBackground2(0.055, 0.10, 0.20);
        renderer->GradientBackgroundOn();
        addSphere(1.0, 0.05, 0.24, 0.39, 1.0);
        addSphere(1.0 + scene.atmosphereHeight_m / C::earthRadius_m,
                  0.20, 0.75, 0.95, 0.12);
        addMagneticLines(scene);
        if (animated) addAnimatedActors(scene);
        else addElectronPaths(scene);

        vtkNew<vtkTextActor> title;
        const std::string caption = animated ?
            "ELECTRON FLIGHT  |  " + std::to_string(scene.paths.size()) + " markers" :
            "ELECTRON PRECIPITATION  |  " +
            std::to_string(scene.paths.size()) + " markers  |  " +
            std::to_string(scene.precipitationCount()) + " atmosphere entries  |  " +
            std::to_string(scene.countOutcome(ParticleOutcome::Mirrored)) + " mirrored";
        title->SetInput(caption.c_str());
        title->SetDisplayPosition(24, height - 48);
        title->GetTextProperty()->SetFontSize(std::min(22, width / 32));
        title->GetTextProperty()->SetColor(0.93, 0.97, 1.0);
        title->GetTextProperty()->BoldOn();
        renderer->AddViewProp(title);
        vtkNew<vtkTextActor> key;
        key->SetInput(animated ?
            "Gold: entry trail   Magenta: mirror trail   Cyan: flying   Orange: entered" :
            "Gold: precipitated   Magenta: mirrored   Cyan: injection   Orange: atmosphere entry   Blue: field");
        key->SetDisplayPosition(24, 22);
        key->GetTextProperty()->SetFontSize(std::min(16, width / 50));
        key->GetTextProperty()->SetColor(0.77, 0.87, 0.98);
        renderer->AddViewProp(key);

        Vector3 focal_m;
        for (const ElectronPath& path : scene.paths) focal_m = focal_m + path.positions_m.front();
        focal_m = focal_m / static_cast<double>(scene.paths.size());
        const Vector3 focus = focal_m / C::earthRadius_m;
        vtkCamera* camera = renderer->GetActiveCamera();
        camera->SetFocalPoint(focus.x, focus.y, focus.z - 0.025);
        camera->SetPosition(focus.x + 0.7, focus.y - 1.1, focus.z + 0.55);
        camera->SetViewUp(0.0, 0.0, 1.0);
        camera->ParallelProjectionOn();
        camera->SetParallelScale(0.40);
        renderer->ResetCameraClippingRange();
    }

    vtkRenderer* getRenderer() const { return renderer; }

    void setCameraScale(double scale_Re) {
        if (!std::isfinite(scale_Re) || scale_Re <= 0.0) {
            throw std::invalid_argument("VtkElectronScene: invalid camera scale");
        }
        renderer->GetActiveCamera()->SetParallelScale(scale_Re);
        renderer->ResetCameraClippingRange();
    }

    void setAnimationFrame(const ElectronSceneData& scene, double time_s) {
        if (!animated || !std::isfinite(time_s) || time_s < 0.0)
            throw std::invalid_argument("VtkElectronScene: invalid animation frame");
        std::array<vtkSmartPointer<vtkPoints>, 4> trailPoints;
        std::array<vtkSmartPointer<vtkCellArray>, 4> trailCells;
        std::array<vtkSmartPointer<vtkPoints>, 4> markerPoints;
        std::array<vtkSmartPointer<vtkCellArray>, 4> markerCells;
        for (std::size_t category = 0; category < 4; ++category) {
            trailPoints[category] = vtkSmartPointer<vtkPoints>::New();
            trailCells[category] = vtkSmartPointer<vtkCellArray>::New();
            markerPoints[category] = vtkSmartPointer<vtkPoints>::New();
            markerCells[category] = vtkSmartPointer<vtkCellArray>::New();
        }
        std::size_t flying = 0;
        std::size_t precipitated = 0;
        std::size_t mirrored = 0;
        for (const auto& path : scene.paths) {
            const auto state = ElectronAnimation::sample(path, time_s);
            if (!state.visible) continue;
            const std::size_t category = pathCategory(path.outcome);
            std::vector<Vector3> trail(path.positions_m.begin(),
                                       path.positions_m.begin() +
                                           static_cast<std::ptrdiff_t>(state.completedPoints));
            if (!state.finished &&
                (trail.back() - state.position_m).norm() > 0.0)
                trail.push_back(state.position_m);
            addLine(trailPoints[category], trailCells[category], trail);
            const std::size_t markerCategory = !state.finished ? 0 :
                path.outcome == ParticleOutcome::Precipitated ? 1 :
                path.outcome == ParticleOutcome::Mirrored ? 2 : 3;
            const vtkIdType pointId = markerPoints[markerCategory]->InsertNextPoint(
                state.position_m.x / C::earthRadius_m,
                state.position_m.y / C::earthRadius_m,
                state.position_m.z / C::earthRadius_m);
            markerCells[markerCategory]->InsertNextCell(1, &pointId);
            if (!state.finished) ++flying;
            else if (path.outcome == ParticleOutcome::Precipitated) ++precipitated;
            else if (path.outcome == ParticleOutcome::Mirrored) ++mirrored;
        }
        for (std::size_t category = 0; category < 4; ++category) {
            animationTrails[category]->SetPoints(trailPoints[category]);
            animationTrails[category]->SetLines(trailCells[category]);
            animationTrails[category]->Modified();
            animationMarkers[category]->SetPoints(markerPoints[category]);
            animationMarkers[category]->SetVerts(markerCells[category]);
            animationMarkers[category]->Modified();
        }
        const std::string clock = "t = " + std::to_string(time_s * 1e6) +
            " us  |  flying " + std::to_string(flying) +
            "  |  entered " + std::to_string(precipitated) +
            "  |  mirrored " + std::to_string(mirrored);
        animationClock->SetInput(clock.c_str());
        setEmissionReveal(ElectronAnimation::arrivedWeightFraction(scene, time_s));
    }

    void setEmissionReveal(double fraction) {
        if (!std::isfinite(fraction) || fraction < 0.0 || fraction > 1.0)
            throw std::invalid_argument("VtkElectronScene: invalid emission reveal");
        for (auto& entry : emissionOpacities) setOpacity(entry, fraction);
    }

    void addEmission(const VolumeEmissionGrid& emission,
                     const std::string& lineId,
                     double red, double green, double blue,
                     double relativeOpacity = 1.0) {
        const auto& grid = emission.configuration();
        if (!std::isfinite(relativeOpacity) || relativeOpacity < 0.0 ||
            relativeOpacity > 1.0) {
            throw std::invalid_argument("VtkElectronScene: invalid relative opacity");
        }
        if (grid.cellsX > static_cast<std::size_t>(std::numeric_limits<int>::max()) ||
            grid.cellsY > static_cast<std::size_t>(std::numeric_limits<int>::max()) ||
            grid.cellsZ > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
            throw std::invalid_argument("VtkElectronScene: emission grid too large");
        }
        vtkNew<vtkImageData> image;
        image->SetDimensions(static_cast<int>(grid.cellsX),
                             static_cast<int>(grid.cellsY),
                             static_cast<int>(grid.cellsZ));
        image->SetOrigin(
            (grid.origin_m.x + grid.cellSize_m.x / 2) / C::earthRadius_m,
            (grid.origin_m.y + grid.cellSize_m.y / 2) / C::earthRadius_m,
            (grid.origin_m.z + grid.cellSize_m.z / 2) / C::earthRadius_m);
        image->SetSpacing(grid.cellSize_m.x / C::earthRadius_m,
                          grid.cellSize_m.y / C::earthRadius_m,
                          grid.cellSize_m.z / C::earthRadius_m);
        image->AllocateScalars(VTK_DOUBLE, 1);
        const std::size_t count = grid.cellsX * grid.cellsY * grid.cellsZ;
        std::vector<double> display(count, 0.0);
        const auto flat = [&](std::size_t x, std::size_t y, std::size_t z) {
            return (z * grid.cellsY + y) * grid.cellsX + x;
        };
        // Photon-conserving kernel smoothing is for the display only. The
        // exported emission grid retains the unsmoothed physical cell values.
        for (std::size_t z = 0; z < grid.cellsZ; ++z) {
            for (std::size_t y = 0; y < grid.cellsY; ++y) {
                for (std::size_t x = 0; x < grid.cellsX; ++x) {
                    const double value = emission.emissivity(lineId, x, y, z);
                    if (value == 0.0) continue;
                    double sum = 0.0;
                    for (int dz = -1; dz <= 1; ++dz) {
                        for (int dy = -2; dy <= 2; ++dy) {
                            for (int dx = -2; dx <= 2; ++dx) {
                                const long long nx = static_cast<long long>(x) + dx;
                                const long long ny = static_cast<long long>(y) + dy;
                                const long long nz = static_cast<long long>(z) + dz;
                                if (nx < 0 || ny < 0 || nz < 0 ||
                                    nx >= static_cast<long long>(grid.cellsX) ||
                                    ny >= static_cast<long long>(grid.cellsY) ||
                                    nz >= static_cast<long long>(grid.cellsZ)) continue;
                                const double radius2 =
                                    (dx * dx + dy * dy) / (1.4 * 1.4) + dz * dz;
                                sum += std::exp(-0.5 * radius2);
                            }
                        }
                    }
                    for (int dz = -1; dz <= 1; ++dz) {
                        for (int dy = -2; dy <= 2; ++dy) {
                            for (int dx = -2; dx <= 2; ++dx) {
                                const long long nx = static_cast<long long>(x) + dx;
                                const long long ny = static_cast<long long>(y) + dy;
                                const long long nz = static_cast<long long>(z) + dz;
                                if (nx < 0 || ny < 0 || nz < 0 ||
                                    nx >= static_cast<long long>(grid.cellsX) ||
                                    ny >= static_cast<long long>(grid.cellsY) ||
                                    nz >= static_cast<long long>(grid.cellsZ)) continue;
                                const double radius2 =
                                    (dx * dx + dy * dy) / (1.4 * 1.4) + dz * dz;
                                display[flat(static_cast<std::size_t>(nx),
                                             static_cast<std::size_t>(ny),
                                             static_cast<std::size_t>(nz))] +=
                                    value * std::exp(-0.5 * radius2) / sum;
                            }
                        }
                    }
                }
            }
        }
        double maximum = 0.0;
        for (std::size_t z = 0; z < grid.cellsZ; ++z) {
            for (std::size_t y = 0; y < grid.cellsY; ++y) {
                for (std::size_t x = 0; x < grid.cellsX; ++x) {
                    const double value = display[flat(x, y, z)];
                    *static_cast<double*>(image->GetScalarPointer(
                        static_cast<int>(x), static_cast<int>(y), static_cast<int>(z))) = value;
                    maximum = std::max(maximum, value);
                }
            }
        }
        if (maximum == 0.0) return;
        vtkNew<vtkSmartVolumeMapper> mapper;
        mapper->SetInputData(image);
        mapper->SetBlendModeToComposite();
        vtkNew<vtkColorTransferFunction> colors;
        colors->AddRGBPoint(0.0, 0.0, 0.0, 0.0);
        colors->AddRGBPoint(maximum * 0.02, red * 0.35, green * 0.35, blue * 0.35);
        colors->AddRGBPoint(maximum, red, green, blue);
        vtkNew<vtkPiecewiseFunction> opacity;
        EmissionOpacity entry{opacity, maximum, relativeOpacity};
        setOpacity(entry, 1.0);
        vtkNew<vtkVolumeProperty> property;
        property->SetColor(colors);
        property->SetScalarOpacity(opacity);
        property->SetScalarOpacityUnitDistance(
            std::min({grid.cellSize_m.x, grid.cellSize_m.y, grid.cellSize_m.z}) /
            C::earthRadius_m);
        property->SetInterpolationTypeToLinear();
        property->ShadeOff();
        vtkNew<vtkVolume> volume;
        volume->SetMapper(mapper);
        volume->SetProperty(property);
        renderer->AddVolume(volume);
        emissionOpacities.push_back(entry);
        renderer->ResetCameraClippingRange();
    }

    void savePng(const std::filesystem::path& path) {
        window->SetOffScreenRendering(1);
        window->Render();
        vtkNew<vtkWindowToImageFilter> capture;
        capture->SetInput(window);
        capture->SetInputBufferTypeToRGB();
        capture->ReadFrontBufferOff();
        capture->Update();
        vtkNew<vtkPNGWriter> writer;
        writer->SetFileName(path.string().c_str());
        writer->SetInputConnection(capture->GetOutputPort());
        writer->Write();
        if (writer->GetErrorCode() != 0) {
            throw std::runtime_error("VtkElectronScene: failed to write PNG");
        }
    }

    void interact() {
        window->SetOffScreenRendering(0);
        vtkNew<vtkRenderWindowInteractor> interactor;
        interactor->SetRenderWindow(window);
        window->Render();
        interactor->Initialize();
        interactor->Start();
    }
};

#endif // AURORA_VTKELECTRONSCENE_HPP
