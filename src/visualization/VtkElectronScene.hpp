#ifndef AURORA_VTKELECTRONSCENE_HPP
#define AURORA_VTKELECTRONSCENE_HPP

#include "ElectronSceneData.hpp"
#include "../emission/VolumeEmissionGrid.hpp"

#include <vtkActor.h>
#include <vtkCamera.h>
#include <vtkCellArray.h>
#include <vtkDoubleArray.h>
#include <vtkColorTransferFunction.h>
#include <vtkGlyph3D.h>
#include <vtkImageData.h>
#include <vtkLookupTable.h>
#include <vtkNew.h>
#include <vtkPNGWriter.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkPolyLine.h>
#include <vtkPiecewiseFunction.h>
#include <vtkProperty.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkScalarBarActor.h>
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
#include <cmath>
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

    static void addLine(vtkPoints* points, vtkCellArray* cells,
                        const std::vector<Vector3>& positions_m,
                        vtkDoubleArray* energy_keV = nullptr,
                        double value_keV = 0.0) {
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
            if (energy_keV) energy_keV->InsertNextValue(value_keV);
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
        vtkNew<vtkPoints> points;
        vtkNew<vtkCellArray> lines;
        vtkNew<vtkDoubleArray> energies;
        energies->SetName("Electron energy (keV)");
        double low_keV = std::numeric_limits<double>::infinity();
        double high_keV = 0.0;
        std::vector<Vector3> sources;
        std::vector<Vector3> impacts;
        const std::size_t markerStride =
            std::max<std::size_t>(1, (scene.paths.size() + 199) / 200);
        for (std::size_t pathIndex = 0; pathIndex < scene.paths.size(); ++pathIndex) {
            const ElectronPath& path = scene.paths[pathIndex];
            if (path.positions_m.empty() || !std::isfinite(path.energy_eV) ||
                path.energy_eV <= 0.0 || !std::isfinite(path.weight_s1) ||
                path.weight_s1 < 0.0) {
                throw std::invalid_argument("VtkElectronScene: invalid electron path");
            }
            const double energy_keV = path.energy_eV / 1000.0;
            low_keV = std::min(low_keV, energy_keV);
            high_keV = std::max(high_keV, energy_keV);
            addLine(points, lines, path.positions_m, energies, energy_keV);
            if (pathIndex % markerStride == 0) {
                sources.push_back(path.positions_m.front());
                if (path.outcome == ParticleOutcome::Precipitated) {
                    impacts.push_back(path.positions_m.back());
                }
            }
        }
        if (lines->GetNumberOfCells() > 0) {
            vtkNew<vtkPolyData> data;
            data->SetPoints(points);
            data->SetLines(lines);
            data->GetPointData()->SetScalars(energies);
            vtkNew<vtkTubeFilter> tubes;
            tubes->SetInputData(data);
            tubes->SetRadius(scene.paths.size() > 300 ? 0.00055 : 0.0012);
            tubes->SetNumberOfSides(6);
            vtkNew<vtkLookupTable> colors;
            colors->SetHueRange(0.52, 0.02);
            colors->SetSaturationRange(0.8, 0.95);
            colors->SetValueRange(0.9, 1.0);
            colors->Build();
            vtkNew<vtkPolyDataMapper> mapper;
            mapper->SetInputConnection(tubes->GetOutputPort());
            mapper->SetLookupTable(colors);
            mapper->SetScalarModeToUsePointData();
            if (low_keV == high_keV) {
                low_keV *= 0.9;
                high_keV *= 1.1;
            }
            mapper->SetScalarRange(low_keV, high_keV);
            vtkNew<vtkActor> actor;
            actor->SetMapper(mapper);
            actor->GetProperty()->SetAmbient(0.5);
            actor->GetProperty()->SetOpacity(scene.paths.size() > 300 ? 0.42 : 1.0);
            renderer->AddActor(actor);
            vtkNew<vtkScalarBarActor> legend;
            legend->SetLookupTable(colors);
            legend->SetTitle("Electron energy (keV)");
            legend->SetNumberOfLabels(4);
            legend->SetWidth(0.09);
            legend->SetHeight(0.35);
            legend->SetPosition(0.87, 0.1);
            legend->GetTitleTextProperty()->SetColor(0.9, 0.95, 1.0);
            legend->GetLabelTextProperty()->SetColor(0.9, 0.95, 1.0);
            renderer->AddViewProp(legend);
        }
        const bool dense = scene.paths.size() > 300;
        addMarkers(sources, dense ? 0.0012 : 0.0022,
                   0.25, 0.95, 1.0, dense ? 0.65 : 1.0);
        addMarkers(impacts, dense ? 0.0016 : 0.0032,
                   1.0, 0.5, 0.16, dense ? 0.65 : 1.0);
    }

public:
    VtkElectronScene(const ElectronSceneData& scene, int width = 1600,
                     int height = 1000) {
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
        addElectronPaths(scene);

        vtkNew<vtkTextActor> title;
        const std::string caption = "ELECTRON PRECIPITATION  |  " +
            std::to_string(scene.paths.size()) + " markers  |  " +
            std::to_string(scene.precipitationCount()) + " impacts";
        title->SetInput(caption.c_str());
        title->SetDisplayPosition(24, height - 48);
        title->GetTextProperty()->SetFontSize(22);
        title->GetTextProperty()->SetColor(0.93, 0.97, 1.0);
        title->GetTextProperty()->BoldOn();
        renderer->AddViewProp(title);
        vtkNew<vtkTextActor> key;
        key->SetInput("Cyan: injection   Orange: atmosphere entry   Blue: dipole field");
        key->SetDisplayPosition(24, 22);
        key->GetTextProperty()->SetFontSize(16);
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
        opacity->AddPoint(0.0, 0.0);
        opacity->AddPoint(maximum * 0.001, 0.0);
        opacity->AddPoint(maximum * 0.05, 0.09 * relativeOpacity);
        opacity->AddPoint(maximum, 0.55 * relativeOpacity);
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
