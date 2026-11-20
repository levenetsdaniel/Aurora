#include "aurora/visualization/VtkViewer.hpp"
#include "aurora/physics/Gyromotion.hpp"
#include <vtkActor.h>
#include <vtkAxis.h>
#include <vtkBillboardTextActor3D.h>
#include <vtkBrush.h>
#include <vtkCallbackCommand.h>
#include <vtkCamera.h>
#include <vtkCellArray.h>
#include <vtkChartXY.h>
#include <vtkCommand.h>
#include <vtkContextActor.h>
#include <vtkContextScene.h>
#include <vtkCoordinate.h>
#include <vtkDoubleArray.h>
#include <vtkErrorCode.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkLookupTable.h>
#include <vtkNew.h>
#include <vtkPNGWriter.h>
#include <vtkPen.h>
#include <vtkPlot.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkScalarBarActor.h>
#include <vtkSmartPointer.h>
#include <vtkSphereSource.h>
#include <vtkTable.h>
#include <vtkTextActor.h>
#include <vtkTextProperty.h>
#include <vtkTubeFilter.h>
#include <vtkWindowToImageFilter.h>
#include <algorithm>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <numbers>
#include <sstream>

namespace aurora::visualization {
namespace {
constexpr double R = constants::earth_radius_m;
using Clock = std::chrono::steady_clock;

bool is_control_key(const std::string& key) {
    return key == "1" || key == "2" || key == "3" || key == "space" || key == "r" || key == "R" ||
        key == "Left" || key == "Right" || key == "plus" || key == "equal" || key == "minus" ||
        key == "s" || key == "S" || key == "q" || key == "Q";
}

vtkSmartPointer<vtkLookupTable> colors() {
    auto table = vtkSmartPointer<vtkLookupTable>::New();
    table->SetNumberOfTableValues(256);
    table->SetRange(-1, 1);
    table->Build();
    const double negative[] = {0.96, 0.47, 0.28};
    const double zero[] = {0.95, 0.95, 0.77};
    const double positive[] = {0.20, 0.88, 0.85};
    for (int i = 0; i < 256; ++i) {
        const double value = -1 + 2 * static_cast<double>(i) / 255;
        const double f = std::abs(value);
        const auto* end = value < 0 ? negative : positive;
        table->SetTableValue(i, zero[0] * (1 - f) + end[0] * f, zero[1] * (1 - f) + end[1] * f,
            zero[2] * (1 - f) + end[2] * f, 1);
    }
    return table;
}

vtkSmartPointer<vtkActor> sphere(double radius, double red, double green, double blue,
    double opacity = 1, bool wireframe = false) {
    vtkNew<vtkSphereSource> source;
    source->SetRadius(radius);
    source->SetThetaResolution(wireframe ? 24 : 128);
    source->SetPhiResolution(wireframe ? 14 : 64);
    source->SetLatLongTessellation(wireframe);
    vtkNew<vtkPolyDataMapper> mapper;
    mapper->SetInputConnection(source->GetOutputPort());
    auto actor = vtkSmartPointer<vtkActor>::New();
    actor->SetMapper(mapper);
    actor->GetProperty()->SetColor(red, green, blue);
    actor->GetProperty()->SetOpacity(opacity);
    actor->GetProperty()->SetAmbient(0.32);
    actor->GetProperty()->SetDiffuse(0.68);
    if (wireframe) {
        actor->GetProperty()->SetRepresentationToWireframe();
        actor->GetProperty()->LightingOff();
    }
    return actor;
}

vtkSmartPointer<vtkActor> line_actor(const std::vector<Vec3>& positions,
    const std::vector<double>& scalars, double radius, vtkLookupTable* lookup) {
    vtkNew<vtkPoints> points;
    points->SetDataTypeToDouble();
    vtkNew<vtkCellArray> cells;
    cells->InsertNextCell(static_cast<int>(positions.size()));
    vtkNew<vtkDoubleArray> values;
    values->SetName("p_parallel / |p|");
    for (std::size_t i = 0; i < positions.size(); ++i) {
        points->InsertNextPoint(positions[i].x, positions[i].y, positions[i].z);
        cells->InsertCellPoint(static_cast<vtkIdType>(i));
        values->InsertNextValue(scalars.empty() ? 0 : scalars[i]);
    }
    vtkNew<vtkPolyData> data;
    data->SetPoints(points);
    data->SetLines(cells);
    data->GetPointData()->SetScalars(values);
    vtkNew<vtkPolyDataMapper> mapper;
    if (radius > 0) {
        vtkNew<vtkTubeFilter> tube;
        tube->SetInputData(data);
        tube->SetRadius(radius);
        tube->SetNumberOfSides(8);
        tube->CappingOn();
        mapper->SetInputConnection(tube->GetOutputPort());
    } else
        mapper->SetInputData(data);
    if (lookup) {
        mapper->SetLookupTable(lookup);
        mapper->SetScalarRange(-1, 1);
        mapper->SetScalarModeToUsePointData();
    } else
        mapper->ScalarVisibilityOff();
    auto actor = vtkSmartPointer<vtkActor>::New();
    actor->SetMapper(mapper);
    actor->GetProperty()->SetLineWidth(1.2f);
    actor->GetProperty()->SetAmbient(0.6);
    return actor;
}

vtkSmartPointer<vtkTextActor> label(
    vtkRenderer* renderer, const std::string& text, double x, double y, int size = 15) {
    auto actor = vtkSmartPointer<vtkTextActor>::New();
    actor->SetInput(text.c_str());
    actor->GetPositionCoordinate()->SetCoordinateSystemToNormalizedViewport();
    actor->SetPosition(x, y);
    actor->GetTextProperty()->SetFontFamilyToArial();
    actor->GetTextProperty()->SetFontSize(size);
    actor->GetTextProperty()->SetColor(0.86, 0.91, 0.96);
    actor->GetTextProperty()->SetVerticalJustificationToTop();
    renderer->AddViewProp(actor);
    return actor;
}

vtkSmartPointer<vtkTable> two_point_table(double x0, double x1, double y0, double y1) {
    auto table = vtkSmartPointer<vtkTable>::New();
    vtkNew<vtkDoubleArray> x;
    vtkNew<vtkDoubleArray> y;
    x->SetName("x");
    y->SetName("y");
    x->InsertNextValue(x0);
    x->InsertNextValue(x1);
    y->InsertNextValue(y0);
    y->InsertNextValue(y1);
    table->AddColumn(x);
    table->AddColumn(y);
    return table;
}

struct Chart {
    vtkSmartPointer<vtkChartXY> chart;
    vtkSmartPointer<vtkTable> cursor;
};

Chart make_chart(vtkRenderer* renderer, const DemoCase& demo, bool altitude) {
    vtkNew<vtkContextActor> actor;
    auto chart = vtkSmartPointer<vtkChartXY>::New();
    actor->GetScene()->AddItem(chart);
    renderer->AddActor(actor);
    chart->SetAutoSize(true);
    chart->SetShowLegend(false);
    chart->SetTitle(altitude ? "Altitude above Earth" : "Parallel momentum: reflection at zero");
    chart->GetTitleProperties()->SetColor(0.88, 0.92, 0.97);
    chart->GetTitleProperties()->SetFontSize(16);
    chart->GetBackgroundBrush()->SetColor(8, 17, 27, 255);
    const double duration_ms = demo.trajectory.back().time_s * 1000;
    double maximum_height = 1;
    for (const auto& s : demo.trajectory)
        maximum_height = std::max(maximum_height, s.altitude_km);
    const double low = altitude ? 0 : -1.1;
    const double high = altitude ? maximum_height * 1.08 : 1.1;
    for (int index : {vtkAxis::LEFT, vtkAxis::BOTTOM}) {
        auto* axis = chart->GetAxis(index);
        axis->GetLabelProperties()->SetColor(0.70, 0.79, 0.86);
        axis->GetLabelProperties()->SetFontSize(13);
        axis->GetTitleProperties()->SetColor(0.77, 0.85, 0.92);
        axis->GetTitleProperties()->SetFontSize(14);
        axis->GetGridPen()->SetColor(29, 45, 59, 255);
        axis->GetPen()->SetColor(85, 109, 128, 255);
        axis->SetBehavior(vtkAxis::FIXED);
        axis->SetNumberOfTicks(5);
        axis->SetPrecision(3);
    }
    chart->GetAxis(vtkAxis::BOTTOM)->SetTitle("Physical time (ms)");
    chart->GetAxis(vtkAxis::BOTTOM)->SetRange(0, std::max(duration_ms, 1e-9));
    chart->GetAxis(vtkAxis::LEFT)->SetTitle(altitude ? "Altitude (km)" : "p_parallel / |p|");
    chart->GetAxis(vtkAxis::LEFT)->SetRange(low, high);
    vtkNew<vtkTable> data;
    vtkNew<vtkDoubleArray> times;
    vtkNew<vtkDoubleArray> quantities;
    times->SetName("time_ms");
    quantities->SetName("value");
    for (const auto& s : demo.trajectory) {
        times->InsertNextValue(s.time_s * 1000);
        quantities->InsertNextValue(altitude ? s.altitude_km : s.parallel_fraction);
    }
    data->AddColumn(times);
    data->AddColumn(quantities);
    auto* curve = chart->AddPlot(vtkChart::LINE);
    curve->SetInputData(data, 0, 1);
    if (altitude)
        curve->SetColor(244, 192, 99);
    else
        curve->SetColor(55, 219, 211);
    curve->SetWidth(2.0);
    auto reference =
        two_point_table(0, duration_ms, altitude ? demo.parameters.atmosphere_height_m / 1000 : 0,
            altitude ? demo.parameters.atmosphere_height_m / 1000 : 0);
    auto* reference_plot = chart->AddPlot(vtkChart::LINE);
    reference_plot->SetInputData(reference, 0, 1);
    reference_plot->SetColor(115, 145, 161, 180);
    reference_plot->SetWidth(1);
    for (const auto& event : demo.events) {
        auto marker =
            two_point_table(event.state.time_s * 1000, event.state.time_s * 1000, low, high);
        auto* plot = chart->AddPlot(vtkChart::LINE);
        plot->SetInputData(marker, 0, 1);
        plot->SetColor(247, 142, 100, 150);
        plot->SetWidth(1);
    }
    auto cursor = two_point_table(0, 0, low, high);
    auto* cursor_plot = chart->AddPlot(vtkChart::LINE);
    cursor_plot->SetInputData(cursor, 0, 1);
    cursor_plot->SetColor(236, 247, 255, 225);
    cursor_plot->SetWidth(1.5);
    return {chart, cursor};
}

class Viewer {
  public:
    explicit Viewer(const std::vector<DemoCase>& cases, bool offscreen = false) : cases_(cases) {
        if (cases.empty())
            throw std::invalid_argument("No visualization cases");
        window_ = vtkSmartPointer<vtkRenderWindow>::New();
        window_->SetSize(1500, 960);
        window_->SetWindowName("Aurora | Electron dynamics microproject");
        window_->SetMultiSamples(0);
        window_->SetAlphaBitPlanes(1);
        if (offscreen) {
            window_->SetOffScreenRendering(1);
            window_->SetShowWindow(false);
        }
        global_ = renderer(0, 0, 0.64, 1);
        altitude_ = renderer(0.64, 0.67, 1, 1);
        momentum_ = renderer(0.64, 0.34, 1, 0.67);
        local_ = renderer(0.64, 0, 1, 0.34);
        altitude_->InteractiveOff();
        momentum_->InteractiveOff();
        interactor_ = vtkSmartPointer<vtkRenderWindowInteractor>::New();
        interactor_->SetRenderWindow(window_);
        vtkNew<vtkInteractorStyleTrackballCamera> style;
        interactor_->SetInteractorStyle(style);
        lookup_ = colors();
        select(0);
    }

    void select(std::size_t index) {
        if (index >= cases_.size())
            return;
        selected_ = index;
        time_ = 0;
        playing_ = true;
        last_tick_ = Clock::now();
        global_->RemoveAllViewProps();
        local_->RemoveAllViewProps();
        altitude_->RemoveAllViewProps();
        momentum_->RemoveAllViewProps();
        const auto& demo = current();
        const auto prepared = prepare_experiment(demo.parameters);
        global_->AddActor(sphere(1, 0.055, 0.18, 0.31));
        global_->AddActor(sphere(1.0004, 0.16, 0.31, 0.43, 0.5, true));
        const double atmospheric_radius = 1 + demo.parameters.atmosphere_height_m / R;
        global_->AddActor(sphere(atmospheric_radius, 0.12, 0.85, 0.69, 0.12));
        global_->AddActor(sphere(atmospheric_radius, 0.10, 0.54, 0.47, 0.18, true));
        const Vec3 radial = normalized(prepared.initial.position_m);
        const Vec3 axis = prepared.field.axis();
        const Vec3 sideways = cross(axis, radial);
        for (double scale : {0.75, 1.0, 1.25})
            for (int azimuth = 0; azimuth < 4; ++azimuth) {
                const double req = norm(prepared.initial.position_m) * scale;
                if (req <= atmospheric_radius * R)
                    continue;
                const double phi = azimuth * std::numbers::pi / 2;
                auto line = dipole_field_line(prepared.field,
                    (radial * std::cos(phi) + sideways * std::sin(phi)) * req,
                    atmospheric_radius * R);
                for (auto& p : line)
                    p = p / R;
                auto actor = line_actor(line, {}, 0, nullptr);
                actor->GetProperty()->SetColor(0.24, 0.43, 0.58);
                actor->GetProperty()->SetOpacity(0.35);
                global_->AddActor(actor);
            }
        std::vector<Vec3> path;
        std::vector<double> values;
        for (const auto& s : demo.trajectory) {
            path.push_back(s.position_m / R);
            values.push_back(s.parallel_fraction);
        }
        global_->AddActor(line_actor(path, values, 0.0045, lookup_));
        auto start = sphere(0.018, 0.55, 0.72, 1.0);
        const auto p0 = demo.trajectory.front().position_m / R;
        start->SetPosition(p0.x, p0.y, p0.z);
        global_->AddActor(start);
        for (const auto& event : demo.events) {
            const bool mirror = event.type == EventType::Mirror;
            auto marker = sphere(0.024, 1.0, mirror ? 0.76 : 0.29, mirror ? 0.25 : 0.27);
            const auto p = event.state.position_m / R;
            marker->SetPosition(p.x, p.y, p.z);
            global_->AddActor(marker);
            vtkNew<vtkBillboardTextActor3D> text;
            const char* event_label = mirror
                ? "Mirror"
                : (event.type == EventType::InnerBoundary ? "Atmosphere entry" : "Outer boundary");
            text->SetInput(event_label);
            text->SetPosition(p.x + 0.045, p.y - 0.025, p.z + 0.035);
            text->GetTextProperty()->SetFontSize(16);
            text->GetTextProperty()->SetColor(1.0, 0.81, 0.53);
            global_->AddActor(text);
        }
        particle_ = sphere(0.014, 0.98, 0.99, 1.0);
        particle_->GetProperty()->LightingOff();
        global_->AddActor(particle_);
        vtkNew<vtkScalarBarActor> bar;
        bar->SetLookupTable(lookup_);
        bar->SetTitle("p_parallel / |p|");
        bar->SetOrientationToHorizontal();
        bar->SetNumberOfLabels(3);
        bar->SetPosition(0.045, 0.068);
        bar->SetWidth(0.42);
        bar->SetHeight(0.065);
        bar->GetTitleTextProperty()->SetFontSize(14);
        bar->GetLabelTextProperty()->SetFontSize(13);
        global_->AddViewProp(bar);
        std::ostringstream title;
        title << "AURORA / electron dynamics\n"
              << selected_ + 1 << "  |  " << demo.parameters.energy_keV << " keV  |  pitch "
              << demo.parameters.pitch_deg << " deg  |  loss cone " << std::fixed
              << std::setprecision(2)
              << prepared.loss_cone.critical_angle_rad * 180 / std::numbers::pi << " deg";
        label(global_, title.str(), 0.035, 0.975, 22);
        readout_ = label(global_, "", 0.035, 0.885, 15);
        label(global_,
            "Coordinates: Earth radii | magnetic field lines: blue\nLine width and particle/event markers are symbolic",
            0.035, 0.185, 13);
        label(global_,
            "1 / 2 / 3  experiments     SPACE  play/pause     R  restart\nLeft / Right  seek     +/-  speed     S  save PNG     Q  close",
            0.035, 0.045, 13);
        height_chart_ = make_chart(altitude_, demo, true);
        momentum_chart_ = make_chart(momentum_, demo, false);
        build_local(prepared);
        const double hemisphere =
            !demo.events.empty() && demo.events.front().state.position_m.z < 0 ? -1 : 1;
        auto* camera = global_->GetActiveCamera();
        camera->SetPosition(5, -8, 3.5 * hemisphere);
        camera->SetFocalPoint(0.35, 0, 0.25 * hemisphere);
        camera->SetViewUp(0, 0, 1);
        global_->ResetCamera();
        camera->Zoom(1.12);
        global_->ResetCameraClippingRange();
        update();
    }

    void set_fraction(double fraction) {
        playing_ = false;
        time_ = std::clamp(fraction, 0.0, 1.0) * current().trajectory.back().time_s;
        update();
    }

    void capture(const std::filesystem::path& path) {
        if (path.has_parent_path())
            std::filesystem::create_directories(path.parent_path());
        window_->Render();
        vtkNew<vtkWindowToImageFilter> image;
        image->SetInput(window_);
        image->SetInputBufferTypeToRGB();
        image->ReadFrontBufferOff();
        image->Update();
        vtkNew<vtkPNGWriter> writer;
        writer->SetFileName(path.string().c_str());
        writer->SetInputConnection(image->GetOutputPort());
        writer->Write();
        if (writer->GetErrorCode() != vtkErrorCode::NoError || !std::filesystem::exists(path) ||
            std::filesystem::file_size(path) == 0)
            throw std::runtime_error("Cannot save VTK screenshot");
    }

    void interact() {
        vtkNew<vtkCallbackCommand> callback;
        callback->SetClientData(this);
        callback->SetCallback([](vtkObject* caller, unsigned long event, void* data, void*) {
            auto* viewer = static_cast<Viewer*>(data);
            try {
                if (event == vtkCommand::TimerEvent)
                    viewer->tick();
                else
                    viewer->key(vtkRenderWindowInteractor::SafeDownCast(caller)->GetKeySym());
            } catch (const std::exception& error) {
                std::cerr << "VTK: " << error.what() << '\n';
            }
        });
        interactor_->AddObserver(vtkCommand::TimerEvent, callback);
        interactor_->AddObserver(vtkCommand::KeyPressEvent, callback);

        // Keep VTK's built-in character shortcuts from handling our playback keys again.
        vtkNew<vtkCallbackCommand> characters;
        characters->SetClientData(characters);
        characters->SetCallback([](vtkObject* caller, unsigned long, void* data, void*) {
            auto* command = static_cast<vtkCallbackCommand*>(data);
            auto* input = vtkRenderWindowInteractor::SafeDownCast(caller);
            command->SetAbortFlag(is_control_key(input->GetKeySym()) ? 1 : 0);
        });
        interactor_->AddObserver(vtkCommand::CharEvent, characters, 1.0f);

        window_->Render();
        interactor_->Initialize();
        interactor_->CreateRepeatingTimer(33);
        last_tick_ = Clock::now();
        interactor_->Start();
    }

  private:
    vtkSmartPointer<vtkRenderer> renderer(double x0, double y0, double x1, double y1) {
        auto r = vtkSmartPointer<vtkRenderer>::New();
        r->SetViewport(x0, y0, x1, y1);
        r->SetBackground(0.025, 0.045, 0.07);
        window_->AddRenderer(r);
        return r;
    }

    const DemoCase& current() const {
        return cases_.at(selected_);
    }

    void build_local(const PreparedExperiment& e) {
        const auto& samples = current().gyro_segment;
        if (samples.size() < 2) {
            label(local_, "Dense gyro segment was not saved", 0.05, 0.93);
            return;
        }
        const Vec3 origin = e.initial.position_m;
        const Vec3 b = normalized(e.field.at(origin, 0));
        const Vec3 radial = normalized(origin);
        const Vec3 side = cross(b, radial);
        std::vector<Vec3> positions;
        std::vector<double> values;
        double span = 1;
        for (const auto& s : samples) {
            const Vec3 delta = s.position_m - origin;
            const Vec3 local{dot(delta, radial), dot(delta, side), dot(delta, b)};
            positions.push_back(local);
            values.push_back(s.parallel_fraction);
            span = std::max(span, norm(local));
        }
        const double gyroradius = larmor_radius_m(e.initial.momentum_kg_m_s, e.field.at(origin, 0));
        const double line_radius = std::max({0.3, gyroradius * 0.035, span * 0.002});
        local_->AddActor(line_actor(positions, values, line_radius, lookup_));
        const double axis_length = span * 0.12;
        const Vec3 ends[] = {{axis_length, 0, 0}, {0, axis_length, 0}, {0, 0, axis_length}};
        const char* names[] = {"e1", "e2", "B"};
        for (std::size_t i = 0; i < 3; ++i) {
            auto axis = line_actor({Vec3{}, ends[i]}, {}, 0, nullptr);
            axis->GetProperty()->SetColor(
                i == 0 ? 0.9 : 0.3, i == 1 ? 0.9 : 0.3, i == 2 ? 0.9 : 0.3);
            local_->AddActor(axis);

            vtkNew<vtkBillboardTextActor3D> caption;
            caption->SetInput(names[i]);
            caption->SetPosition(ends[i].x, ends[i].y, ends[i].z);
            caption->GetTextProperty()->SetFontSize(12);
            local_->AddActor(caption);
        }
        std::ostringstream text;
        text << "Gyromotion / initial " << std::fixed << std::setprecision(1)
             << samples.back().time_s * 1e6 << " us\nr_L = " << gyroradius
             << " m | all axes in meters";
        label(local_, text.str(), 0.045, 0.96, 14);
        local_->GetActiveCamera()->SetPosition(span * 3, -span * 5, span * 2);
        local_->GetActiveCamera()->SetViewUp(0, 0, 1);
        local_->ResetCamera();
        local_->GetActiveCamera()->Zoom(0.78);
        local_->ResetCameraClippingRange();
    }

    void update() {
        const auto& demo = current();
        const auto sample = sample_at_time(demo.trajectory, time_);
        const Vec3 p = sample.position_m / R;
        particle_->SetPosition(p.x, p.y, p.z);
        for (auto* chart : {&height_chart_, &momentum_chart_}) {
            auto* x = vtkDoubleArray::SafeDownCast(chart->cursor->GetColumn(0));
            x->SetValue(0, time_ * 1000);
            x->SetValue(1, time_ * 1000);
            x->Modified();
            chart->cursor->Modified();
        }
        std::ostringstream text;
        text << std::fixed << std::setprecision(3) << "t = " << time_ * 1000 << " / "
             << demo.trajectory.back().time_s * 1000 << " ms"
             << "    " << (playing_ ? "PLAY" : "PAUSED") << "  x" << speed_
             << "\nAltitude = " << sample.altitude_km
             << " km    p_parallel/|p| = " << sample.parallel_fraction
             << "\nRecorded events: " << demo.events.size()
             << "    Stop: " << outcome_name(demo.outcome) << "\n"
             << std::scientific << std::setprecision(2) << "Energy drift = " << demo.energy_error
             << "    mu variation = " << demo.mu_error;
        readout_->SetInput(text.str().c_str());
    }

    void tick() {
        const auto now = Clock::now();
        const double elapsed = std::chrono::duration<double>(now - last_tick_).count();
        last_tick_ = now;
        if (!playing_)
            return;
        const double duration = current().trajectory.back().time_s;
        time_ = std::min(duration, time_ + elapsed * duration / 12.0 * speed_);
        if (time_ == duration)
            playing_ = false;
        update();
        window_->Render();
    }

    void key(const std::string& key) {
        if (key.size() == 1 && key[0] >= '1' && key[0] <= '3')
            select(static_cast<std::size_t>(key[0] - '1'));
        else if (key == "space") {
            playing_ = !playing_;
            last_tick_ = Clock::now();
        } else if (key == "r" || key == "R") {
            time_ = 0;
            playing_ = true;
            last_tick_ = Clock::now();
        } else if (key == "Left" || key == "Right") {
            playing_ = false;
            time_ = std::clamp(
                time_ + (key == "Left" ? -0.02 : 0.02) * current().trajectory.back().time_s, 0.0,
                current().trajectory.back().time_s);
        } else if (key == "plus" || key == "equal")
            speed_ = std::min(8.0, speed_ * 2);
        else if (key == "minus")
            speed_ = std::max(0.125, speed_ / 2);
        else if (key == "s" || key == "S") {
            const auto stamp = std::chrono::duration_cast<std::chrono::milliseconds>(
                Clock::now().time_since_epoch())
                                   .count();
            const auto path = current().directory / ("scene-" + std::to_string(stamp) + ".png");
            capture(path);
            std::cout << "PNG: " << path.string() << '\n';
        } else if (key == "q" || key == "Q") {
            interactor_->TerminateApp();
            return;
        }
        update();
        window_->Render();
    }

    const std::vector<DemoCase>& cases_;
    std::size_t selected_{};
    double time_{}, speed_{1};
    bool playing_{true};
    Clock::time_point last_tick_{Clock::now()};
    vtkSmartPointer<vtkRenderWindow> window_;
    vtkSmartPointer<vtkRenderWindowInteractor> interactor_;
    vtkSmartPointer<vtkRenderer> global_, altitude_, momentum_, local_;
    vtkSmartPointer<vtkLookupTable> lookup_;
    vtkSmartPointer<vtkActor> particle_;
    vtkSmartPointer<vtkTextActor> readout_;
    Chart height_chart_, momentum_chart_;
};
}

void show_demo(const std::vector<DemoCase>& cases) {
    Viewer(cases).interact();
}

void render_demo_snapshot(const std::vector<DemoCase>& cases, std::size_t selected_case,
    double time_fraction, const std::filesystem::path& path) {
    Viewer viewer(cases, true);
    viewer.select(selected_case);
    viewer.set_fraction(time_fraction);
    viewer.capture(path);
}
} // namespace aurora::visualization
