#include "aurora/diagnostics/UniformFieldSolution.hpp"
#include "aurora/io/TrajectoryCsv.hpp"
#include "aurora/physics/Constants.hpp"
#include "aurora/physics/Relativity.hpp"
#include "aurora/simulation/TrajectoryRunner.hpp"
#include <algorithm>
#include <functional>
#include <iostream>
#include <limits>
#include <map>
#include <numbers>
#include <sstream>
#include <string>
#include <vector>

namespace {
using namespace aurora;
void require(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}
void near(double actual, double expected, double tolerance, const std::string& label) {
    if (!std::isfinite(actual) || std::abs(actual-expected)>tolerance) {
        std::ostringstream error;
        error.precision(17);
        error << label << ": actual=" << actual << " expected=" << expected << " tolerance=" << tolerance;
        throw std::runtime_error(error.str());
    }
}
template<class F> void rejects(F&& action, const std::string& label) {
    try { action(); }
    catch (const std::exception&) { return; }
    throw std::runtime_error("Expected rejection: " + label);
}
void vectors() {
    const Vec3 a{1,2,3}, b{-2,4,1};
    near(dot(a,b),9,0,"dot");
    const Vec3 c=cross(a,b);
    near(norm(c-Vec3{-10,-7,8}),0,0,"cross");
    near(dot(c,a),0,0,"orthogonality a");
    near(dot(c,b),0,0,"orthogonality b");
    near(norm(normalized(a)),1,1e-15,"unit vector");
    near(norm(Vec3{3e200,4e200,0})/1e200,5,1e-15,"scaled norm");
    rejects([]{ (void)normalized({}); },"zero normalization");
}
void relativity() {
    for (double eV : {1e-12,1e-6,1.0,1e4,511e3,1e8,1e12}) {
        const double energy=constants::eV_to_J(eV);
        const Vec3 p=momentum_from_energy(energy,{1,2,-3});
        near(kinetic_energy_J(p)/energy,1,2e-15,"energy round trip");
        near(lorentz_factor(p),1+energy/constants::electron_rest_energy_J,
             2e-15*(1+energy/constants::electron_rest_energy_J),"gamma from energy");
        const auto v=velocity_from_momentum(p);
        require(norm(v)<constants::speed_of_light_m_s,"subluminal speed");
        if (eV<=1e8) near(norm(momentum_from_velocity(v)-p)/norm(p),0,1e-10,"velocity round trip");
    }
    const auto slow=momentum_from_velocity({100,0,0});
    near(kinetic_energy_J(slow)/(0.5*constants::electron_mass_kg*10000),1,1e-12,"nonrelativistic limit");
    near(kinetic_energy_J({}),0,0,"rest energy is not kinetic energy");
    near(norm(momentum_from_energy(0,{1,0,0})),0,0,"stationary electron");
}
void zero_field() {
    const ZeroField field;
    const ParticleState initial{{1,-3,5},momentum_from_energy(constants::eV_to_J(4e6),{1,2,3}),2e-6};
    const auto summary=run_trajectory(initial,field,{1e-9,1000,1});
    const Vec3 expected=initial.position_m+velocity_from_momentum(initial.momentum_kg_m_s)*1e-6;
    near(norm(summary.final_state.position_m-expected),0,1e-10,"straight trajectory");
    near(norm(summary.final_state.momentum_kg_m_s-initial.momentum_kg_m_s),0,0,"constant momentum");
    near(summary.final_state.time_s,3e-6,1e-21,"time");
    near(summary.max_relative_energy_error,0,0,"zero-field energy");
    const auto reference=uniform_field_solution(initial,{},3e-6);
    near(norm(reference.position_m-expected),0,1e-12,"zero-field analytic limit");

    std::vector<std::size_t> seen;
    run_trajectory(initial,field,{1e-9,10,3},[&](std::size_t step,const ParticleState&){seen.push_back(step);});
    require(seen==std::vector<std::size_t>{0,3,6,9,10},"sampling includes final incomplete stride");
    seen.clear();
    const auto empty=run_trajectory(initial,field,{1e-9,0,3},[&](std::size_t step,const ParticleState&){seen.push_back(step);});
    require(seen==std::vector<std::size_t>{0} && empty.completed_steps==0,"zero-step run");
}
void uniform_field() {
    // Independent reference for B along z: initially px>0, py=0 -> py grows for q_e<0.
    const Vec3 b{0,0,3e-5};
    const ParticleState start{{},momentum_from_energy(constants::eV_to_J(10000),{1,0,0}),0};
    const double period=gyroperiod_s(start.momentum_kg_m_s,b);
    const double radius=larmor_radius_m(start.momentum_kg_m_s,b);
    const auto quarter=uniform_field_solution(start,b,period/4);
    near(quarter.position_m.x/radius,1,1e-14,"analytic quarter-turn x");
    near(quarter.position_m.y/radius,1,1e-14,"analytic quarter-turn y");
    require(quarter.momentum_kg_m_s.y>0,"analytic electron rotation sign");
    const auto numerical=run_trajectory(start,UniformField(b),{period/1000,250,1});
    near(numerical.final_state.position_m.x/radius,1,2e-5,"numerical quarter-turn x");
    near(numerical.final_state.position_m.y/radius,1,2e-5,"numerical quarter-turn y");
    require(numerical.final_state.momentum_kg_m_s.y>0,"numerical electron rotation sign");

    for (Vec3 field_vector : {b,-b,Vec3{1e-5,2e-5,-3e-5}}) {
        const UniformField field(field_vector);
        for (double energy_eV : {1.0,1e4,1e7}) {
            const ParticleState initial{{2,-1,3},momentum_from_energy(constants::eV_to_J(energy_eV),{2,-1,4}),0};
            const double dt=gyroperiod_s(initial.momentum_kg_m_s,field_vector)/1000;
            double max_r_error=0, max_p_error=0;
            const auto result=run_trajectory(initial,field,{dt,1000,1},[&](std::size_t,const ParticleState& s) {
                const auto exact=uniform_field_solution(initial,field_vector,s.time_s);
                max_r_error=std::max(max_r_error,norm(s.position_m-exact.position_m));
                max_p_error=std::max(max_p_error,norm(s.momentum_kg_m_s-exact.momentum_kg_m_s));
            });
            near(max_r_error/larmor_radius_m(initial.momentum_kg_m_s,field_vector),0,3e-5,"helical position");
            near(max_p_error/norm(initial.momentum_kg_m_s),0,3e-5,"helical momentum");
            near(result.max_relative_energy_error,0,1e-11,"helical energy");
        }
    }
    const ParticleState parallel{{},momentum_from_energy(constants::eV_to_J(1e5),{0,0,-1}),0};
    const auto along=run_trajectory(parallel,UniformField(b),{period/100,100,1});
    near(norm(along.final_state.position_m-velocity_from_momentum(parallel.momentum_kg_m_s)*period),0,1e-10,"parallel motion");
    const auto rest=run_trajectory({{1,2,3},{},0},UniformField(b),{period/100,100,1});
    near(norm(rest.final_state.position_m-Vec3{1,2,3}),0,0,"rest in B");
    near(rest.max_relative_energy_error,0,0,"stationary energy diagnostic");

    const RelativisticBorisPusher pusher(period/200);
    const UniformField field(b);
    const auto synchronized=pusher.synchronize(pusher.initialize(start,field),field);
    near(norm(synchronized.momentum_kg_m_s-start.momentum_kg_m_s)/norm(start.momentum_kg_m_s),0,1e-15,"half-step initialization inverse");
}
void convergence() {
    const Vec3 b{0,0,3e-5};
    const UniformField field(b);
    const ParticleState initial{{},momentum_from_energy(constants::eV_to_J(10000),{std::sqrt(3.0),0,1}),0};
    const double end=0.7*gyroperiod_s(initial.momentum_kg_m_s,b);
    const auto reference=uniform_field_solution(initial,b,end);
    double previous_r=0,previous_p=0;
    for (std::size_t steps : {40u,80u,160u,320u}) {
        const auto result=run_trajectory(initial,field,{end/static_cast<double>(steps),steps,1});
        const double r_error=norm(result.final_state.position_m-reference.position_m);
        const double p_error=norm(result.final_state.momentum_kg_m_s-reference.momentum_kg_m_s);
        if (previous_r>0) {
            near(std::log2(previous_r/r_error),2,0.1,"second-order position");
            near(std::log2(previous_p/p_error),2,0.1,"second-order momentum");
        }
        previous_r=r_error; previous_p=p_error;
    }
}
void long_run() {
    const Vec3 b{0,0,3e-5};
    const UniformField field(b);
    const ParticleState initial{{},momentum_from_energy(constants::eV_to_J(10000),{std::sqrt(3.0),0,1}),0};
    const double period=gyroperiod_s(initial.momentum_kg_m_s,b);
    const auto result=run_trajectory(initial,field,{period/200,20000,137});
    near(result.max_relative_energy_error,0,1e-9,"100-period energy bound");
    require(result.max_speed_m_s<constants::speed_of_light_m_s,"100-period subluminal speed");
    const auto reference=uniform_field_solution(initial,b,result.final_state.time_s);
    const double relative_p_error=norm(result.final_state.momentum_kg_m_s-reference.momentum_kg_m_s)/norm(initial.momentum_kg_m_s);
    require(relative_p_error>0.01 && relative_p_error<0.06,"finite phase error despite conserved energy");
}
void invalid_input() {
    const double nan=std::numeric_limits<double>::quiet_NaN();
    const double inf=std::numeric_limits<double>::infinity();
    for (double dt : {0.0,-1.0,nan,inf}) rejects([&]{RelativisticBorisPusher p(dt);},"invalid dt");
    rejects([&]{UniformField f({inf,0,0});},"infinite B");
    rejects([&]{(void)momentum_from_energy(-1,{1,0,0});},"negative energy");
    rejects([&]{(void)momentum_from_energy(inf,{1,0,0});},"infinite energy");
    rejects([&]{(void)momentum_from_energy(1,{});},"undefined direction");
    rejects([&]{(void)momentum_from_velocity({constants::speed_of_light_m_s,0,0});},"speed c");
    rejects([&]{(void)momentum_from_velocity({2*constants::speed_of_light_m_s,0,0});},"superluminal speed");
    rejects([&]{(void)velocity_from_momentum({nan,0,0});},"nan momentum");
    rejects([&]{(void)gyroperiod_s({1,0,0},{});},"period at B=0");
    rejects([&]{(void)larmor_radius_m({1,0,0},{});},"radius at B=0");
    const ZeroField zero;
    rejects([&]{run_trajectory({{nan,0,0},{},0},zero,{1e-9,1,1});},"nan position");
    rejects([&]{run_trajectory({{}, {},nan},zero,{1e-9,1,1});},"nan time");
    rejects([&]{run_trajectory({{}, {},0},zero,{1e-9,1,0});},"zero stride");
    rejects([&]{run_trajectory({{}, {},0},zero,{1e308,10,1});},"overflowing duration");
    rejects([&]{run_trajectory({{}, {},1e30},zero,{1e-9,1,1});},"unresolvable time increment");
    class BadField final : public MagneticField {
        Vec3 at(Vec3,double) const override { return {std::numeric_limits<double>::quiet_NaN(),0,0}; }
    } bad;
    rejects([&]{run_trajectory({{}, {},0},bad,{1e-9,1,1});},"bad field result");
    RelativisticBorisPusher first(1e-9),second(2e-9);
    auto state=first.initialize({{}, {},0},zero);
    rejects([&]{second.advance(state,zero);},"staggered step change");
}
void csv() {
    std::ostringstream out;
    write_trajectory_header(out);
    const ParticleState initial{{1,2,3},momentum_from_energy(constants::eV_to_J(10000),{1,0,0}),0};
    write_trajectory_row(out,0,initial);
    std::istringstream input(out.str());
    std::string header,line;
    std::getline(input,header); std::getline(input,line);
    require(header.starts_with("step,time_s,x_m,y_m,z_m"),"CSV header and SI units");
    require(std::count(header.begin(),header.end(),',')==12,"CSV column count");
    std::vector<double> values;
    std::istringstream row(line);
    std::string token;
    while (std::getline(row,token,',')) values.push_back(std::stod(token));
    require(values.size()==13,"CSV row column count");
    near(values[5],initial.momentum_kg_m_s.x,0,"CSV double round trip");
    near(values[12],kinetic_energy_J(initial.momentum_kg_m_s),0,"CSV kinetic energy");
    std::ostringstream failed;
    failed.setstate(std::ios::badbit);
    rejects([&]{write_trajectory_header(failed);},"CSV write failure");
}
} // namespace

int main(int argc,char** argv) {
    const std::map<std::string,std::function<void()>> suites{
        {"vectors",vectors},{"relativity",relativity},{"zero_field",zero_field},
        {"uniform_field",uniform_field},{"convergence",convergence},{"long_run",long_run},
        {"invalid_input",invalid_input},{"csv",csv}};
    try {
        if (argc!=2 || !suites.contains(argv[1])) throw std::invalid_argument("Supply a known test suite name");
        suites.at(argv[1])();
        std::cout << "PASS: " << argv[1] << '\n';
        return 0;
    } catch (const std::exception& e) { std::cerr << "FAIL: " << e.what() << '\n'; return 1; }
}
