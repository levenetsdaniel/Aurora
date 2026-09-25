cmake_minimum_required(VERSION 4.1)
project(Aurora)

set(CMAKE_CXX_STANDARD 20)

add_executable(Aurora main.cpp src/math/Vector3.hpp src/physics/MagneticField.hpp src/physics/DipoleField.hpp)

enable_testing()
add_executable(aurora_boris_test tests/boris_integrator_test.cpp)
add_test(NAME boris_integrator COMMAND aurora_boris_test)
add_executable(aurora_diagnostics_test tests/diagnostics_test.cpp)
add_test(NAME diagnostics COMMAND aurora_diagnostics_test)
add_executable(aurora_simulation_engine_test tests/simulation_engine_test.cpp)
add_test(NAME simulation_engine COMMAND aurora_simulation_engine_test)
