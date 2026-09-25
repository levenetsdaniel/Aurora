cmake_minimum_required(VERSION 4.1)
project(Aurora)

set(CMAKE_CXX_STANDARD 20)

add_executable(Aurora main.cpp src/math/Vector3.hpp src/physics/MagneticField.hpp src/physics/DipoleField.hpp)

enable_testing()
add_executable(aurora_boris_test tests/boris_integrator_test.cpp)
add_test(NAME boris_integrator COMMAND aurora_boris_test)
