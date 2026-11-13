cmake_minimum_required(VERSION 3.20)
project(Aurora VERSION 0.1.0 LANGUAGES CXX)
set(CMAKE_CXX_EXTENSIONS OFF)

option(AURORA_ENABLE_SANITIZERS "Enable address and undefined behaviour sanitizers" OFF)
option(BUILD_TESTING "Build physics tests and the validation report tool" OFF)
include(CTest)
include(cmake/SourceSnapshot.cmake)

add_library(aurora_core
    src/fields/DipoleField.cpp
    src/physics/Relativity.cpp
    src/physics/ParticleProperties.cpp
    src/physics/Gyromotion.cpp
    src/physics/LossCone.cpp
    src/numerics/RelativisticBorisPusher.cpp
    src/simulation/TrajectoryRunner.cpp
    src/simulation/AtmosphereExperiment.cpp
    src/simulation/EventLocator.cpp)
target_include_directories(aurora_core PUBLIC include)
target_compile_features(aurora_core PUBLIC cxx_std_20)
set_target_properties(aurora_core PROPERTIES CXX_EXTENSIONS OFF)

add_library(aurora_build_options INTERFACE)
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(aurora_build_options INTERFACE -Wall -Wextra -Wpedantic -Wconversion -Wshadow)
    if(AURORA_ENABLE_SANITIZERS)
        target_compile_options(aurora_build_options INTERFACE -fsanitize=address,undefined -fno-omit-frame-pointer)
        target_link_options(aurora_build_options INTERFACE -fsanitize=address,undefined)
    endif()
elseif(MSVC)
    target_compile_options(aurora_build_options INTERFACE /W4 /permissive-)
endif()
target_link_libraries(aurora_core PUBLIC aurora_build_options)

add_library(aurora_io src/io/TrajectoryCsv.cpp src/io/EventCsv.cpp src/io/ExperimentFiles.cpp)
target_include_directories(aurora_io PRIVATE "${PROJECT_BINARY_DIR}/generated")
target_compile_definitions(aurora_io PRIVATE AURORA_BUILD_TYPE="$<CONFIG>")
target_link_libraries(aurora_io PUBLIC aurora_core)
add_executable(Aurora main.cpp)
target_link_libraries(Aurora PRIVATE aurora_io)

if(BUILD_TESTING)
    add_executable(aurora_week3_report tests/week3_experiments.cpp)
    target_link_libraries(aurora_week3_report PRIVATE aurora_io)
    add_executable(aurora_week3_tests tests/week3_tests.cpp)
    target_link_libraries(aurora_week3_tests PRIVATE aurora_io)
    foreach(suite IN ITEMS classification parameters hemispheres reproducibility)
        add_test(NAME ${suite} COMMAND aurora_week3_tests ${suite})
        set_tests_properties(${suite} PROPERTIES WORKING_DIRECTORY "${PROJECT_BINARY_DIR}")
    endforeach()
    add_executable(aurora_replay tests/replay_experiment.cpp)
    target_link_libraries(aurora_replay PRIVATE aurora_io)
    add_executable(aurora_loss_cone_tests tests/loss_cone_tests.cpp)
    target_link_libraries(aurora_loss_cone_tests PRIVATE aurora_core)
    add_test(NAME loss_cone COMMAND aurora_loss_cone_tests)
    add_executable(aurora_atmosphere_tests tests/atmosphere_tests.cpp)
    target_link_libraries(aurora_atmosphere_tests PRIVATE aurora_io)
    add_test(NAME atmosphere_entry COMMAND aurora_atmosphere_tests)
    # Analytical trajectories are test references, separate from the simulation library.
    add_library(aurora_test_reference tests/reference/UniformFieldSolution.cpp)
    target_include_directories(aurora_test_reference PUBLIC tests)
    target_link_libraries(aurora_test_reference PUBLIC aurora_core)
    add_executable(aurora_week2_report tests/week2_validation.cpp)
    target_link_libraries(aurora_week2_report PRIVATE aurora_core)
    add_executable(aurora_week2_tests tests/week2_tests.cpp)
    target_link_libraries(aurora_week2_tests PRIVATE aurora_io)
    foreach(suite IN ITEMS dipole particle_properties events limits mirrors dipole_convergence)
        add_test(NAME ${suite} COMMAND aurora_week2_tests ${suite})
    endforeach()
    add_executable(aurora_validation_report tests/validation_report.cpp)
    target_link_libraries(aurora_validation_report PRIVATE aurora_test_reference)
    add_executable(aurora_tests tests/week1_tests.cpp)
    target_link_libraries(aurora_tests PRIVATE aurora_io aurora_test_reference)
    foreach(suite IN ITEMS vectors relativity zero_field uniform_field convergence long_run invalid_input csv)
        add_test(NAME ${suite} COMMAND aurora_tests ${suite})
    endforeach()
    file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/test-results/example)
    add_test(NAME example_run COMMAND Aurora)
    set_tests_properties(example_run PROPERTIES
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/test-results/example)
    add_test(NAME validation_report COMMAND aurora_validation_report
        --output ${CMAKE_CURRENT_BINARY_DIR}/test-results/validation.md)
endif()
