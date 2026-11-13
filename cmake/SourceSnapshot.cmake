# Keep an immutable source snapshot for each compiled version, including dirty working trees.
file(GLOB_RECURSE AURORA_SOURCES CONFIGURE_DEPENDS RELATIVE "${PROJECT_SOURCE_DIR}"
    "${PROJECT_SOURCE_DIR}/include/*.hpp" "${PROJECT_SOURCE_DIR}/src/*.cpp"
    "${PROJECT_SOURCE_DIR}/tests/*.cpp" "${PROJECT_SOURCE_DIR}/tests/*.hpp"
    "${PROJECT_SOURCE_DIR}/cmake/*.cmake" "${PROJECT_SOURCE_DIR}/cmake/*.in")
list(APPEND AURORA_SOURCES main.cpp CMakeLists.txt)
list(REMOVE_DUPLICATES AURORA_SOURCES)
list(SORT AURORA_SOURCES)
set(AURORA_SOURCE_RECORDS "")
foreach(relative_path IN LISTS AURORA_SOURCES)
    set(absolute_path "${PROJECT_SOURCE_DIR}/${relative_path}")
    set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS "${absolute_path}")
    file(SHA256 "${absolute_path}" file_hash)
    string(APPEND AURORA_SOURCE_RECORDS "${file_hash}  ${relative_path}\n")
endforeach()
string(SHA256 AURORA_SOURCE_SHA256 "${AURORA_SOURCE_RECORDS}")
set(AURORA_SNAPSHOT_DIR "${PROJECT_BINARY_DIR}/source-snapshots/${AURORA_SOURCE_SHA256}")
foreach(relative_path IN LISTS AURORA_SOURCES)
    configure_file("${PROJECT_SOURCE_DIR}/${relative_path}" "${AURORA_SNAPSHOT_DIR}/${relative_path}" COPYONLY)
endforeach()
file(WRITE "${AURORA_SNAPSHOT_DIR}/source.sha256" "${AURORA_SOURCE_RECORDS}")
execute_process(COMMAND git rev-parse HEAD WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}"
    OUTPUT_VARIABLE AURORA_GIT_COMMIT OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
if(NOT AURORA_GIT_COMMIT)
    set(AURORA_GIT_COMMIT unknown)
endif()
execute_process(COMMAND git status --porcelain --untracked-files=normal WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}"
    OUTPUT_VARIABLE git_status OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
if(git_status)
    set(AURORA_GIT_STATE dirty)
else()
    set(AURORA_GIT_STATE clean_or_unavailable)
endif()
string(TOUPPER "${CMAKE_BUILD_TYPE}" build_type_upper)
set(AURORA_RECORDED_FLAGS "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_${build_type_upper}}; sanitizers=${AURORA_ENABLE_SANITIZERS}")
configure_file("${PROJECT_SOURCE_DIR}/cmake/BuildInfo.hpp.in"
    "${PROJECT_BINARY_DIR}/generated/aurora/BuildInfo.hpp" @ONLY)
