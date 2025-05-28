cmake_minimum_required(VERSION 3.1.0)
project(shm3d_bindings)

add_subdirectory(deps/signed-heat-3d)
add_subdirectory(deps/nanobind)

# Try to import all Python components potentially needed by nanobind
find_package(Python 3.8
  REQUIRED COMPONENTS Interpreter Development.Module
  OPTIONAL_COMPONENTS Development.SABIModule)

nanobind_add_module(shm3d_bindings
  src/cpp/solvers.cpp
)

include_directories(shm3d_bindings ${CMAKE_CURRENT_SOURCE_DIR}/src/cpp)
include_directories(shm3d_bindings ${CMAKE_CURRENT_SOURCE_DIR}/deps/signed-heat-3d/include)

# Include directories from the submodules of signed-heat-3d... must be a better way?
include_directories(shm3d_bindings ${CMAKE_CURRENT_SOURCE_DIR}/deps/signed-heat-3d/deps/geometry-central/include)
include_directories(shm3d_bindings ${CMAKE_CURRENT_SOURCE_DIR}/deps/signed-heat-3d/deps/polyscope/include)
include_directories(shm3d_bindings ${CMAKE_CURRENT_SOURCE_DIR}/deps/signed-heat-3d/deps/libigl/include)
include_directories(shm3d_bindings ${CMAKE_CURRENT_SOURCE_DIR}/deps/signed-heat-3d/deps/tetgen)

# Include targets the submodules of signed-heat-3d... must be a better way?
target_link_libraries(shm3d_bindings PRIVATE signed-heat-3d geometry-central polyscope tetgen)

install(TARGETS shm3d_bindings LIBRARY DESTINATION .)
