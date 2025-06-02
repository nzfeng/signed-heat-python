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

# Include directories from all the submodules of signed-heat-3d... there's not toooo many of them
include_directories(shm3d_bindings ${CMAKE_CURRENT_SOURCE_DIR}/deps/signed-heat-3d/deps/geometry-central/include)
include_directories(shm3d_bindings ${CMAKE_CURRENT_SOURCE_DIR}/deps/signed-heat-3d/deps/libigl/include)
include_directories(shm3d_bindings ${CMAKE_CURRENT_SOURCE_DIR}/deps/signed-heat-3d/deps/tetgen)
include_directories(shm3d_bindings ${CMAKE_CURRENT_SOURCE_DIR}/deps/signed-heat-3d/deps/amgcl)

add_definitions(-UAMGCL_NO_BOOST)

# Compile C++ code from signed-heat-3d
set(SHM3D_SRCS
  deps/signed-heat-3d/src/signed_heat_3d.cpp
  deps/signed-heat-3d/src/signed_heat_grid_solver.cpp
  deps/signed-heat-3d/src/signed_heat_tet_solver.cpp
)
add_library(signed-heat-3d STATIC ${SHM3D_SRCS})
target_link_libraries(signed-heat-3d PRIVATE geometry-central tetgen amgcl::amgcl)

# Include targets from the submodules of signed-heat-3d, necessary for the cpp files in this project
target_link_libraries(shm3d_bindings PRIVATE signed-heat-3d geometry-central amgcl::amgcl) # confused why amgcl is needed here as well...

install(TARGETS shm3d_bindings LIBRARY DESTINATION .)
