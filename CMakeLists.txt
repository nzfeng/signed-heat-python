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

add_definitions(-UAMGCL_NO_BOOST) # explicitly un-define AMGCL's "ANGCL_NO_BOOST" variable, so we do indeed use Boost

# Download Boost at CMake-time, if the user doesn't have it already.

#set(BOOST_INCLUDE_LIBRARIES property_tree)
#set(BOOST_ENABLE_CMAKE ON)
#include(FetchContent)
#FetchContent_Declare(
#  Boost
#  GIT_REPOSITORY https://github.com/boostorg/boost.git
#  GIT_TAG boost-1.80.0
#  GIT_SHALLOW TRUE
#)
#FetchContent_MakeAvailable(Boost)


if(POLICY CMP0167)
  cmake_policy(SET CMP0167 NEW)
endif()
find_package(Boost REQUIRED)
target_include_directories(shm3d_bindings PRIVATE ${Boost_INCLUDE_DIRS})

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
