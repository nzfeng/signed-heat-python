cmake_minimum_required(VERSION 3.5)
set(CMAKE_POLICY_VERSION_MINIMUM 3.5)
project(shm3d_bindings)

# geometry-central
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

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

add_definitions(-UAMGCL_NO_BOOST) # explicitly un-define AMGCL's "ANGCL_NO_BOOST" variable, so we do indeed use Boost

target_link_libraries(shm3d_bindings PRIVATE signed-heat-3d)

install(TARGETS shm3d_bindings LIBRARY DESTINATION .)
