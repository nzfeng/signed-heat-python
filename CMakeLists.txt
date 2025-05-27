cmake_minimum_required(VERSION 3.1.0)
project(shm3d_bindings)

add_subdirectory(deps/signed-heat-3d)
add_subdirectory(deps/nanobind)

nanobind_add_module(shm3d_bindings
  src/cpp/solvers.cpp
)

include_directories(shm3d_bindings ${CMAKE_CURRENT_SOURCE_DIR}/src/cpp)
#include_directories(shm3d_bindings ${CMAKE_CURRENT_SOURCE_DIR}/deps/signed-heat-3d/include)

target_link_libraries(shm3d_bindings PRIVATE signed-heat-3d)

install(TARGETS shm3d_bindings LIBRARY DESTINATION .)
