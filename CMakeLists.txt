cmake_minimum_required(VERSION 3.5)
set(CMAKE_POLICY_VERSION_MINIMUM 3.5)
project(shm3d_bindings)

# Copied from signed-heat-3d (and hence gc-polyscope-project-template)
if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang" OR "${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")

  # using Clang (linux or apple) or GCC
  message("Using clang/gcc compiler flags")
  SET(CMAKE_CXX_STANDARD 11)
  SET(BASE_CXX_FLAGS "-std=c++11 -Wall -Wextra")
  SET(DISABLED_WARNINGS " -Wno-unused-parameter -Wno-unused-variable -Wno-unused-function -Wno-deprecated-declarations -Wno-missing-braces -Wno-unused-label -Wno-register -Wno-sign-compare -Wno-unknown-pragmas -Wno-unused-result -Wno-narrowing -Wno-unused-but-set-variable -Wno-unused-lambda-capture -Wno-unused-local-typedefs")
  SET(TRACE_INCLUDES " -H -Wno-error=unused-command-line-argument")

  if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
    message("Setting clang-specific options")
    SET(BASE_CXX_FLAGS "${BASE_CXX_FLAGS} -ferror-limit=3 -Wfatal-errors -fcolor-diagnostics")
    SET(CMAKE_CXX_FLAGS_DEBUG          "-g3 -fsanitize=address -fno-limit-debug-info")
  elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    message("Setting gcc-specific options")
    SET(BASE_CXX_FLAGS "${BASE_CXX_FLAGS} -fmax-errors=5")
    SET(CMAKE_CXX_FLAGS_DEBUG          "-g3")
    SET(DISABLED_WARNINGS "${DISABLED_WARNINGS} -Wno-maybe-uninitialized -Wno-format-zero-length -Wno-unused-but-set-parameter -Wno-unused-but-set-variable")
  endif()

  SET(CMAKE_CXX_FLAGS "${BASE_CXX_FLAGS} ${DISABLED_WARNINGS} -std=c++11")

  include(CheckCXXCompilerFlag)
  unset(COMPILER_SUPPORTS_MARCH_NATIVE CACHE)
  CHECK_CXX_COMPILER_FLAG(-march=native COMPILER_SUPPORTS_MARCH_NATIVE)
  if(COMPILER_SUPPORTS_MARCH_NATIVE)
    SET(CMAKE_CXX_FLAGS_RELEASE        "-O3 -march=native -DNDEBUG -std=c++11")
  else()
    SET(CMAKE_CXX_FLAGS_RELEASE        "-O3 -mcpu=apple-m1 -DNDEBUG -std=c++11") # Apple M1
  endif()

elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
  # using Visual Studio C++
  message("Using Visual Studio compiler flags")
  set(BASE_CXX_FLAGS "${BASE_CXX_FLAGS} /W4")
  set(BASE_CXX_FLAGS "${BASE_CXX_FLAGS} /MP") # parallel build
  SET(DISABLED_WARNINGS "${DISABLED_WARNINGS} /wd\"4267\"")  # ignore conversion to smaller type (fires more aggressively than the gcc version, which is annoying)
  SET(DISABLED_WARNINGS "${DISABLED_WARNINGS} /wd\"4244\"")  # ignore conversion to smaller type (fires more aggressively than the gcc version, which is annoying)
  SET(DISABLED_WARNINGS "${DISABLED_WARNINGS} /wd\"4305\"")  # ignore truncation on initialization
  SET(CMAKE_CXX_FLAGS "${BASE_CXX_FLAGS} ${DISABLED_WARNINGS}")
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /MD")
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /MDd")

  add_definitions(/D "_CRT_SECURE_NO_WARNINGS")
  add_definitions(-DNOMINMAX)
  add_definitions(-D_USE_MATH_DEFINES)
else()
  # unrecognized
  message( FATAL_ERROR "Unrecognized compiler [${CMAKE_CXX_COMPILER_ID}]" )
endif()

# set location of eigen for geometry-central
set(GC_EIGEN_LOCATION "${CMAKE_CURRENT_SOURCE_DIR}/deps/eigen" CACHE PATH "my path")

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

# Include directories from all the submodules of signed-heat-3d... there's not toooo many of them
include_directories(shm3d_bindings ${CMAKE_CURRENT_SOURCE_DIR}/deps/signed-heat-3d/deps/geometry-central/include)
include_directories(shm3d_bindings ${CMAKE_CURRENT_SOURCE_DIR}/deps/signed-heat-3d/deps/libigl/include)
include_directories(shm3d_bindings ${CMAKE_CURRENT_SOURCE_DIR}/deps/signed-heat-3d/deps/tetgen)
include_directories(shm3d_bindings ${CMAKE_CURRENT_SOURCE_DIR}/deps/signed-heat-3d/deps/amgcl)

add_definitions(-UAMGCL_NO_BOOST) # explicitly un-define AMGCL's "ANGCL_NO_BOOST" variable, so we do indeed use Boost

# Make sure we have Boost (for AMGCL)
if(POLICY CMP0167)
  cmake_policy(SET CMP0167 NEW)
endif()
find_package(Boost)
if(NOT Boost_FOUND)
	message("Boost not found")
	message(STATUS "Downloading and extracting Boost library sources. This will take some time...")
	set(BOOST_INCLUDE_LIBRARIES property_tree)
	set(BOOST_ENABLE_CMAKE ON)
	Set(FETCHCONTENT_QUIET FALSE) # print downloading progress
	include(FetchContent)
	FetchContent_Declare(
	  Boost
	  GIT_REPOSITORY https://github.com/boostorg/boost.git
	  GIT_TAG boost-1.73.0
	  GIT_SHALLOW TRUE
	)
	FetchContent_MakeAvailable(Boost)
	find_package(Boost)
	message(" ${Boost_LIBRARIES}")
	message(" ${Boost_INCLUDE_DIRS}")
	include_directories(${Boost_INCLUDE_DIRS})
endif()

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
