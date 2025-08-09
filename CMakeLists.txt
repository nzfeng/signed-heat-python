cmake_minimum_required(VERSION 3.5)
set(CMAKE_POLICY_VERSION_MINIMUM 3.5)
project(shm3d_bindings)

# geometry-central
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

add_definitions(-UAMGCL_NO_BOOST) # explicitly un-define AMGCL's "ANGCL_NO_BOOST" variable, so we do indeed use Boost

# Make sure we have Boost (for AMGCL)
find_package(Boost)
if(NOT Boost_FOUND)
	message("Boost not found")
	message(STATUS "Downloading and extracting Boost library sources. This will take some time...")
	set(BOOST_INCLUDE_LIBRARIES system filesystem program_options property_tree)
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

	# amgcl only needs property_tree, but in Boost every file is seemingly included by every other file :((
	file(GLOB V_GLOB LIST_DIRECTORIES true ${boost_SOURCE_DIR}/libs/*)
	foreach(item ${V_GLOB})
	  if(IS_DIRECTORY ${item})
	    include_directories(${item}/include/)
	  endif()
	endforeach()

	set(BOOST_ROOT "${boost_SOURCE_DIR}")
	set(BOOST_INCLUDEDIR "${boost_SOURCE_DIR}")
	set(BOOST_LIBRARYDIR "${boost_SOURCE_DIR}")
	set(Boost_FOUND TRUE)
endif()

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

target_link_libraries(shm3d_bindings PRIVATE signed-heat-3d)

install(TARGETS shm3d_bindings LIBRARY DESTINATION .)
