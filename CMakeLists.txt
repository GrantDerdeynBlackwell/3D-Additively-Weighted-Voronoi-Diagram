cmake_minimum_required(VERSION 3.10)

set (CMAKE_CXX_STANDARD 17)
set (CMAKE_CXX_FLAGS "-Wall -Wextra")
set (CMAKE_CXX_FLAGS_DEBUG "-g -O0")
set (CMAKE_CXX_FLAGS_RELWITHDEBUGINFO "-g -O3")
set (CMAKE_CXX_FLAGS_RELEASE "-O3")

# set the project name
project(awVd_ray LANGUAGES CXX)

# boost program_options uses std11 strings
add_compile_definitions(_GLIBCXX_USE_CXX11_ABI=1)

find_package( Boost COMPONENTS program_options REQUIRED)

find_package( Eigen3 3.3 REQUIRED NO_MODULE)

if ( NOT Boost_FOUND )

  message(STATUS "This project requires the Boost library, and will not be compiled.")

  return()

endif()

# include for local directory

# include for local package
include_directories(../ESBTL-1.0-beta01/include)
include_directories(../overlap)
include_directories(../cgal-5.4/CGAL_Core/include)

# add the executable
add_executable(vdos
  ${PROJECT_SOURCE_DIR}/src/main.cpp
  ${PROJECT_SOURCE_DIR}/src/power.cpp
  ${PROJECT_SOURCE_DIR}/src/awVd_ray.cpp
  ${PROJECT_SOURCE_DIR}/src/compute_bisectors.cpp
  ${PROJECT_SOURCE_DIR}/src/ray_intersections.cpp
  ${PROJECT_SOURCE_DIR}/src/voronoi_mesh_details.cpp
  ${PROJECT_SOURCE_DIR}/src/bfs_for_vert.cpp)

target_link_libraries(vdos PUBLIC CGAL::CGAL Eigen3::Eigen PRIVATE Boost::program_options)
