cmake_minimum_required(VERSION 3.10)

set (CMAKE_CXX_STANDARD 17)
set (CMAKE_CXX_FLAGS "-Wall -Wextra")
set (CMAKE_CXX_FLAGS_DEBUG "-g -O0")
set (CMAKE_CXX_FLAGS_RELWITHDEBUGINFO "-g -O3")
set (CMAKE_CXX_FLAGS_RELEASE "-O3")

# set the project name
project(awVd_ray LANGUAGES CXX CUDA)

find_package( Boost REQUIRED )

find_package( Eigen3 3.3 REQUIRED NO_MODULE)

# CGAL and its components
find_package( CGAL QUIET COMPONENTS Qt5 )

add_definitions(-DCGAL_USE_BASIC_VIEWER)

if ( NOT Boost_FOUND )

  message(STATUS "This project requires the Boost library, and will not be compiled.")

  return()

endif()

#add_library(awVd_ray SHARED
#            awVd_ray.cpp
#            compute_bisectors.cpp
#            ray_intersections.cpp
#            voronoi_mesh_details.cpp)

#target_compile_features(awVd_ray PUBLIC cxx_std_17)

#set_target_properties(awVd_ray PROPERTIES
#                            CUDA_SEPARABLE_COMPILATION ON)

# include for local directory

# include for local package
include_directories(/home/e5-2690/builds/ESBTL-1.0-beta01/include)
include_directories(/home/e5-2690/builds/eigen-3.3.9/include)
include_directories(/home/e5-2690/builds/CGAL-5.3/include)

# add the executable
add_executable(vdos
  ${PROJECT_SOURCE_DIR}/src/main.cpp
  ${PROJECT_SOURCE_DIR}/src/power.cpp
  ${PROJECT_SOURCE_DIR}/src/awVd_ray.cpp
  ${PROJECT_SOURCE_DIR}/src/compute_bisectors.cpp
  ${PROJECT_SOURCE_DIR}/src/ray_intersections.cpp
  ${PROJECT_SOURCE_DIR}/src/voronoi_mesh_details.cpp
  ${PROJECT_SOURCE_DIR}/src/bfs_for_vert.cpp)

target_link_libraries(vdos PUBLIC CGAL::CGAL Eigen3::Eigen PRIVATE )
