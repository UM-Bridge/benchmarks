cmake_minimum_required (VERSION 3.26)

project(OndoMathX/Ocean)

#Specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

#Compile options
set(CMAKE_CXX_FLAGS "-O3")

#Include the main source directory
include_directories(${CMAKE_SOURCE_DIR}/Sources)

#Include the OndoMathX sources
include_directories(${CMAKE_SOURCE_DIR}/../core/Sources)

#Include Eigen
include_directories(${CMAKE_SOURCE_DIR}/../core/External/eigen-3.4.0)

#Include Eigen
include_directories(${CMAKE_SOURCE_DIR}/../../../umbridge/lib)

#Include OpenMP
find_package(OpenMP REQUIRED)
include_directories(${OpenMP_CXX_INCLUDE_DIRS})

#Link OpenMP
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
link_libraries(OpenMP::OpenMP_CXX)
 

#Make directory for data outputs
make_directory(${CMAKE_BINARY_DIR}/Data)

#Subdirectory#######################

add_subdirectory (Sources/)

####################################

 