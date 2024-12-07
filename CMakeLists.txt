cmake_minimum_required(VERSION 3.10)
project(KF-GINS)

set(CMAKE_CXX_STANDARD 14)

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)
set(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/lib)

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
    add_definitions( -D_USE_MATH_DEFINES=True)
    add_compile_options("$<$<C_COMPILER_ID:MSVC>:/utf-8>")
    add_compile_options("$<$<CXX_COMPILER_ID:MSVC>:/utf-8>")
    add_compile_options(-std=c++14 -pthread -g  -ggdb -O0 -Wall   -fstack-protector-all -march=native -fnon-call-exceptions)
else()
    set(CMAKE_CXX_FLAGS "-Wno-unused-result")
endif()

file(GLOB_RECURSE SOURCE
        src/new_kf.cpp
        src/fileio/fileloader.cc
        src/fileio/filesaver.cc
        src/kf-gins/gi_engine.cpp
        src/kf-gins/insmech.cpp
        src/kf-gins/drawError.cpp
        src/kf-gins/satellite.cpp
        src/kf-gins/DataProcessing.cpp
        )
message(STATUS "All source files: ${SOURCE}")
include_directories(${PROJECT_SOURCE_DIR}/src)
add_executable(${PROJECT_NAME} ${SOURCE})

# python3
set(Python3_INCLUDE_DIRS "/usr/include/python3.10")
set(Python3_LIBRARIES "/usr/lib/x86_64-linux-gnu/libpython3.10.so")
include_directories(${Python3_INCLUDE_DIRS})
link_directories(${Python3_LIBRARIES})
target_link_libraries(${PROJECT_NAME} ${Python3_LIBRARIES})

# matplotlib-cpp
include_directories("/home/sxh/matplotlib-cpp")

# Eigen3
include_directories(ThirdParty/eigen-3.3.9)

# yaml-cpp-0.7.0
add_subdirectory(ThirdParty/yaml-cpp-0.7.0)
target_link_libraries(${PROJECT_NAME} yaml-cpp)

# abseil
set(ABSL_PROPAGATE_CXX_STD true)
add_subdirectory(ThirdParty/abseil-cpp-20220623.1)
target_link_libraries(${PROJECT_NAME}
        absl::strings
        absl::str_format
        absl::time)

#pybind11
set(Pybind11_INCLUDE_DIRS "/usr/local/include/pybind11")
find_package(pybind11 REQUIRED)
target_link_libraries(${PROJECT_NAME} pybind11::module)
target_include_directories(${PROJECT_NAME} PRIVATE "/usr/local/include/pybind11")
