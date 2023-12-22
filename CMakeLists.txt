cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(cifTools)

# Common settings
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_BUILD_TYPE Release)

if(LINUX)
    set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} "-fPIC")
endif()

add_definitions(-DSYSTEM_VERBOSE)
set(CMAKE_VERBOSE_MAKEFILE ON)

include_directories(
    ./
    gemmi/
)

# Specify the output directory for executables
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# cifAlignment executable
add_executable(cifAlignment
    src/c++/cifAlignment.cpp
)

# cifConverter executable
add_executable(cifConverter
    src/c++/cifConverter.cpp
)
