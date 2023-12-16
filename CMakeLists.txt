cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(cifTools)

# Common settings
set(CMAKE_CXX_STANDARD 17)
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
    libs/
    gemmi/
)

# cifAlignment executable
add_executable(cifAlignment
    cifAlignment.cpp
)

# cifConverter executable
add_executable(cifConverter
    cifConverter.cpp
    libs/cParameters.cpp
    libs/cTimer.cpp
    libs/cInfo.cpp
)