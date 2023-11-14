cmake_minimum_required(VERSION 3.26)
project(TOSN)

set(CMAKE_CXX_STANDARD 14)

include_directories(src)
include_directories(src/core)

add_executable(TOSN
        src/core/algorithms.cpp
        src/core/algorithms.h
        src/core/deployment.cpp
        src/core/deployment.h
        src/core/sensor.cpp
        src/core/sensor.h
        src/input.cpp
        src/input.h
        src/main.cpp)
