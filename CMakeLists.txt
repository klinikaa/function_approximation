cmake_minimum_required(VERSION 3.26)
project(function_approximation)

set(CMAKE_CXX_STANDARD 17)

add_executable(function_approximation main.cpp
        approximation.cpp)
