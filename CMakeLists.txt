cmake_minimum_required(VERSION 3.26)
project(function_approximation)

set(CMAKE_CXX_STANDARD 17)

add_executable(function_approximation main.cpp
        approximation.cpp
        approximation.h
        graph.cpp
        graph.h
        ../matplotlib-cpp/matplotlibcpp.h)

add_definitions(-D MS_WIN64)

include_directories(C:/Users/chesn/AppData/Local/Programs/Python/Python311/include 
        C:/Users/chesn/AppData/Local/Programs/Python/Python311/Lib/site-packages/numpy/core/include/
        C:/Users/chesn/AppData/Local/Programs/Python/Python311/)

#TODO change default main path