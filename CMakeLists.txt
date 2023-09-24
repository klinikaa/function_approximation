cmake_minimum_required(VERSION 3.26)
project(function_approximation)

include_directories(/home/cleanyco/Downloads/matplotlib-cpp-master/, /home/cleanyco/sciplot/)

add_executable(function_approximation main.cpp
        approximation.cpp
        approximation.h
        deviation.cpp
        deviation.h
        table.cpp
        table.h
        util.h
        util.cpp
        process.cpp
        process.h
        graph.cpp
        graph.h)
