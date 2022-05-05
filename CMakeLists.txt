cmake_minimum_required(VERSION 3.1)
project(OpenDP)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0")
set(CMAKE_C_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0")

add_executable(OpenDP
        src/main.cpp
        src/parser.cpp
        src/parser_helper.cpp
        src/agent.cpp
        src/assign.cpp
        src/check_legal.cpp
        src/utility.cpp
        src/place.cpp
        src/rtree.cpp
        src/state_calc.cpp)