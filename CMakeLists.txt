cmake_minimum_required(VERSION 3.13)
project(FHN)

SET(CMAKE_CXX_FLAGS "-lm -lfftw3")
SET(CMAKE_C_FLAGS   "-lm -lfftw3")

set(CMAKE_CXX_STANDARD 17)

add_executable(FHN main.cpp src/ode.h src/ode.cpp src/ArgParse.cpp src/ArgParse.h src/Argument.cpp src/Argument.h src/npy.h)
target_link_libraries(FHN fftw3)