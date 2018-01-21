cmake_minimum_required (VERSION 2.8.11)

project(Wavelet)

set(CMAKE_CXX_FLAGS "-std=c++11 -g -Wall -fopenmp")

set(SOURCES
        src/Matrix.cpp
        src/Vector.cpp
        src/WaveletNetwork.cpp
        src/Wavelon.cpp
        src/WaveletFunction.cpp
        src/Utils.cpp
        )

add_executable(waveletnet src/main.cpp ${SOURCES})
