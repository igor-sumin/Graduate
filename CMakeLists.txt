cmake_minimum_required(VERSION 3.20)
project(GraduateWork-2022)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -fopenmp")

include_directories(SweepMethod)

add_executable(GraduateWork-2022
        SweepMethod/main/interfaces/AbstractSweepMethod.h
        SweepMethod/main/impl/Instrumental.cpp
        SweepMethod/main/interfaces/Instrumental.h
        SweepMethod/main/Application.cpp
        SweepMethod/main/impl/parallel/ParallelSweepMethod.cpp
        SweepMethod/main/interfaces/parallel/ParallelSweepMethod.h
        SweepMethod/main/impl/serial/SerialSweepMethod.cpp
        SweepMethod/main/interfaces/serial/SerialSweepMethod.h
        SweepMethod/main/impl/parallel/ParallelInstrumental.cpp
        SweepMethod/main/interfaces/parallel/ParallelInstrumental.h
        SweepMethod/main/impl/serial/SerialInstrumental.cpp
        SweepMethod/main/interfaces/serial/SerialInstrumental.h
        SweepMethod/test/ParallelAlgorithmComponentTest.h
        SweepMethod/test/SerialAlgorithmComponentTest.h
        SweepMethod/test/common/Profiler.h
        SweepMethod/test/common/TestRunner.h
        SweepMethod/test/common/BaseComponentTest.h
        SweepMethod/test/common/MeasureTimeComponentTest.h
        SweepMethod/main/interfaces/Task.h)