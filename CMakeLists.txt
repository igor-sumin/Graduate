cmake_minimum_required(VERSION 3.20)
project(GraduateWork-2022)

set(CMAKE_CXX_STANDARD 14)

include_directories(SweepMethod)

add_executable(GraduateWork-2022
        SweepMethod/main/interfaces/AbstractSweepMethod.h
        SweepMethod/main/impl/tools/Instrumental.cpp
        SweepMethod/main/interfaces/tools/Instrumental.h
        SweepMethod/main/main.cpp
        SweepMethod/main/impl/parallel/ParallelSweepMethod.cpp
        SweepMethod/main/interfaces/parallel/ParallelSweepMethod.h
        SweepMethod/test/Profiler.h
        SweepMethod/main/impl/serial/SerialSweepMethod.cpp
        SweepMethod/main/interfaces/serial/SerialSweepMethod.h
        SweepMethod/test/TestRunner.h
        SweepMethod/test/UnitTests.h)
