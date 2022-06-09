cmake_minimum_required(VERSION 3.20)
project(GraduateWork-2022)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -fopenmp")

include_directories(graduatework)

add_executable(GraduateWork-2022
        graduatework/algo/interfaces/AbstractSweepMethod.h
        graduatework/algo/impl/Instrumental.cpp
        graduatework/algo/interfaces/Instrumental.h
        graduatework/Application.cpp
        graduatework/algo/impl/parallel/ParallelSweepMethod.cpp
        graduatework/algo/interfaces/parallel/ParallelSweepMethod.h
        graduatework/algo/impl/serial/SerialSweepMethod.cpp
        graduatework/algo/interfaces/serial/SerialSweepMethod.h
        graduatework/algo/impl/parallel/ParallelInstrumental.cpp
        graduatework/algo/interfaces/parallel/ParallelInstrumental.h
        graduatework/algo/impl/serial/SerialInstrumental.cpp
        graduatework/algo/interfaces/serial/SerialInstrumental.h
        graduatework/test/ParallelAlgorithmComponentTest.h
        graduatework/test/SerialAlgorithmComponentTest.h
        graduatework/test/common/Profiler.h
        graduatework/test/common/TestRunner.h
        graduatework/test/common/BaseComponentTest.h
        graduatework/test/common/InstrumentalComponentTest.h
        graduatework/task/interfaces/Task.h
        graduatework/task/interfaces/System.h
        graduatework/task/interfaces/tools/Parameters.h
        graduatework/task/interfaces/tools/Area.h
        graduatework/task/interfaces/tools/Grid.h
        graduatework/task/interfaces/InitConditions.h
        graduatework/fstreams/interfaces/FWork.h
        graduatework/task/interfaces/constants/AppConstants.h)