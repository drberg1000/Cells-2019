cmake_minimum_required(VERSION 3.15)
project(Cells)

set(CMAKE_CXX_STANDARD 14)

include_directories(src)

add_executable(Cells
        src/cartesian.c
        src/cartesian.h
        src/Cells.c
        src/Cells.h
        src/mersenne.c
        src/mersenne.h
        src/Queue.c
        src/Queue.h
        src/Simulation.c)
target_link_libraries(Cells m)

add_executable(Test
        src/cartesian.c
        src/cartesian.h
        src/Cells-test.c
        src/Cells.c
        src/Cells.h
        src/mersenne.c
        src/mersenne.h
        src/Queue.c
        src/Queue.h)
target_link_libraries(Test m)
