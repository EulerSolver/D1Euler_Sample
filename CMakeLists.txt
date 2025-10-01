cmake_minimum_required(VERSION 3.10)

if(CMAKE_SYSTEM_NAME STREQUAL "Windows")
    set(CMAKE_GENERATOR "MinGW Makefiles" CACHE INTERNAL "")
endif()


project(ProjectName)


set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)


set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g")
set(CMAKE_CXX_FLAGS_RELEASE "-O2 -DNDEBUG")


include_directories(${PROJECT_SOURCE_DIR}/headers)


set(SOURCES
    src/InitialCDs.cpp #
    src/boundaryCDs.cpp #
    src/SolverEuler.cpp #
    src/TVDRK3.cpp #
    src/WENO5JS.cpp #
    src/CharactMatrixCompute.cpp #
    src/GLF_FVS.cpp #
    src/FourMethodsForTests.cpp #
    src/PhysicsValueAndControl.cpp #
    src/Main.cpp   #
    src/RunEuler.cpp   #
    src/outputToFile.cpp    #
)

add_executable(exe ${SOURCES})