CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

#set(CMAKE_VERBOSE_MAKEFILE ON)
#set(CMAKE_EXE_LINKER_FLAGS "-t") # Enable to show linked libraries

PROJECT(Stable_Fluids)

SET(OpenCV_STATIC OFF)
FIND_PACKAGE(OpenCV REQUIRED)
INCLUDE_DIRECTORIES(${OpenCV_INCLUDE_DIRS})

SET(OpenBLAS_DIR "/usr/local/opt/openblas/lib/cmake/openblas")
FIND_PACKAGE(OpenBLAS REQUIRED)
#INCLUDE_DIRECTORIES(${OpenBLAS_INCLUDE_DIRS}) # Currently not used

FIND_PACKAGE(GSL REQUIRED 2.0)
INCLUDE_DIRECTORIES(${GSL_INCLUDE_DIRS})

set(Boost_USE_STATIC_LIBS        ON) # only find static libs
set(Boost_USE_MULTITHREADED      ON)
find_package(Boost COMPONENTS system thread REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

find_package(Threads REQUIRED)

ADD_EXECUTABLE(stable_fluids
        main.cpp main.h
        image.h
        DiffusionSolver.cpp DiffusionSolver.h
        TransportSolver.cpp TransportSolver.h
        ProjectSolver.cpp ProjectSolver.h)

# Replace GSL_LIBRARY with GLS_LIBRARIES to disable OpenBLAS
TARGET_LINK_LIBRARIES(stable_fluids ${OpenCV_LIBS} ${GSL_LIBRARY} ${OpenBLAS_LIBRARIES} ${Boost_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
