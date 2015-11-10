CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

SET(OpenCV_STATIC OFF)
FIND_PACKAGE(OpenCV REQUIRED)

PROJECT(Stable_Fluid)

INCLUDE_DIRECTORIES(${OpenCV_INCLUDE_DIRS})
ADD_EXECUTABLE(solver solver.cpp addSource addSource.cpp diffuse diffuse.cpp transport transport.cpp conserveMass conserveMass.cpp dissipate dissipate.cpp)