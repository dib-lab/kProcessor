cmake_minimum_required (VERSION 3.0.0)

project (MQF)

set(CMAKE_BUILD_TYPE Release)

add_definitions(-std=c++14)



include_directories("include")

file(GLOB SOURCES "src/*cpp")

set(PROFILE "-pg -no-pie") # for bug in gprof.
set(DEBUG "")
set(OPT "-Ofast")

set(CMAKE_CXX_FLAGS "-fPIC -Wall ${DEBUG} ${PROFILE} ${OPT} -m64 -Wno-unused-result -Wno-strict-aliasing -Wno-unused-function -fopenmp")


add_definitions(${CMAKE_CXX_FLAGS})

add_library(lMQF STATIC ${SOURCES})

