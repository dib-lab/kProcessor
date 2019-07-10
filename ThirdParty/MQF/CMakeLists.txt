cmake_minimum_required (VERSION 3.0.0)

project (MQF)

set(CMAKE_BUILD_TYPE Release)

add_definitions(-std=c++14)









file(GLOB SOURCES "src/*cpp")

set(PROFILE "-pg -no-pie") # for bug in gprof.
set(DEBUG "")
set(OPT "-Ofast")

set(CMAKE_CXX_FLAGS "-fPIC -Wall ${DEBUG} ${PROFILE} ${OPT} -m64 -Wno-unused-result -Wno-strict-aliasing -Wno-unused-function -fopenmp  -fpermissive")


add_definitions(${CMAKE_CXX_FLAGS})



## stxxl

# disallow in-source builds
if("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")
  message(SEND_ERROR "In-source builds are not allowed.")
endif("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")

# enable warnings (always good)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -W -Wall")

# include the STXXL library in my-project
add_subdirectory(stxxl)

# apply STXXL CXXFLAGS to our configuration
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${STXXL_CXX_FLAGS}")

# add STXXL includes path
include_directories(${STXXL_INCLUDE_DIRS})

get_cmake_property(_variableNames VARIABLES)
list (SORT _variableNames)
foreach (_variableName ${_variableNames})
    message(STATUS "${_variableName}=${${_variableName}}")
endforeach()


include_directories("include")


add_library(lMQF STATIC ${SOURCES})
target_link_libraries(lMQF ${STXXL_LIBRARIES})
