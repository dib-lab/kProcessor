# Still needs 
## [ ] Build & Release options
## [ ] Build with debug option
## [ ] Build test cases option
## [ ] Build and generate static library option


# Minimum cmake version
cmake_minimum_required (VERSION 3.0.0)

set (SEQAN_INCLUDE_PATH "${CMAKE_CURRENT_LIST_DIR}/ThirdParty/seqan/include")
#set(BUILD_STATIC "NO" CACHE BOOL "")
#set(... CACHE ...)

option(MyOption "BUILD_STATIC" OFF)

# Name of project and that it is C/C++.
project (kProcessor C CXX)

# ----------------------------------------------------------------------------
# Seqan building
# ----------------------------------------------------------------------------

find_package (ZLIB)

set (SeqAn_DIR "ThirdParty/seqan/share/cmake/seqan")

find_package (SeqAn REQUIRED)

set(SEQAN_INCLUDE_DIRS "${CMAKE_CURRENT_LIST_DIR}/ThirdParty/seqan/include")

include_directories (${SEQAN_INCLUDE_DIRS})

add_definitions (${SEQAN_DEFINITIONS})

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

message(STATUS "DEF: ${SEQAN_DEFINITIONS}")

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -std=c++14 -fPIC -fopenmp -W -Wall -pedantic")

message(STATUS "FLAGS: ${CMAKE_CXX_FLAGS}")


set(LDFLAGS "-lrt -lpthread -lbz2 -lz")

include_directories("ThirdParty/CLI")
include_directories("include")
include_directories("ThirdParty/MQF/include")


file(GLOB SOURCES "src/*.cpp" "lib/KmerCounter/*.cpp" "lib/KmerDecoder/*.cpp" "lib/HashUtils/*.cpp")

add_subdirectory(ThirdParty/MQF/)


if(BUILD_STATIC)
    add_library(kProcessor ${SOURCES})    
else()
    add_executable(kProcessor ${SOURCES})
endif()

target_link_libraries(kProcessor lMQF)
target_link_libraries (kProcessor ${SEQAN_LIBRARIES})





