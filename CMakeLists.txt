# Minimum cmake version
cmake_minimum_required (VERSION 3.0.0)

set (SEQAN_INCLUDE_PATH "${CMAKE_CURRENT_LIST_DIR}/ThirdParty/seqan/include")

# Name of project and that it is C/C++.
project (kProcessor C CXX)

# ----------------------------------------------------------------------------
# Seqan building
# ----------------------------------------------------------------------------

find_package (ZLIB)
set (SeqAn_DIR "ThirdParty/seqan/share/cmake/seqan")
find_package (SeqAn REQUIRED)
include_directories (${SEQAN_INCLUDE_DIRS})
add_definitions (${SEQAN_DEFINITIONS})
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DSEQAN_HAS_ZLIB=1 -DSEQAN_HAS_BZIP2=1  -DSEQAN_HAS_OPENMP=1 -Wl,--whole-archive -Wl,--no-whole-archive")

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -std=c++14 -fPIC -fopenmp -W -Wall -pedantic")

set(LDFLAGS "-lrt -lpthread -lbz2 -lz")

LINK_DIRECTORIES(${CMAKE_BINARY_DIR}/MQF)

include_directories("ThirdParty/CLI")
include_directories("include")
include_directories("KmerDecoder")
include_directories("KmerCounter")
include_directories("HashUtils")
include_directories("Utils")
include_directories("MQF")


file(GLOB SOURCES "src/*.cpp")



add_executable(kProcessor ${SOURCES})


