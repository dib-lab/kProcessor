# Still needs
## [ ] Build & Release options
## [ ] Build with debug option
## [ ] Build test cases option
## [ ] Build and generate static library option
## create a cmake file for build seqan
## change the static linking options to something cmake


# Minimum cmake version
cmake_minimum_required (VERSION 3.0.0)
# Name of project and that it is C/C++.
project (kProcessor C CXX)

# IMPORTANT == -FPIC #SDSL-LINKING
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

enable_testing()


include_directories("ThirdParty/CLI")
include_directories("ThirdParty/MQF/include")
include_directories("ThirdParty/json")
include_directories("ThirdParty/sdsl-lite/include")
include_directories("ThirdParty/kmerDecoder/include")
include_directories("ThirdParty/kmerDecoder/lib/parallel_hashmap")
include_directories("ThirdParty/cereal/include")
include_directories("ThirdParty/ntCard/include")
include_directories("ThirdParty/mum-hash")

set (CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS} -Wall -Wextra -std=c++17 -fPIC  -W -Wall -pedantic -O0 -g -lprofiler")

set (CMAKE_CXX_FLAGS_RELEASE  "${CMAKE_CXX_FLAGS} -Wall -Wextra -std=c++17 -fPIC  -W -Wall -pedantic -O3 -Ofast")

set (CMAKE_EXE_LINKER_FLAGS_DEBUG "${CMAKE_EXE_LINKER_FLAGS_DEBUG} -lprofiler")

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

add_subdirectory(ThirdParty/MQF/)
add_subdirectory(ThirdParty/sdsl-lite/lib/)
add_subdirectory(ThirdParty/kmerDecoder)

include(ExternalProject)
ExternalProject_Add(
  libntcard_a
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/ThirdParty/ntCard
        CONFIGURE_COMMAND  ${CMAKE_CURRENT_SOURCE_DIR}/ThirdParty/ntCard/autogen.sh && ${CMAKE_CURRENT_SOURCE_DIR}/ThirdParty/ntCard/configure --prefix=${CMAKE_CURRENT_SOURCE_DIR}/ThirdParty/ntCard/  CXXFLAGS=-fPIC
  BUILD_COMMAND make
  BUILD_IN_SOURCE 1
)
set(ntcardLib "${CMAKE_CURRENT_SOURCE_DIR}/ThirdParty/ntCard/libntcard.a")


set(LIBRARY_NAME
    kProcessor
)
file(GLOB kProcessorLibSrc "src/*.cpp" "src/kDataFrames/*.cpp" "src/KmerDecoder/*.cpp" "src/Utils/*.cpp")


add_library(${LIBRARY_NAME} STATIC
    ${kProcessorLibSrc}
    )


target_link_libraries(kProcessor MQF)
target_link_libraries (kProcessor ${SEQAN_LIBRARIES})
target_link_libraries (kProcessor ${ntcardLib}  -ldl)
#set_property(TARGET sdsl PROPERTY POSITION_INDEPENDENT_CODE ON)
target_link_libraries(kProcessor sdsl)
target_link_libraries(kProcessor kmerDecoder)

target_include_directories(${LIBRARY_NAME} PRIVATE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/${LIBRARY_NAME}>
    $<INSTALL_INTERFACE:include/${LIBRARY_NAME}>
)

target_include_directories(${LIBRARY_NAME} PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/${LIBRARY_NAME}>
    $<INSTALL_INTERFACE:include/${LIBRARY_NAME}>
)


set(LDFLAGS "-lrt -lpthread -lbz2 -lz")


add_subdirectory(tests/)
add_subdirectory(apps/)


#### Building Documentation with Doxygen

if(BUILD_DOCUMENTATION)
find_package(Doxygen)
option(BUILD_DOCUMENTATION "Create and install the HTML based API
documentation (requires Doxygen)" ${DOXYGEN_FOUND})

    if(NOT DOXYGEN_FOUND)
         message(FATAL_ERROR "Doxygen is needed to build the documentation.")
    endif()

    set(doxyfile_in ${CMAKE_CURRENT_SOURCE_DIR}/kProcessor.doxygen)
    set(doxyfile ${CMAKE_CURRENT_BINARY_DIR}/doxyfile)

    configure_file(${doxyfile_in} ${doxyfile} @ONLY)

    message("Doxygen build started.")

    add_custom_target(doc
                      COMMAND ${DOXYGEN_EXECUTABLE} ${doxyfile_in}
                      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/
                      COMMENT "Generating API documentation with Doxygen"
                      VERBATIM)

    install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/html DESTINATION     share/doc)
endif()
