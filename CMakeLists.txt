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


enable_testing()


include_directories("ThirdParty/CLI")
include_directories("ThirdParty/MQF/include")
include_directories("ThirdParty/json")
include_directories("ThirdParty/sdsl-lite/include")


add_subdirectory(ThirdParty/MQF/)
add_subdirectory(ThirdParty/sdsl-lite/)


set(LIBRARY_NAME
    kProcessor
)
file(GLOB kProcessorLibSrc "src/*.cpp" "src/KmerDecoder/*.cpp"
"src/HashUtils/*.cpp" "src/Utils/*.cpp")


add_library(${LIBRARY_NAME} STATIC
    ${kProcessorLibSrc}
    )
target_link_libraries(kProcessor lMQF)
target_link_libraries (kProcessor ${SEQAN_LIBRARIES})
target_link_libraries(kProcessor sdsl)

target_include_directories(${LIBRARY_NAME} PRIVATE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/${LIBRARY_NAME}>
    $<INSTALL_INTERFACE:include/${LIBRARY_NAME}>
)

target_include_directories(${LIBRARY_NAME} PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/${LIBRARY_NAME}>
    $<INSTALL_INTERFACE:include/${LIBRARY_NAME}>
)




# ----------------------------------------------------------------------------
# Seqan building
# ----------------------------------------------------------------------------
set (SEQAN_INCLUDE_PATH "${CMAKE_CURRENT_LIST_DIR}/ThirdParty/seqan/include")
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







add_subdirectory(tests/)
add_subdirectory(apps/)


#### Building Documentation with Doxygen
find_package(Doxygen)
option(BUILD_DOCUMENTATION "Create and install the HTML based API
documentation (requires Doxygen)" ${DOXYGEN_FOUND})

if(BUILD_DOCUMENTATION)
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
