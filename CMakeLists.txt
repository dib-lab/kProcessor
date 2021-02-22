# Still needs
## [ ] Build & Release options
## [ ] Build with debug option
## [ ] Build test cases option
## [ ] Build and generate static library option
## change the static linking options to something cmake


# Minimum cmake version
cmake_minimum_required (VERSION 2.6.4)
# Name of project and that it is C/C++.
project (kProcessor C CXX)


#set(CMAKE_CXX_COMPILER /usr/bin/g++-10)
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
include_directories("ThirdParty/Blight")
include_directories("ThirdParty/KMC/kmc_api")


add_subdirectory(ThirdParty/MQF/)
add_subdirectory(ThirdParty/sdsl-lite/)
add_subdirectory(ThirdParty/kmerDecoder)



#set (CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS} -Wall -Wextra -std=c++17 -fPIC  -W -Wall -pedantic -O0 -g -lprofiler")

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -std=c++17 -fPIC  -W -Wall -pedantic  -O3 -Ofast")
set (CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS} -Wall -Wextra -std=c++17 -fPIC  -W -Wall -pedantic -O0 -g ")

set (CMAKE_CXX_FLAGS_RELEASE  "${CMAKE_CXX_FLAGS} -Wall -Wextra -std=c++17 -fPIC  -W -Wall -pedantic -O3 -Ofast")

set (CMAKE_EXE_LINKER_FLAGS_DEBUG "${CMAKE_EXE_LINKER_FLAGS_DEBUG}")
#set (CMAKE_EXE_LINKER_FLAGS_DEBUG "${CMAKE_EXE_LINKER_FLAGS_DEBUG} -lprofiler")

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()


include(ExternalProject)
ExternalProject_Add(
  libntcard_a
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/ThirdParty/ntCard
        CONFIGURE_COMMAND  ${CMAKE_CURRENT_SOURCE_DIR}/ThirdParty/ntCard/autogen.sh && ${CMAKE_CURRENT_SOURCE_DIR}/ThirdParty/ntCard/configure --prefix=${CMAKE_CURRENT_SOURCE_DIR}/ThirdParty/ntCard/  CXXFLAGS=-fPIC
  BUILD_COMMAND make
  BUILD_IN_SOURCE 1
)
set(ntcardLib "${CMAKE_CURRENT_SOURCE_DIR}/ThirdParty/ntCard/libntcard.a")


ExternalProject_Add(
        libblight_a
        SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/ThirdParty/Blight
        CONFIGURE_COMMAND ls
        BUILD_COMMAND make
        INSTALL_COMMAND ls
        BUILD_IN_SOURCE 1
)
set(blightLib "${CMAKE_CURRENT_SOURCE_DIR}/ThirdParty/Blight/libblight.a")


ExternalProject_Add(
        libKMCAPI_a
        SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/ThirdParty/KMC/kmc_api
        CONFIGURE_COMMAND ls
        BUILD_COMMAND make -C ../ kmc_api/kmc_file.o  kmc_api/kmer_api.o  kmc_api/mmer.o && ar -rc libKMCAPI.so  kmc_file.o  kmer_api.o  mmer.o && ranlib libKMCAPI.so
        INSTALL_COMMAND ls
        BUILD_IN_SOURCE 1
)

set(kmcLib "${CMAKE_CURRENT_SOURCE_DIR}/ThirdParty/KMC/kmc_api/libKMCAPI.so")

set(LIBRARY_NAME
    kProcessor
)

set(kProcessorLibSrc
        ${CMAKE_CURRENT_SOURCE_DIR}/src/algorithms.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/batchQuery.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/colorTable.cpp

        ${CMAKE_CURRENT_SOURCE_DIR}/src/kDataFrames/colored_kDataFrame.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/kDataFrames/defaultColumn.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/kDataFrames/kDataFrameBlight.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/kDataFrames/kDataFrameBMQF.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/kDataFrames/kDataFrame.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/kDataFrames/kDataFrameMAP.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/kDataFrames/kDataFrameMQF.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/kDataFrames/kDataFramePHMAP.cpp

        ${CMAKE_CURRENT_SOURCE_DIR}/src/Utils/utils.cpp
    )

add_library(${LIBRARY_NAME} STATIC
    ${kProcessorLibSrc}
    )


target_link_libraries(kProcessor MQF)
# TO BE REMOVED TODO V2
# target_link_libraries (kProcessor ${SEQAN_LIBRARIES})
target_link_libraries (kProcessor ${ntcardLib}  -ldl)
target_link_libraries (kProcessor ${blightLib})
target_link_libraries (kProcessor ${kmcLib})
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
add_subdirectory(usecases/kDifferentialExpression/)


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
