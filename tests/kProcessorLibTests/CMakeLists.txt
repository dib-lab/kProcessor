file(GLOB SRCS *.cpp)

ADD_EXECUTABLE(testKprocessorLib ${SRCS})

TARGET_LINK_LIBRARIES(testKprocessorLib kProcessor lMQF ${SEQAN_LIBRARIES} libgtest)

add_test(NAME testKprocessorLib
         COMMAND testKprocessorLib)

file(COPY ../testData/test.noN.fastq DESTINATION ${CMAKE_CURRENT_BINARY_DIR})         
