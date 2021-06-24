file(GLOB SRCS *.cpp)

ADD_EXECUTABLE(testKprocessorLib ${SRCS})



TARGET_LINK_LIBRARIES(testKprocessorLib kProcessor MQF libgtest ${ntcardLib} )

add_test(NAME testKprocessorLib
         COMMAND testKprocessorLib)

file(COPY ../testData/test.noN.fastq DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY ../testData/test2.noN.fastq DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY ../testData/test3.noN.fastq DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY ../testData/test1.fa DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY ../testData/test1.fa.names DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY ../testData/test.noN.fastq.ntCardRes DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY ../testData/test2.noN.fastq.ntCardRes DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY ../testData/test.noN.fastq.21.kmc_pre DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY ../testData/test.noN.fastq.21.kmc_suf DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY ../testData/test.noN.fastq.31.kmc_pre DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY ../testData/test.noN.fastq.31.kmc_suf DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY ../testData/colorColumn.mixVector DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY ../testData/.stxxl DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
