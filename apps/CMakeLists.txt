#file(GLOB SRCS *.cpp)

#ADD_EXECUTABLE(kProcessorApp ${SRCS})
#
#TARGET_LINK_LIBRARIES(kProcessorApp kProcessor MQF)


add_executable(loadFromTSV loadFromTSV.cpp)

target_link_libraries (loadFromTSV LINK_PUBLIC kProcessor MQF ${SEQAN_LIBRARIES}  ${ntcardLib})


add_executable(loadFromKMC loadFromKMC.cpp)

target_link_libraries (loadFromKMC LINK_PUBLIC kProcessor MQF ${SEQAN_LIBRARIES}  ${ntcardLib})

add_executable(convert convert.cpp)
target_link_libraries (convert LINK_PUBLIC kProcessor MQF ${SEQAN_LIBRARIES}  ${ntcardLib})

add_executable(createForest createForest.cpp)
target_link_libraries (createForest LINK_PUBLIC kProcessor MQF ${SEQAN_LIBRARIES}  ${ntcardLib})


add_executable(index index.cpp)
target_link_libraries (index LINK_PUBLIC kProcessor MQF ${SEQAN_LIBRARIES}  ${ntcardLib})

add_executable(mergeIndex mergeIndex.cpp)
target_link_libraries (mergeIndex LINK_PUBLIC kProcessor MQF ${SEQAN_LIBRARIES}  ${ntcardLib})

add_executable(partitionReads_part1 partitionReads_part1.cpp)
target_link_libraries (partitionReads_part1 LINK_PUBLIC kProcessor MQF ${SEQAN_LIBRARIES}  ${ntcardLib})



add_executable(testIndex testIndex.cpp)
target_link_libraries (testIndex LINK_PUBLIC kProcessor MQF ${SEQAN_LIBRARIES}  ${ntcardLib})


add_executable(dumpMain dumpMain.cpp)
target_link_libraries (dumpMain LINK_PUBLIC kProcessor MQF ${SEQAN_LIBRARIES}  ${ntcardLib})

add_executable(optimizeIndex optimizeIndex.cpp)
target_link_libraries (optimizeIndex LINK_PUBLIC kProcessor MQF ${SEQAN_LIBRARIES}  ${ntcardLib})

add_executable(convertColorsAndtestIndex convertColorsAndtestIndex.cpp)
target_link_libraries (convertColorsAndtestIndex LINK_PUBLIC kProcessor MQF ${SEQAN_LIBRARIES}  ${ntcardLib})


add_executable(playground playground.cpp)
target_link_libraries (playground LINK_PUBLIC kProcessor MQF ${SEQAN_LIBRARIES}  ${ntcardLib})
