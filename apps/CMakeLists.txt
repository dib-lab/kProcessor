file(GLOB SRCS *.cpp)

#ADD_EXECUTABLE(kProcessorApp ${SRCS})
#
#TARGET_LINK_LIBRARIES(kProcessorApp kProcessor MQF)


add_executable(loadFromTSV loadFromTSV.cpp)

target_link_libraries (loadFromTSV LINK_PUBLIC kProcessor MQF ${SEQAN_LIBRARIES}  ${ntcardLib})



add_executable(index index.cpp)
target_link_libraries (index LINK_PUBLIC kProcessor MQF ${SEQAN_LIBRARIES}  ${ntcardLib})



add_executable(testIndex testIndex.cpp)
target_link_libraries (testIndex LINK_PUBLIC kProcessor MQF ${SEQAN_LIBRARIES}  ${ntcardLib})

