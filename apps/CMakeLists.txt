file(GLOB SRCS *.cpp)

ADD_EXECUTABLE(kProcessorApp ${SRCS})

TARGET_LINK_LIBRARIES(kProcessorApp kProcessor MQF)
