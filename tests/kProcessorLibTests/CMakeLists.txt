file(GLOB SRCS *.cpp)

ADD_EXECUTABLE(testKprocessorLib ${SRCS})

TARGET_LINK_LIBRARIES(testKprocessorLib kProcessor  libgtest)

add_test(NAME testKprocessorLib
         COMMAND testKprocessorLib)
