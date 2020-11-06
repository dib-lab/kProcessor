add_executable(kDifferntialExpression main.cpp)
target_link_libraries (kDifferntialExpression LINK_PUBLIC kProcessor MQF ${SEQAN_LIBRARIES} )


file(COPY Snakefile DESTINATION ${CMAKE_CURRENT_BINARY_DIR})