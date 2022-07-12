add_executable(kUpsetPlot main.cpp)
target_link_libraries (kUpsetPlot LINK_PUBLIC kProcessor MQF ${SEQAN_LIBRARIES} )


file(COPY Snakefile DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY env.yaml DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY project_config.yaml  DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY config.yaml  DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY sample_table.csv  DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY subsample_table.csv  DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY downloadData.sh  DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
