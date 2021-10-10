add_executable(kCorticall main.cpp)
target_link_libraries (kCorticall LINK_PUBLIC kProcessor MQF ${SEQAN_LIBRARIES} )


file(COPY Snakefile DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY env.yaml DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY project_config.yaml  DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY config.yaml  DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY sample_table.csv  DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY subsample_table.csv  DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY download.smk  DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY urls  DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
