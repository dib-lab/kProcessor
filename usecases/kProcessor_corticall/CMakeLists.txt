add_executable(cleanDBG clean.cpp)
target_link_libraries (cleanDBG LINK_PUBLIC kProcessor MQF ${SEQAN_LIBRARIES} )

add_executable(addLinks addLinks.cpp)
target_link_libraries (addLinks LINK_PUBLIC kProcessor MQF ${SEQAN_LIBRARIES} )

add_executable(detectNovelKmers detectNovelKmers.cpp)
target_link_libraries (detectNovelKmers LINK_PUBLIC kProcessor MQF ${SEQAN_LIBRARIES}  )


file(COPY Snakefile DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY SnakefileCorticall DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY env.yaml DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY project_config.yaml  DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY config.yaml  DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY sample_table.csv  DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY subsample_table.csv  DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY download.smk  DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY urls  DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
