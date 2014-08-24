
add_executable(laplace laplace/laplace.cpp)

include_directories(${UMFPACK_INCLUDES})
include_directories (${MiniFPM_SOURCE_DIR}/src) 
target_link_libraries(laplace mfpm ${UMFPACK_LIBRARIES} ${LAPACK_LIBRARIES})

install(TARGETS laplace
        RUNTIME DESTINATION ${MiniFPM_SOURCE_DIR}/examples/laplace
       )

add_executable(mixer mixer/mixer.cpp)

include_directories(${UMFPACK_INCLUDES})
include_directories (${MiniFPM_SOURCE_DIR}/src) 
target_link_libraries(mixer mfpm ${UMFPACK_LIBRARIES} ${LAPACK_LIBRARIES})

install(TARGETS mixer
        RUNTIME DESTINATION ${MiniFPM_SOURCE_DIR}/examples/mixer
       )

