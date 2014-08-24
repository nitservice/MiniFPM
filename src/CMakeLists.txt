set(main_SRCS
    Array.cpp
    domainIntegration.cpp
    matrixUmfpack.cpp
    mfpmException.cpp
    mfpmMatrix.cpp
    mfpmParameters.cpp
    operators.cpp
    optimisation.cpp
    points.cpp
    pointSet.cpp
    postResult.cpp
    primitive_operators.cpp
    refinement.cpp
    scalarField.cpp
    vectorField.cpp
    writeDebugVTK.cpp
   )

set(adolc_SRCS
	Array_AD.cpp
	mfpmMatrix_AD.cpp
)
	

set(all_SRCS
    ${main_SRCS}
   )

set(incl_DIRS
	${UMFPACK_INCLUDES}
)
	
   
if(ADOLC_FOUND)
	set(all_SRCS
		${all_SRCS}
		${adolc_SRCS}
	)
	set(incl_DIRS
		${incl_DIRS}
		${ADOLC_INCLUDE_DIRS}
	)
endif(ADOLC_FOUND)


add_library(mfpm STATIC ${all_SRCS})

include_directories(${incl_DIRS})

install(TARGETS mfpm
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib
       )

INSTALL_FILES(/include .h)
