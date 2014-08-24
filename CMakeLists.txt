project(MiniFPM)

cmake_minimum_required(VERSION 2.6)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/modules/")
if(NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Set to Release or Debug" FORCE)
endif(NOT CMAKE_BUILD_TYPE)

find_package(umfpack REQUIRED)
find_package(LAPACK REQUIRED)
find_package(ADOLC)

add_subdirectory(src)
add_subdirectory(examples)

