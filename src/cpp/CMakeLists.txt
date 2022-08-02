include(CMakePrintHelpers)

cmake_minimum_required(VERSION 3.15.0)
project (test)

enable_language(Fortran)


set(EXEC sqnm)
file(GLOB SOURCES *.f90 *.cpp)

find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})
cmake_print_variables(EIGEN3_INCLUDE_DIR)

set(libraries ${libraries} ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})


# compiles the files defined by SOURCES to generante the executable defined by EXEC
add_executable(${EXEC} ${SOURCES})
target_link_libraries(${EXEC} ${libraries})
