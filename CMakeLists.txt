
CMAKE_MINIMUM_REQUIRED(VERSION 2.8.8)

PROJECT(stenseal)

# include deal.II, but don't use their macros to retain control
FIND_PACKAGE(deal.II 8.0 REQUIRED
  HINTS ${DEAL_II_DIR} ../ ../../ $ENV{DEAL_II_DIR}
  )

STRING(TOUPPER "${CMAKE_BUILD_TYPE}" uppercase_CMAKE_BUILD_TYPE)
IF (deal.II_FOUND)
    INCLUDE_DIRECTORIES(${DEAL_II_INCLUDE_DIRS})
    IF (${uppercase_CMAKE_BUILD_TYPE} STREQUAL "DEBUG")
        LINK_DIRECTORIES(${DEAL_II_LIBRARIES_DEBUG})
        SET(DEAL_II_LIB
            ${DEAL_II_LIBRARIES_DEBUG}
            )
    ELSE()
        LINK_DIRECTORIES(${DEAL_II_LIBRARIES_RELEASE})
        SET(DEAL_II_LIB
            ${DEAL_II_LIBRARIES_RELEASE}
            )
        ENDIF()
    MESSAGE("-- DEAL.II include files found at ${DEAL_II_INCLUDE_DIRS}")
ELSE()
    MESSAGE("-- DEAL.II package not found")
ENDIF()

# we crave modern C++
SET (CMAKE_CXX_STANDARD 14)

# we also crave optimal optimizations
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -march=native -mtune=native")


# set up testing
ENABLE_TESTING()
MACRO(STENSEAL_ADD_TEST _target)
  SET(test_build_target ctest_build_${_target})
  SET(test_run_target ctest_run_${_target})
  ADD_TEST(${test_build_target} "${CMAKE_COMMAND}" --build ${CMAKE_BINARY_DIR} --target ${_target})
  ADD_TEST(${test_run_target} ${_target})
  SET_TESTS_PROPERTIES ( ${test_run_target} PROPERTIES DEPENDS ${test_build_target})
ENDMACRO(STENSEAL_ADD_TEST)

# add our include directory to search path
INCLUDE_DIRECTORIES(include/)

# set up the library
SET (LIBNAME stenseal)

FILE(GLOB lib_srcs ${CMAKE_SOURCE_DIR}/source/*.cc)

ADD_LIBRARY(${LIBNAME} ${lib_srcs} )
TARGET_LINK_LIBRARIES(${LIBNAME} ${DEAL_II_LIB})


#------------------------------------------------------------------------------
# Applications
#------------------------------------------------------------------------------

# a wave equation solver
ADD_EXECUTABLE(wave-equation examples/wave-equation.cc)
TARGET_LINK_LIBRARIES(wave-equation ${LIBNAME})

# 1d upwind operator benchmarks
ADD_EXECUTABLE(benchmarks_1d examples/benchmarks_1d.cc)
TARGET_LINK_LIBRARIES(benchmarks_1d ${LIBNAME})

#------------------------------------------------------------------------------
# Tests
#------------------------------------------------------------------------------

# test vector initialization in geometry
ADD_EXECUTABLE(test_initialize_vector tests/test_initialize_vector.cc)
TARGET_LINK_LIBRARIES(test_initialize_vector ${LIBNAME})
STENSEAL_ADD_TEST(test_initialize_vector)

# test for upwind first derivative in 1D
ADD_EXECUTABLE(test_1d_first_derivative tests/test_1d_first_derivative.cc)
TARGET_LINK_LIBRARIES(test_1d_first_derivative ${LIBNAME})
STENSEAL_ADD_TEST(test_1d_first_derivative)

# test for upwind laplace in 1D
ADD_EXECUTABLE(test_1d_upwind_laplace tests/test_1d_upwind_laplace.cc)
TARGET_LINK_LIBRARIES(test_1d_upwind_laplace ${LIBNAME})
STENSEAL_ADD_TEST(test_1d_upwind_laplace)

# test for compact laplace in 1D
ADD_EXECUTABLE(test_1d_compact_laplace tests/test_1d_compact_laplace)
TARGET_LINK_LIBRARIES(test_1d_compact_laplace ${LIBNAME})
STENSEAL_ADD_TEST(test_1d_compact_laplace)

# test for upwind laplace in 1D for general geometry
ADD_EXECUTABLE(test_1d_general_upwind_laplace tests/test_1d_general_upwind_laplace.cc)
TARGET_LINK_LIBRARIES(test_1d_general_upwind_laplace ${LIBNAME})
STENSEAL_ADD_TEST(test_1d_general_upwind_laplace)

# test for compact laplace in 1D for general geometry
ADD_EXECUTABLE(test_1d_general_compact_laplace tests/test_1d_general_compact_laplace.cc)
TARGET_LINK_LIBRARIES(test_1d_general_compact_laplace ${LIBNAME})
STENSEAL_ADD_TEST(test_1d_general_compact_laplace)

# test for coefficient matrix upwind laplance
ADD_EXECUTABLE(test_1d_upwind_laplace_matrix tests/test_1d_upwind_laplace_matrix.cc)
TARGET_LINK_LIBRARIES(test_1d_upwind_laplace_matrix ${LIBNAME})
STENSEAL_ADD_TEST(test_1d_upwind_laplace_matrix)

# test for compact laplace in 1D for general geometry
ADD_EXECUTABLE(test_1d_compact_laplace_matrix tests/test_1d_compact_laplace_matrix.cc)
TARGET_LINK_LIBRARIES(test_1d_compact_laplace_matrix ${LIBNAME})
STENSEAL_ADD_TEST(test_1d_compact_laplace_matrix)

# test for quadrature in 1D 
ADD_EXECUTABLE(test_1d_quadrature tests/test_1d_quadrature.cc)
TARGET_LINK_LIBRARIES(test_1d_quadrature ${LIBNAME})
STENSEAL_ADD_TEST(test_1d_quadrature)