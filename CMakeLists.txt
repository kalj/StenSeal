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
TARGET_LINK_LIBRARIES(wave-equation stenseal)

ADD_EXECUTABLE(matrix-build examples/matrix-build.cc)
TARGET_LINK_LIBRARIES(matrix-build stenseal)

#------------------------------------------------------------------------------
# Tests
#------------------------------------------------------------------------------

# test for upwind first derivative in 1D
ADD_EXECUTABLE(test_1d_first_derivative tests/test_1d_first_derivative.cc)
TARGET_LINK_LIBRARIES(test_1d_first_derivative stenseal)
STENSEAL_ADD_TEST(test_1d_first_derivative)

# test for upwind laplace in 1D
ADD_EXECUTABLE(test_1d_upwind_laplace tests/test_1d_upwind_laplace.cc)
TARGET_LINK_LIBRARIES(test_1d_upwind_laplace stenseal)
STENSEAL_ADD_TEST(test_1d_upwind_laplace)

# test for compact laplace in 1D
ADD_EXECUTABLE(test_1d_compact_laplace tests/test_1d_compact_laplace)
TARGET_LINK_LIBRARIES(test_1d_compact_laplace stenseal)
STENSEAL_ADD_TEST(test_1d_compact_laplace)
