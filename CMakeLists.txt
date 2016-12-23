CMAKE_MINIMUM_REQUIRED(VERSION 2.8.8)

PROJECT(stenseal)

# include deal.II, but don't use their macros to retain control
FIND_PACKAGE(deal.II 8.0 REQUIRED
  HINTS ${DEAL_II_DIR} ../ ../../ $ENV{DEAL_II_DIR}
  )

IF (deal.II_FOUND)
    INCLUDE_DIRECTORIES(${DEAL_II_INCLUDE_DIRS})
    IF (${CMAKE_BUILD_TYPE} STREQUAL "DEBUG")
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

# we crave 'modern' C++
SET (CMAKE_CXX_STANDARD 11)

# add our include directory to search path
INCLUDE_DIRECTORIES(include/)

# set up the library
SET (LIBNAME stenseal)

FILE(GLOB lib_srcs ${CMAKE_SOURCE_DIR}/source/*.cc)

ADD_LIBRARY(${LIBNAME} ${lib_srcs} )
TARGET_LINK_LIBRARIES(${LIBNAME} ${DEAL_II_LIB})


# Tests
ADD_EXECUTABLE(test_1d_upwind_block_operator tests/test_1d_upwind_block_operator.cc)
TARGET_LINK_LIBRARIES(test_1d_upwind_block_operator stenseal)

# an example using the library
ADD_EXECUTABLE(example example.cc)
TARGET_LINK_LIBRARIES(example stenseal)

# a wave equation solver
ADD_EXECUTABLE(wave-equation examples/wave-equation.cc)
TARGET_LINK_LIBRARIES(wave-equation stenseal)
