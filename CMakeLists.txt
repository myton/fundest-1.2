# fundest CMake file; see http://www.cmake.org and 
#                        http://www.insightsoftwareconsortium.org/wiki/index.php/CMake_Tutorial

CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
PROJECT(FUNDEST)

# set default build mode for single-configuration generator
IF(UNIX OR CYGWIN)
  IF(NOT CMAKE_BUILD_TYPE)
    SET(CMAKE_BUILD_TYPE Release CACHE STRING
      "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel" FORCE)
  ENDIF(NOT CMAKE_BUILD_TYPE)
ENDIF(UNIX OR CYGWIN)

IF(CMAKE_GENERATOR MATCHES "^Visual Studio *")
  # get rid of CRT warnings
  ADD_DEFINITIONS(-D_CRT_SECURE_NO_WARNINGS)
ENDIF(CMAKE_GENERATOR MATCHES "^Visual Studio *")


# levmar library
SET(LEVMAR_INCDIR "/usr/local/include/levmar" CACHE PATH "Path to LEVMAR library header")
SET(LEVMAR_LIBDIR "/usr/local/lib/liblevmar.a" CACHE PATH "Path to LEVMAR library")

message(status "LEVMAR_INCDIR: " ${LEVMAR_INCDIR})
message(status "LEVMAR_LIBDIR: " ${LEVMAR_LIBDIR})

# lapack/blas
SET(HAVE_LAPACK 1 CACHE BOOL "Do we have LAPACK/BLAS?")
# the directory where the lapack/blas/f2c libraries reside
SET(LAPACKBLAS_DIR "/usr/lib" CACHE PATH "Path to lapack/blas libraries")
#SET(LAPACKBLAS_DIR "C:\Program Files\Intel\MKL\10.2.2.025\ia32\lib" CACHE PATH "Path to lapack/blas libraries") # MKL
SET(NEED_F2C 1 CACHE BOOL "Do we need either f2c or F77/I77?")
SET(HAVE_PLASMA 0 CACHE BOOL "Do we have PLASMA parallel linear algebra library?")
IF(HAVE_PLASMA)
 SET(PLASMA_DIR "/usr/local/PLASMA" CACHE PATH "Path to PLASMA root")
ENDIF(HAVE_PLASMA)
OPTION(BUILD_DEMO "Build demo program?" TRUE)

# actual names for the lapack/blas/f2c libraries
SET(LAPACKBLAS_LIB_NAMES "lapack;blas" CACHE STRING "The name of the lapack & blas libraries")
#SET(LAPACKBLAS_LIB_NAMES "mkl_solver_sequential;mkl_intel_c;mkl_sequential;mkl_core" CACHE STRING "The name of the lapack libraries") # MKL
IF(NEED_F2C)
  SET(F2C_LIB_NAME f2c CACHE STRING "The name of the f2c or F77/I77 library")
  # f2c is sometimes equivalent to libF77 & libI77
  #SET(F2C_LIB_NAME "libF77;libI77" CACHE STRING "The name of the f2c or F77/I77 library")
ELSE(NEED_F2C)
  SET(F2C_LIB_NAME "" CACHE STRING "The name of the f2c or F77/I77 library")
ENDIF(NEED_F2C)

# actual names for the PLASMA libraries
IF(HAVE_PLASMA)
  SET(PLASMA_LIB_NAMES "plasma;coreblas;quark;lapacke" CACHE STRING "The names of the PLASMA libraries")
ENDIF(HAVE_PLASMA)

SET(HAVE_HOMEST 0 CACHE BOOL "Do we have HOMEST installed? (for GRIC only)")
IF(HAVE_HOMEST)
 SET(HOMEST_INCDIR "/usr/local/homest-1.4" CACHE PATH "Path to HOMEST library header")
 SET(HOMEST_LIBDIR "/usr/local/homest-1.4" CACHE PATH "Path to HOMEST library")
 ADD_DEFINITIONS(-DHAVE_HOMEST)
ENDIF(HAVE_HOMEST)


########################## NO CHANGES BEYOND THIS POINT ##########################

INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR})
INCLUDE_DIRECTORIES(${LEVMAR_INCDIR})
INCLUDE_DIRECTORIES(${HOMEST_INCDIR})
IF(HAVE_PLASMA)
 INCLUDE_DIRECTORIES(${PLASMA_DIR}/include ${PLASMA_DIR}/quark)
ENDIF(HAVE_PLASMA)

# fundest library source files
ADD_LIBRARY(fundest STATIC
  calc_fund_coeffs.c fundest.c lqs.c buckets.c norm.c linalg.c gric.c misc.c
  compiler.h calc_fund_coeffs.h fundest.h lqs.h util.h
)

# demo program
IF(BUILD_DEMO)
  IF(HAVE_HOMEST)
   SET(LIBS homest)
   LINK_DIRECTORIES(${HOMEST_LIBDIR})
  ENDIF(HAVE_HOMEST)

  SET(LIBS ${LIBS} fundest levmar)

  LINK_DIRECTORIES(${CMAKE_BINARY_DIR}) # location of the fundest library
  LINK_DIRECTORIES(${LAPACKBLAS_DIR})
  LINK_DIRECTORIES(${LEVMAR_LIBDIR})

# libraries the demo depends on
  IF(HAVE_PLASMA)
    LINK_DIRECTORIES(${PLASMA_DIR}/lib)
    SET(LIBS ${LIBS} ${PLASMA_LIB_NAMES})
  ENDIF(HAVE_PLASMA)

  IF(HAVE_LAPACK)
    IF(NEED_F2C)
      SET(LIBS ${LIBS} ${LAPACKBLAS_LIB_NAMES} ${F2C_LIB_NAME})
    ELSE(NEED_F2C)
      SET(LIBS ${LIBS} ${LAPACKBLAS_LIB_NAMES})
    ENDIF(NEED_F2C)
  ENDIF(HAVE_LAPACK)

  ADD_EXECUTABLE(fundest_demo fundest_demo.c fundest.h levmar.h)
  TARGET_LINK_LIBRARIES(fundest_demo ${LIBS} -lm)
  MESSAGE(STATUS "fundest_demo will be linked against ${LIBS}")

# make sure that the library is built before the demo
  ADD_DEPENDENCIES(fundest_demo fundest)
ENDIF(BUILD_DEMO)

#ADD_SUBDIRECTORY(matlab)
