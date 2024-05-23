cmake_minimum_required(VERSION 3.22)

project (quicksched C CXX)
message(STATUS "CMake system name is ${CMAKE_SYSTEM_NAME}")

# ------------------------------------------------------------------------------
# Policies - use the latest of everything
# ------------------------------------------------------------------------------
cmake_policy(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})

#------------------------------------------------------------------------------
# Options
#------------------------------------------------------------------------------
OPTION (BUILD_SHARED_LIBS "Build shared libraries." ON)

#------------------------------------------------------------------------------
# Setup Output directories
#------------------------------------------------------------------------------
SET (CMAKE_RUNTIME_OUTPUT_DIRECTORY
    ${PROJECT_BINARY_DIR}/bin CACHE PATH "Single Directory for all Executables."
)
SET (CMAKE_LIBRARY_OUTPUT_DIRECTORY
    ${PROJECT_BINARY_DIR}/bin CACHE PATH "Single Directory for all Libraries"
)
SET (CMAKE_ARCHIVE_OUTPUT_DIRECTORY
    ${PROJECT_BINARY_DIR}/bin CACHE PATH "Single Directory for all static libraries."
)

#------------------------------------------------------------------------------
# Custom cmake modules should go here
#------------------------------------------------------------------------------
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")

#------------------------------------------------------------------------------
# Enable testing
#------------------------------------------------------------------------------
enable_testing()
include(CTest)

#------------------------------------------------------------------------------
# It's a C project, but we will add hpx/pika c++ support 
#------------------------------------------------------------------------------
set(CMAKE_CXX_STANDARD 20)

#------------------------------------------------------------------------------
# MPI 
#------------------------------------------------------------------------------
find_package(MPI REQUIRED)
message(STATUS "MPI found version ${MPI_CXX_VERSION}")
list(APPEND QUICKSCHED_LIBS MPI::MPI_C)

# ------------------------------------------------------------------------------
# Threads
# ------------------------------------------------------------------------------
find_package(Threads REQUIRED)
set(HAVE_PTHREAD 1)
message(STATUS "Threads found version ${Threads_VERSION}")
list(APPEND QUICKSCHED_LIBS Threads::Threads)

# ------------------------------------------------------------------------------
# OpenMP
# ------------------------------------------------------------------------------
find_package(OpenMP COMPONENTS C QUIET)
if(OpenMP_C_FOUND)
  set(HAVE_OPENMP 1)
  message(STATUS "OpenMP found version ${OpenMP_C_VERSION}")
  list(APPEND QUICKSCHED_LIBS OpenMP::OpenMP_C)
endif()

# ------------------------------------------------------------------------------
# lapack/blas
# ------------------------------------------------------------------------------
set(BLA_VENDOR OpenBLAS)
find_package(BLAS REQUIRED)
message(STATUS "BLAS found version ${BLAS_VERSION}")
list(APPEND QUICKSCHED_LIBS BLAS::BLAS)

#------------------------------------------------------------------------------
# Generate "config.h" from "cmake/config.h.in"
#------------------------------------------------------------------------------
include(CheckEnv)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cmake/config.h.in
  ${CMAKE_CURRENT_BINARY_DIR}/config.h)
 
#------------------------------------------------------------------------------
# Generate "version.h" from "src/version.h.in"
#------------------------------------------------------------------------------
find_package(Git)

if(GIT_EXECUTABLE)
  # git commit hash macro
  execute_process(
    COMMAND ${GIT_EXECUTABLE} describe --abbrev=8  --always --tags --dirty
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    OUTPUT_VARIABLE GIT_REVISION
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )
endif()

if(GIT_REVISION STREQUAL "")
  set(GIT_REVISION 0.0.0-unknown)
  message(WARNING "Failed to determine version from Git tags. Using default version \"${FOO_VERSION}\".")
endif()

configure_file(src/version.h.in ${CMAKE_CURRENT_BINARY_DIR}/version.h @ONLY)  

#------------------------------------------------------------------------------
# main library
#------------------------------------------------------------------------------
add_library(libquicksched SHARED 
  ./src/queue.c
  ./src/qsched.c
  ./src/atomic.h
  ./src/task.h
  ./src/queue.h
  ./src/qsched.h
  ./src/cycle.h
  ./src/quicksched.h
  ./src/error.h
  ./src/res.h
  ./src/lock.h
)
target_include_directories(libquicksched PUBLIC 
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_BINARY_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/src
)
target_compile_definitions(libquicksched PUBLIC TIMERS HAVE_CONFIG_H HAVE_PTHREAD)
target_link_libraries(libquicksched PUBLIC ${QUICKSCHED_LIBS} m) 
target_include_directories(libquicksched PUBLIC
  ${CMAKE_CURRENT_BINARY_DIR}
  ${CMAKE_CURRENT_BINARY_DIR}/src
)

#------------------------------------------------------------------------------
# Tests
#------------------------------------------------------------------------------
set(TEST_SRCS
  test_bh_sorted.c
  test_qr.c
  test_bh.c
  test_fmm_sorted.c
  test.c
)
set(TEST_OMP_SRCS
  test_qr_ompss.c
)

if (0 AND OpenMP_C_FOUND)
  list(APPEND TEST_SRCS ${TEST_OMP_SRCS})
endif()

#------------------------------------------------------------------------------
# create exes for tests
#------------------------------------------------------------------------------
foreach(test ${TEST_SRCS})
  string(REPLACE ".c" "" test_name ${test})
  set(test_name qs_${test_name})

  add_executable(${test_name} ./examples/${test})
  set_target_properties(${test_name} PROPERTIES OUTPUT_NAME ${test_name})

  target_link_libraries(${test_name} PRIVATE libquicksched ${QUICKSCHED_LIBS}) 
  target_include_directories(${test_name} PRIVATE 
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_BINARY_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/src
  )

  add_test(NAME ${test_name} COMMAND ${test_name})
endforeach()
