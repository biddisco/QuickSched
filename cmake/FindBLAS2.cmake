include(Utilities)

# ------------------------------------------------------------------------------
# Find blas using default cmake logic
# ------------------------------------------------------------------------------
find_package(BLAS REQUIRED)
message(STATUS "BLAS found version ${BLAS_VERSION}")
# print_target_properties(BLAS::BLAS)

# ------------------------------------------------------------------------------
# Try to add include directories to make life easier
# ------------------------------------------------------------------------------
get_target_property(BLAS_LIBS BLAS::BLAS INTERFACE_LINK_LIBRARIES)
foreach(lib ${BLAS_LIBS})
  get_filename_component(LIB_DIR "${lib}" DIRECTORY)
  foreach(test_path "include" "../include" "../../include")
    unset(BLAS_INCLUDE_DIRS)
    # message("Checking ${LIB_DIR}/${test_path}")
    find_path(BLAS_INCLUDE_DIRS
      NAMES "cblas.h" "mkl_cblas.h" "mkl_cblas_64.h"
      PATHS "${LIB_DIR}/${test_path}"
      NO_DEFAULT_PATH NO_CACHE)
    if (BLAS_INCLUDE_DIRS)
      # message("GOT:${BLAS_INCLUDE_DIRS}")
      break()
    endif()
  endforeach()
  if (BLAS_INCLUDE_DIRS)
    break()
  endif()
endforeach()

if (BLAS_INCLUDE_DIRS)
  set_target_properties(BLAS::BLAS PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "${BLAS_INCLUDE_DIRS}"
  )
  message(STATUS "BLAS includes FOUND: ${BLAS_INCLUDE_DIRS}")
else()
  message(STATUS "BLAS includes NOT-FOUND: ${BLAS_INCLUDE_DIRS}")
endif()
