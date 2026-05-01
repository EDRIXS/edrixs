cmake_minimum_required(VERSION 3.17.3 FATAL_ERROR) # 3.17 > for Python3_SOABI

project (EDRIXS C Fortran)

if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING "Build type configuration" FORCE)
  message(STATUS "Setting default build type: ${CMAKE_BUILD_TYPE}")
endif()

message(STATUS "CMAKE_COMMAND: ${CMAKE_COMMAND}")
message(STATUS "CMAKE_BUILD_TYPE: ${CMAKE_BUILD_TYPE}")

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

find_package(MPI REQUIRED)

find_package(LAPACK REQUIRED)

option(EDRIXS_PYINTERFACE "Build python interface" OFF)

# see https://gitlab.kitware.com/cmake/cmake/-/issues/21779
if (CMAKE_VERSION VERSION_LESS 3.20)
  if(LAPACK_FOUND)
    set(_lapack_libs "${LAPACK_LIBRARIES}")
    if(_lapack_libs AND TARGET BLAS::BLAS)
      # remove the ${BLAS_LIBRARIES} from the interface and replace it
      # with the BLAS::BLAS target
      list(REMOVE_ITEM _lapack_libs "${BLAS_LIBRARIES}")
      list(APPEND _lapack_libs BLAS::BLAS)
    endif()
    if(_lapack_libs)
      set_target_properties(LAPACK::LAPACK PROPERTIES
        INTERFACE_LINK_LIBRARIES "${_lapack_libs}"
      )
    endif()
    unset(_lapack_libs)
  endif()
endif()

# use a custom find module because arpack-ng config module doesn't provide full paths to libs
# see https://github.com/opencollab/arpack-ng/pull/311
find_package(arpack MODULE REQUIRED COMPONENTS serial parallel)


if (EDRIXS_PY_INTERFACE)
  message(STATUS "Building python interface")
  find_package(Python3 3.7 REQUIRED COMPONENTS Interpreter Development NumPy)

  # Grab the variables from a local Python installation
  # F2PY headers
  execute_process(
    COMMAND "${Python3_EXECUTABLE}"
    -c "import numpy.f2py; print(numpy.f2py.get_include())"
    OUTPUT_VARIABLE F2PY_INCLUDE_DIR
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )

  message(STATUS "Python3_INCLUDE_DIRS: ${Python3_INCLUDE_DIRS}")
  message(STATUS "F2PY_INCLUDE_DIR: ${F2PY_INCLUDE_DIR}")
  message(STATUS "Python3_NumPy_INCLUDE_DIRS: ${Python3_NumPy_INCLUDE_DIRS}")
endif()

add_subdirectory(src)

# Optional Doxygen target: cmake --build . --target docs
find_package(Doxygen OPTIONAL_COMPONENTS dot)
if(DOXYGEN_FOUND)
    add_custom_target(docs
        COMMAND ${DOXYGEN_EXECUTABLE} ${CMAKE_SOURCE_DIR}/Doxyfile
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        COMMENT "Generating Fortran API docs with Doxygen"
        VERBATIM
    )
endif()
