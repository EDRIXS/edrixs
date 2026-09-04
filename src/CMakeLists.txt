set(src
  m_constants.f90
  m_control.f90
  m_types.f90
  m_global.f90
  utils.f90
  io.f90
  fock.f90
  spmv.f90
  linsys.f90
  ham.f90
  full_diag.f90
  lanczos.f90
  arpack.f90
  ed_driver.f90
  xas_driver.f90
  rixs_driver.f90
  opavg_driver.f90)

add_library(edrixs STATIC ${src})
target_link_libraries(edrixs PUBLIC MPI::MPI_Fortran LAPACK::LAPACK PARPACK::PARPACK)

install(TARGETS edrixs DESTINATION lib)

set(progs ed xas rixs opavg)
foreach(prog ${progs})
  add_executable(${prog}.x ${prog}_main.f90)
  target_link_libraries(${prog}.x PUBLIC edrixs)
  install(TARGETS ${prog}.x DESTINATION bin)
endforeach(prog)

if (EDRIXS_PY_INTERFACE)
  # Vars
  set(f2py_module_name fedrixs)
  set(fortran_src_file "${CMAKE_CURRENT_SOURCE_DIR}/pyapi.f90")
  set(f2py_module_c "${f2py_module_name}module.c")
  set(generated_module_file "${f2py_module_name}.${Python3_SOABI}")

  message(STATUS "generated_module_file: ${generated_module_file}")

  add_custom_command(
    OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/${f2py_module_c}"
    COMMAND ${Python3_EXECUTABLE} -m "numpy.f2py"
                     "${fortran_src_file}"
                     -m ${f2py_module_name}
                     --lower # Important
    DEPENDS ${fortran_src_file} # Fortran source
  )

  # Set up target
  add_library(${f2py_module_name} SHARED
    "${CMAKE_CURRENT_BINARY_DIR}/${f2py_module_c}" # Generated
    "${F2PY_INCLUDE_DIR}/fortranobject.c" # From NumPy
    "${fortran_src_file}" # Fortran source(s)
    )

  target_include_directories(${f2py_module_name} PUBLIC ${Python3_INCLUDE_DIRS} ${F2PY_INCLUDE_DIR} ${Python3_NumPy_INCLUDE_DIRS})
  target_link_libraries(${f2py_module_name} PUBLIC edrixs)

  if (LINUX)
    set_target_properties(
        ${f2py_module_name}
        PROPERTIES
            PREFIX ""
            OUTPUT_NAME ${generated_module_file}
            LINKER_LANGUAGE C
      )
  elseif(APPLE)
    set_target_properties(
      ${f2py_module_name}
      PROPERTIES
          PREFIX ""
          OUTPUT_NAME ${generated_module_file}
          SUFFIX ".so"
          LINKER_LANGUAGE C
          LINK_FLAGS "-undefined dynamic_lookup"
    )
  endif()
endif()
