set(DEPENDENCY_SOURCES
    DLSODE.f
)

add_library(fieldline_fortran_obj OBJECT
    fieldline_tracer_module.f90
    ${DEPENDENCY_SOURCES}
)

set_source_files_properties(${DEPENDENCY_SOURCES} PROPERTIES
    COMPILE_FLAGS "${F77_FLAGS_RELAX}"
)

target_include_directories(fieldline_fortran_obj PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
)
