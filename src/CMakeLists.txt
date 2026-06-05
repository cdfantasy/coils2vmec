# Build Python extension via f90wrap + f2py.

find_program(F90WRAP_EXECUTABLE NAMES f90wrap)
if(F90WRAP_EXECUTABLE)
    set(F90WRAP_CMD ${F90WRAP_EXECUTABLE})
else()
    set(F90WRAP_CMD ${Python_EXECUTABLE} -m f90wrap)
endif()

set(F2PY_CMD ${Python_EXECUTABLE} -m numpy.f2py)
set(PY_PACKAGE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/coils2vmec)
set(FORTRAN_DIR ${CMAKE_CURRENT_SOURCE_DIR}/fortran)
set(FORTRAN_BINARY_DIR ${CMAKE_BINARY_DIR}/src/fortran)

set(WRAP_MODULE_NAME fieldline_tracer)
set(EXT_MODULE_NAME _fieldline_tracer)

set(F90WRAP_F90 ${CMAKE_CURRENT_BINARY_DIR}/f90wrap_fieldline_tracer_module.f90)
set(F90WRAP_PY ${CMAKE_CURRENT_BINARY_DIR}/${WRAP_MODULE_NAME}.py)
set(F2PY_STAMP ${CMAKE_CURRENT_BINARY_DIR}/f2py-build.stamp)
set(FIX_WRAPPER_SCRIPT ${CMAKE_CURRENT_BINARY_DIR}/fix_f90wrap_import.py)

file(WRITE ${FIX_WRAPPER_SCRIPT} "from pathlib import Path\nimport sys\n\npath = Path(sys.argv[1])\ntext = path.read_text()\ntext = text.replace('import _fieldline_tracer', 'from . import _fieldline_tracer')\npath.write_text(text)\n")

add_custom_command(
    OUTPUT ${F90WRAP_F90} ${F90WRAP_PY}
    COMMAND ${F90WRAP_CMD} -m ${WRAP_MODULE_NAME} ${FORTRAN_DIR}/fieldline_tracer_module.f90 --verbose
    COMMAND ${Python_EXECUTABLE} ${FIX_WRAPPER_SCRIPT} ${F90WRAP_PY}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    DEPENDS ${FORTRAN_DIR}/fieldline_tracer_module.f90
    COMMENT "Generating f90wrap interfaces"
)

add_custom_command(
    OUTPUT ${F2PY_STAMP}
    COMMAND ${CMAKE_COMMAND} -E env FC=${CMAKE_Fortran_COMPILER} ${F2PY_CMD}
            -c
            --f90flags=-I${FORTRAN_DIR}\ -I${FORTRAN_BINARY_DIR}
            --opt=-O3
            -m ${EXT_MODULE_NAME}
            ${F90WRAP_F90}
            $<TARGET_OBJECTS:fieldline_fortran_obj>
    COMMAND ${CMAKE_COMMAND} -E touch ${F2PY_STAMP}
        COMMAND_EXPAND_LISTS
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    DEPENDS ${F90WRAP_F90} fieldline_fortran_obj
    COMMENT "Building Python extension with f2py"
)

add_custom_target(coils2vmec_python_ext ALL DEPENDS ${F2PY_STAMP})

# Install pure-Python package modules from source tree.
install(
    FILES
    ${PY_PACKAGE_DIR}/__init__.py
    ${PY_PACKAGE_DIR}/fieldline.py
    ${PY_PACKAGE_DIR}/iota.py
    ${PY_PACKAGE_DIR}/lcfs.py
    ${PY_PACKAGE_DIR}/descur.py
    ${PY_PACKAGE_DIR}/descur_python.py
    ${PY_PACKAGE_DIR}/utils.py
    ${PY_PACKAGE_DIR}/plotting.py
    DESTINATION coils2vmec
)

# Install generated extension + generated helper module into package directory.
install(
    DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/
    DESTINATION coils2vmec
    FILES_MATCHING
    PATTERN "_fieldline_tracer*.so"
    PATTERN "_fieldline_tracer*.pyd"
    PATTERN "_fieldline_tracer*.dylib"
    PATTERN "fieldline_tracer.py"
)
