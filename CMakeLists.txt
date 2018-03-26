cmake_minimum_required(VERSION 3.3 FATAL_ERROR)  # 3.3: EP_add(DEPENDS interfacelib

project(qcdb)
set(qcdb_AUTHORS      "Lori A. Burns")
set(qcdb_DESCRIPTION  "Quantum Chemistry Common Driver and Databases")
set(qcdb_URL          "https://github.com/qcdb/qcdb")
set(qcdb_LICENSE      "BSD 3-clause")

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

############################  Options: Build How?  #############################
include(psi4OptionsTools)
option_with_default(CMAKE_BUILD_TYPE "Build type (Release or Debug)" Release)
option_with_default(CMAKE_INSTALL_LIBDIR "Directory to which libraries installed" lib)
option_with_default(PYMOD_INSTALL_LIBDIR "Location within CMAKE_INSTALL_LIBDIR to which python modules are installed" /)
option_with_default(CMAKE_INSTALL_MESSAGE "Specify verbosity of installation messages" LAZY)

########################  Process & Validate Options  ##########################
include(autocmake_safeguards)
include(custom_color_messages)

if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set(CMAKE_INSTALL_PREFIX "/usr/local/qcdb" CACHE PATH "Install path" FORCE)
endif()
message(STATUS "qcdb install: ${CMAKE_INSTALL_PREFIX}")

set(Python_ADDITIONAL_VERSIONS 3.7 3.6 3.5)  # adjust with CMake minimum FindPythonInterp
find_package(PythonLibsNew 2.7 REQUIRED)
message(STATUS "${Cyan}Found Python ${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR}${ColourReset}: ${PYTHON_EXECUTABLE} (found version ${PYTHON_VERSION_STRING})")

################################  Main Project  ################################
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

add_subdirectory(docs)

# <<<  Version  >>>
add_custom_target(update_version ALL
                  COMMAND ${PYTHON_EXECUTABLE} qcdb/versioner.py 
                                               --metaout ${CMAKE_CURRENT_BINARY_DIR}/metadata.py
                                               --cmakeout ${CMAKE_CURRENT_BINARY_DIR}/metadata.cmake
                  COMMAND cmake -DWTO="${CMAKE_CURRENT_BINARY_DIR}/${CMAKECONFIG_INSTALL_DIR}"
                                -DPN="qcdb"
                                -P ${CMAKE_CURRENT_BINARY_DIR}/metadata.cmake
                  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                  COMMENT "Generating version info")
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/metadata.py
        DESTINATION ${CMAKE_INSTALL_LIBDIR}${PYMOD_INSTALL_LIBDIR}/qcdb)

# <<<  Install  >>>

set(PN ${PROJECT_NAME})

install(DIRECTORY qcdb
        DESTINATION ${CMAKE_INSTALL_LIBDIR}${PYMOD_INSTALL_LIBDIR}
        USE_SOURCE_PERMISSIONS
        FILES_MATCHING PATTERN "*.py")

install(DIRECTORY share/qcdb/basis
        DESTINATION ${CMAKE_INSTALL_DATADIR}/${PN}
        USE_SOURCE_PERMISSIONS
        FILES_MATCHING PATTERN "*.gbs")

install(DIRECTORY tests
        DESTINATION ${CMAKE_INSTALL_LIBDIR}${PYMOD_INSTALL_LIBDIR}/qcdb
        USE_SOURCE_PERMISSIONS
        FILES_MATCHING PATTERN "*.py")

configure_file(qcdb/__init__.py qcdb/__init__.py @ONLY)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/qcdb/__init__.py
        DESTINATION ${CMAKE_INSTALL_LIBDIR}${PYMOD_INSTALL_LIBDIR}/qcdb)

# <<<  Export Config  >>>

# GNUInstallDirs "DATADIR" wrong here; CMake search path wants "share".
set(CMAKECONFIG_INSTALL_DIR "share/cmake/${PN}")
configure_package_config_file(cmake/${PN}Config.cmake.in
                              "${CMAKE_CURRENT_BINARY_DIR}/${PN}Config.cmake"
                              INSTALL_DESTINATION ${CMAKECONFIG_INSTALL_DIR})

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${PN}Config.cmake
              ${CMAKE_CURRENT_BINARY_DIR}/${PN}ConfigVersion.cmake
        DESTINATION ${CMAKECONFIG_INSTALL_DIR})

