# CMake build for Windows, Linux, and macOS
# Updated by Liang-Jun Zhu
# Description: Separated the CMakeLists.txt into two
# zlj@lreis.ac.cn
cmake_minimum_required(VERSION 2.8)
project(TauDEM)
IF (COMMAND CMAKE_POLICY)
    CMAKE_POLICY(SET CMP0015 NEW)
ENDIF (COMMAND CMAKE_POLICY)
################# Install directory ##################
IF (INSTALL_PREFIX)
    SET(CMAKE_INSTALL_PREFIX ${INSTALL_PREFIX})
    SET(INSTALL_DIR ${INSTALL_PREFIX})
ELSE ()
    SET(CMAKE_INSTALL_PREFIX ${CMAKE_SOURCE_DIR})
    SET(INSTALL_DIR ${CMAKE_SOURCE_DIR}/bin)
ENDIF ()
################# Find packages ##################
# MPI and GDAL are required
IF (MSVC)
    # 1. GDAL
    SET(GDAL_DIR $ENV{GDAL_DIR})
    SET(GDAL_DIR $ENV{GDAL_DIR})
    SET(GDAL_INCLUDE_DIR ${GDAL_DIR}\\include)
    SET(GDAL_LIBRARIES ${GDAL_DIR}\\lib)
    SET(GDAL_BIN_DIR ${GDAL_DIR}\\bin)
    SET(GDAL_LIBRARY ${GDAL_LIBRARIES}\\gdal_i.lib)
    # 2. MPI
    IF (ARCH EQUAL 64)
        SET(MPI_LIB $ENV{MSMPI_LIB64})
    ELSE ()
        SET(MPI_LIB $ENV{MSMPI_LIB32})
    ENDIF ()
    SET(MPI_HOME $ENV{MSMPI_BIN})
    SET(MPI_INCLUDE_PATH $ENV{MSMPI_INC})
    SET(MPI_LIBRARIES ${MPI_LIB}\\msmpi.lib)

    MESSAGE(STATUS "Found GDAL: ${GDAL_LIBRARIES}")
    MESSAGE(STATUS "Found MPI: ${MPI_LIBRARIES}")
ELSE ()
    # for Linux and macOS
    # 1. GDAL, CMake will find GDAL and get GDAL_INCLUDE_DIR and GDAL_LIBRARY.
    find_package(GDAL REQUIRED)
    # 2. MPI
    find_package(MPI REQUIRED)
    set(CMAKE_CXX_FLAG "${CMAKE_CXX_FLAG} ${MPI_COMPILE_FLAGS}")
    set(CMAKE_CXX_LINK_FLAGS "${CMAKE_CXX_LINK_FLAGS} ${MPI_LINK_FLAGS}")
    IF (APPLE)
        set(CMAKE_EXE_LINKER_FLAGS "-Wl -headerpad_max_install_names")
    ENDIF ()
ENDIF ()
# Fixup GDAL variable.
IF(WIN32 AND MINGW)
    # On Windows, find_package always selects the release dll named
    # libgdal.dll. This is a symbolic links which trips the linker. Also,
    # for debug builds, we want to link against the debug library.
    # Before fixing this, verify that the variable contains a single string
    # and that the name is libgdal.dll. The CMake module may be fixed at
    # some point.
    LIST(LENGTH GDAL_LIBRARIES NR_GDAL_LIBRARIES)
    IF(NR_GDAL_LIBRARIES EQUAL 1)
        GET_FILENAME_COMPONENT(GDAL_LIBRARY_NAME ${GDAL_LIBRARIES} NAME)

        IF(GDAL_LIBRARY_NAME STREQUAL "libgdal.dll")
            SET(GDAL_VERSION "1.11.5")
            SET(OPTIMIZED_GDAL_LIBRARY ${GDAL_LIBRARIES}.${GDAL_VERSION})
            SET(DEBUG_GDAL_LIBRARY ${GDAL_LIBRARIES}.${GDAL_VERSION})
            STRING(REPLACE gdal gdald DEBUG_GDAL_LIBRARY ${DEBUG_GDAL_LIBRARY})
            SET(GDAL_LIBRARIES
                    optimized ${OPTIMIZED_GDAL_LIBRARY}
                    debug ${DEBUG_GDAL_LIBRARY})
        ENDIF()
    ENDIF()
ENDIF()

MESSAGE(STATUS "GDAL include: ${GDAL_INCLUDE_DIR}")
MESSAGE(STATUS "MPI include: ${MPI_INCLUDE_PATH}")

# Set a default build type for single-configuration
# CMake generators if no build type is set.
IF(NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
   SET(CMAKE_BUILD_TYPE RelWithDebInfo)
ENDIF()

# Select flags.
# refers to https://software.intel.com/en-us/articles/intel-cluster-toolkit-for-linux-error-when-compiling-c-aps-using-intel-mpi-library-compilation-driver-mpiicpc
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DMPICH_IGNORE_CXX_SEEK -DMPICH_SKIP_MPICXX")
IF(MSVC)
    SET(CMAKE_CXX_FLAGS_RELWITHDEBINFO "/O2")
    SET(CMAKE_CXX_FLAGS_RELEASE "/O2 /Zi")
    SET(CMAKE_CXX_FLAGS_DEBUG  "/O2 /Zi")
ELSE()
    SET(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -g")
    SET(CMAKE_CXX_FLAGS_RELEASE "-O2")
    SET(CMAKE_CXX_FLAGS_DEBUG  "-O0 -g")
ENDIF()
################# Packages found ##################
IF(MSVC)
  SET_PROPERTY(GLOBAL PROPERTY USE_FOLDERS ON)
ENDIF()
SET(TAUDEM_SRC ${CMAKE_SOURCE_DIR}/src)
SET(TAUDEM_EXT_SRC ${CMAKE_SOURCE_DIR}/src_ext)

# GDAL_INCLUDE_DIR should after TAUDEM_SRC, since tiffIO.h both existed in the two dirs.
INCLUDE_DIRECTORIES(${TAUDEM_SRC} ${GDAL_INCLUDE_DIR} ${MPI_INCLUDE_PATH})
ADD_SUBDIRECTORY(${TAUDEM_SRC})
ADD_SUBDIRECTORY(${TAUDEM_EXT_SRC})
### For CLion to implement the "make install" command
add_custom_target(install_${PROJECT_NAME}
        $(MAKE) install
        DEPENDS ${EXT_TARGETS} ${MY_TARGETS}
        COMMENT "Installing ${PROJECT_NAME}")
MESSAGE(STATUS "    Compiling TauDEM done")
