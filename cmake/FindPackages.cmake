### GDAL (REQUIRED)
FIND_PACKAGE(GDAL REQUIRED)
IF(NOT GDAL_FOUND)
  MESSAGE(FATAL_ERROR "FATAL: Could not find GDAL!")
ENDIF()

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
        GET_FILENAME_COMPONENT(GDAL_LIBRARY_PATH ${GDAL_LIBRARIES} DIRECTORY)
        FILE(GLOB GDAL_LIBWITH_VER "${GDAL_LIBRARY_PATH}/${GDAL_LIBRARY_NAME}.1.*.*")
        GET_FILENAME_COMPONENT(GDAL_VERSION ${GDAL_LIBWITH_VER} EXT)
        STRING(REPLACE ".dll.1" "1" GDAL_VERSION ${GDAL_VERSION})
        IF(GDAL_LIBRARY_NAME STREQUAL "libgdal.dll")
            SET(OPTIMIZED_GDAL_LIBRARY ${GDAL_LIBRARIES}.${GDAL_VERSION})
            SET(DEBUG_GDAL_LIBRARY ${GDAL_LIBRARIES}.${GDAL_VERSION})
            STRING(REPLACE gdal gdald DEBUG_GDAL_LIBRARY ${DEBUG_GDAL_LIBRARY})
            SET(GDAL_LIBRARIES
                optimized ${OPTIMIZED_GDAL_LIBRARY}
                debug ${DEBUG_GDAL_LIBRARY})
        ENDIF()
    ENDIF()
ENDIF()

### MPI (REQUIRED)
# currently, for mingw64 on Windows, you must specify MPI_LIBRARIES and DMPI_INCLUDE_PATH through CMake command.
IF(NOT DEFINED MPI_LIBRARIES)
  FIND_PACKAGE(MPI REQUIRED)
  SET(CMAKE_CXX_FLAG "${CMAKE_CXX_FLAG} ${MPI_COMPILE_FLAGS}")
  SET(CMAKE_CXX_LINK_FLAGS "${CMAKE_CXX_LINK_FLAGS} ${MPI_LINK_FLAGS}")
  IF (CMAKE_COMPILER_IS_CLANGCXX OR CMAKE_COMPILER_IS_CLANGCC) # using regular Clang or AppleClang
    SET(CMAKE_EXE_LINKER_FLAGS "-Wl -headerpad_max_install_names")
  ELSEIF (GEO_ICC)
    # Select flags for ICC, refers to https://software.intel.com/en-us/articles/intel-cluster-toolkit-for-linux-error-when-compiling-c-aps-using-intel-mpi-library-compilation-driver-mpiicpc
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DMPICH_IGNORE_CXX_SEEK -DMPICH_SKIP_MPICXX")
  ENDIF ()
ENDIF()
