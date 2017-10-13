set(shape_srcs ${TAUDEM_SRC}/ReadOutlets.cpp)
set(common_srcs ${TAUDEM_SRC}/commonLib.cpp ${TAUDEM_SRC}/tiffIO.cpp)

# Added by Liangjun Zhu, 2015/5/21
set(CURVATURE Curvature.cpp Curvaturemn.cpp ${common_srcs})

add_executable(curvature ${CURVATURE})

set(EXT_TARGETS curvature
        )
# third-party programs based on TauDEM framework
foreach (c_target ${EXT_TARGETS})
    target_link_libraries(${c_target} ${MPI_LIBRARIES} ${GDAL_LIBRARIES})
    install(TARGETS ${c_target} DESTINATION ${INSTALL_DIR})
    IF (MSVC)
        SET_PROPERTY(TARGET ${c_target} PROPERTY FOLDER "TauDEM_EXT")
    ENDIF ()
endforeach (c_target ${EXT_TARGETS})


