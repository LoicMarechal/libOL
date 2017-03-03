
##################
# BUILD C EXAMPLES
##################

add_executable(test_libol test_libol1.c)
target_link_libraries(test_libol OL.1 ${libMeshb_C_LIBRARY} ${LINK_LIBRARIES})
install (TARGETS test_libol DESTINATION share/libOL/examples COMPONENT examples)


########################
# BUILD FORTRAN EXAMPLES
########################

if (CMAKE_Fortran_COMPILER)
   add_executable(test_libol_f77 test_libol1.f)
   target_link_libraries(test_libol_f77 OL.1 ${libMeshb_Fortran_LIBRARY})
   install (TARGETS test_libol_f77 DESTINATION share/libOL/examples COMPONENT examples)
endif ()
