
##################
# BUILD C EXAMPLES
##################

add_executable(filtered_projections filtered_projections.c)
target_link_libraries(filtered_projections OL.1 ${libMeshb_LIBRARIES} ${math_LIBRARIES})
install (TARGETS filtered_projections DESTINATION share/libOL/examples COMPONENT examples)

add_executable(test_libol test_libol1.c)
target_link_libraries(test_libol OL.1 ${libMeshb_LIBRARIES} ${math_LIBRARIES})
install (TARGETS test_libol DESTINATION share/libOL/examples COMPONENT examples)

if(LPlib_FOUND)
   add_executable(parallel_libol parallel_libol.c)
   target_link_libraries(parallel_libol OL.1 ${libMeshb_LIBRARIES} ${LPlib_LIBRARIES} ${math_LIBRARIES})
   install (TARGETS parallel_libol DESTINATION share/libOL/examples COMPONENT examples)
endif ()


########################
# BUILD FORTRAN EXAMPLES
########################

if (CMAKE_Fortran_COMPILER)
   add_executable(test_libol_f77 test_libol1.f)
   target_link_directories(test_libol_f77 PRIVATE ${libMeshb_LINK_DIRS})
   target_link_libraries(test_libol_f77 OL.1 gfortran ${libMeshb_Fortran_LIBRARIES})
   install (TARGETS test_libol_f77 DESTINATION share/libOL/examples COMPONENT examples)
endif ()
