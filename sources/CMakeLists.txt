
##########################
# BUILD THE STATIC LIBRARY
##########################

add_library(OL.1 libol1.c)
install (FILES libol1.h libol1.ins DESTINATION include COMPONENT headers)
install (TARGETS OL.1 EXPORT libOL-target DESTINATION lib COMPONENT libraries)
install (EXPORT libOL-target DESTINATION lib/cmake/${PROJECT_NAME})
export  (PACKAGE libOL) 
