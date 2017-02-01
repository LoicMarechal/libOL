
##########################
# BUILD THE STATIC LIBRARY
##########################

add_library(OL.1 libol1.c)
install (TARGETS OL.1 DESTINATION lib COMPONENT libraries)
install (FILES libol1.h libol1.ins DESTINATION include COMPONENT headers)
