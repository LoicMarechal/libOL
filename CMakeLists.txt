
################################
# GLOBAL CONFIGURATION VARIABLES
################################

cmake_minimum_required (VERSION 2.6)
project (libOL)
set (libOL_VERSION_MAJOR 1)
set (libOL_VERSION_MINOR 40)
set (CMAKE_C_FLAGS "-O3")
set (CMAKE_Fortran_FLAGS "-O3")
set (CMAKE_INSTALL_PREFIX "/opt/libOL")
find_library (MESHB_LIB Meshb.7 PATHS /opt/libMeshb/lib)
find_library (MESHBF_LIB Meshbf.7 PATHS /opt/libMeshb/lib)
find_path (MESHB_INC libmeshb7.h PATHS /opt/libMeshb/include)


#######################################
# SET FILES AND DIRECTORIES TO BE BUILT
#######################################

include_directories ("${PROJECT_SOURCE_DIR}/sources")
include_directories ("${MESHB_INC}")
add_subdirectory (sources)
add_subdirectory (examples)
install (FILES LICENSE copyright.txt DESTINATION share/libOL)
install (DIRECTORY sample_meshes DESTINATION share/libOL)


######################################
# SET PACKAGE AND DEPLOYMENT VARIABLES
######################################

include (InstallRequiredSystemLibraries)
set (CPACK_GENERATOR "TGZ")
set (CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/copyright.txt")
set (CPACK_PACKAGE_VERSION_MAJOR "${libOL_VERSION_MAJOR}")
set (CPACK_PACKAGE_VERSION_MINOR "${libOL_VERSION_MINOR}")
set(CPACK_COMPONENTS_ALL applications examples libraries headers)
include (CPack)
