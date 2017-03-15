
################################
# GLOBAL CONFIGURATION VARIABLES
################################

cmake_minimum_required (VERSION 3.0.1)
set (libOL_VERSION_MAJOR 1)
set (libOL_VERSION_MINOR 40)

if (WIN32)
   set (CMAKE_C_FLAGS -DWIN32)
   set (CMAKE_CONFIGURATION_TYPES Release)
else ()
   set (CMAKE_C_FLAGS -O3)
   set (CMAKE_Fortran_FLAGS -O3)
   set (LINK_LIBRARIES m)
endif()

project (libOL VERSION ${libOL_VERSION_MAJOR}.${libOL_VERSION_MINOR} LANGUAGES C)

include (CheckLanguage)
check_language (Fortran)
if(CMAKE_Fortran_COMPILER)
  enable_language(Fortran)
endif()


#######################################
# SET FILES AND DIRECTORIES TO BE BUILT
#######################################

include_directories (${PROJECT_SOURCE_DIR}/sources)
add_subdirectory (sources)

find_package(libMeshb 7)

if(libMeshb_FOUND)
   include_directories (${libMeshb_INCLUDE_DIRS})
   include_directories (${PROJECT_SOURCE_DIR}/examples)
   add_subdirectory (examples)
else()
   message (WARNING "
   libMeshb.7 could not be found.
   Consequently, examples and utilities will not be built.
   You can download it from GitHUB at https://github.com/LoicMarechal/libMeshb")
endif()

install (FILES LICENSE copyright.txt DESTINATION share/libOL)
install (DIRECTORY sample_meshes DESTINATION share/libOL)


######################################
# SET PACKAGE AND DEPLOYMENT VARIABLES
######################################

include (InstallRequiredSystemLibraries)
set (CPACK_GENERATOR TXZ)
set (CPACK_RESOURCE_FILE_LICENSE ${CMAKE_CURRENT_SOURCE_DIR}/copyright.txt)
set (CPACK_PACKAGE_VERSION_MAJOR ${libOL_VERSION_MAJOR})
set (CPACK_PACKAGE_VERSION_MINOR ${libOL_VERSION_MINOR})
set(CPACK_COMPONENTS_ALL applications examples libraries headers)
include (CPack)


######################################
# EXPORT CMAKE PACKAGE INFORMATION
######################################

include(CMakePackageConfigHelpers)
file(WRITE ${PROJECT_BINARY_DIR}/${PROJECT_NAME}Config.cmake
"include(\${CMAKE_CURRENT_LIST_DIR}/libOL-target.cmake)
set(libOL_INCLUDE_DIRS ${CMAKE_INSTALL_PREFIX}/include)
set(libOL_LIBRARIES OL.1)
set(libOL_FOUND TRUE)
")
write_basic_package_version_file(
   "${PROJECT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake"
   VERSION ${PROJECT_VERSION}
   COMPATIBILITY SameMajorVersion)
install(FILES
   "${PROJECT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
   "${PROJECT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake"
   DESTINATION lib/cmake/${PROJECT_NAME})
