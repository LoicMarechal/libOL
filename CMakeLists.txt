
################################
# GLOBAL CONFIGURATION VARIABLES
################################

cmake_minimum_required (VERSION 3.7.2)
set (libOL_VERSION_MAJOR 1)
set (libOL_VERSION_MINOR 60)
project (libOL VERSION ${libOL_VERSION_MAJOR}.${libOL_VERSION_MINOR} LANGUAGES C)

include (CheckLanguage)
check_language (Fortran)
if(CMAKE_Fortran_COMPILER)
  enable_language(Fortran)
  if (CMAKE_Fortran_COMPILER_ID STREQUAL GNU)
     set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -std=legacy")
  endif ()
endif()

if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif ()

if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
   if (WIN32)
      set(CMAKE_INSTALL_PREFIX "$ENV{HOMEPATH}/local" CACHE PATH "..." FORCE)
   else ()
      set(CMAKE_INSTALL_PREFIX "$ENV{HOME}/local" CACHE PATH "..." FORCE)
   endif()
endif ()


#######################################
# SET FILES AND DIRECTORIES TO BE BUILT
#######################################

include_directories (${PROJECT_SOURCE_DIR}/sources)
add_subdirectory (sources)

find_package(LPlib 3)
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

install (FILES LICENSE.txt copyright.txt DESTINATION share/libOL)
install (DIRECTORY sample_meshes DESTINATION share/libOL)


######################################
# SET PACKAGE AND DEPLOYMENT VARIABLES
######################################

if (CPACK)
   include (InstallRequiredSystemLibraries)
   set (CPACK_GENERATOR TXZ)
   set (CPACK_RESOURCE_FILE_LICENSE ${CMAKE_CURRENT_SOURCE_DIR}/copyright.txt)
   set (CPACK_PACKAGE_VERSION_MAJOR ${libOL_VERSION_MAJOR})
   set (CPACK_PACKAGE_VERSION_MINOR ${libOL_VERSION_MINOR})
   set(CPACK_COMPONENTS_ALL applications examples libraries headers)
   include (CPack)
endif ()


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
