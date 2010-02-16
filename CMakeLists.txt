cmake_minimum_required (VERSION 2.6)
project (Phenomatrixpp)

set (Phenomatrixpp_PACKAGE_STRING "Phenomatrix++ 0.1")
set (Phenomatrixpp_PACKAGE_BUGREPORT "john.woods@marcottelab.org")
set (Phenomatrixpp_SOURCES main.cpp adjacency_list.cpp avgmindist.cpp distance.cpp euclidean.cpp genephene.cpp hypergeometric.cpp knearest.cpp marshall.cpp mindist.cpp oracle.cpp partialbayes.cpp type_shield.cpp utilities.cpp)
message("PROJECT_NAME=${PROJECT_NAME}")
message("PROJECT_SOURCE_DIR=${PROJECT_SOURCE_DIR}")
message("Phenomatrixpp_SOURCE_DIR=${Phenomatrixpp_SOURCE_DIR}")

# Does this system provide the sqrt function?
include (${CMAKE_ROOT}/Modules/CheckFunctionExists.cmake)
check_function_exists (sqrt HAVE_SQRT)

# Include a config file
configure_file (
  "${PROJECT_SOURCE_DIR}/config.h.in"
  "${PROJECT_BINARY_DIR}/config.h"
  )
  
include_directories("${PROJECT_BINARY_DIR}")

find_package(Boost 1.41.0 COMPONENTS system filesystem program_options REQUIRED)
MESSAGE("Boost information:") 
MESSAGE("  Boost_INCLUDE_DIRS: ${Boost_INCLUDE_DIRS}") 
MESSAGE("  Boost_LIBRARIES: ${Boost_LIBRARIES}") 
MESSAGE("  Boost_LIBRARY_DIRS: ${Boost_LIBRARY_DIRS}") 

add_executable(phenomatrixpp ${Phenomatrixpp_SOURCES})

target_link_libraries(phenomatrixpp
  ${Boost_FILESYSTEM_LIBRARY}
  ${Boost_SYSTEM_LIBRARY}
  ${Boost_PROGRAM_OPTIONS_LIBRARY}
)

