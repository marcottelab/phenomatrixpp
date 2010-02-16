cmake_minimum_required (2.6)
project (Phenomatrixpp)

set (Phenomatrixpp_PACKAGE_STRING "Phenomatrix++ 0.1")
set (Phenomatrixpp_PACKAGE_BUGREPORT "john.woods@marcottelab.org")


# Does this system provide the sqrt function?
include (${CMAKE_ROOT}/Modules/CheckFunctionExists.cmake)
check_function_exists (sqrt HAVE_SQRT)

# Include a config file
configure_file (
  "${PROJECT_SOURCE_DIR}/config.h.in"
  "${PROJECT_BINARY_DIR}/config.h"
  )
  
include_directories("${PROJECT_BINARY_DIR}")

find_package(Boost COMPONENTS system filesystem program_options REQUIRED)
target_link_libraries(Phenomatrixpp
  ${Boost_FILESYSTEM_LIBRARY}
  ${Boost_SYSTEM_LIBRARY}
  ${Boost_PROGRAM_OPTIONS_LIBRARY}
)

add_executable(phenomatrixpp main.cpp)
