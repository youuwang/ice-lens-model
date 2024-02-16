cmake_minimum_required(VERSION 3.1)
project(ice_model)
enable_language(CXX)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# flags for compilation (options)
set(CMAKE_CXX_FLAGS "-Wall -Wextra -pedantic -std=c++17 " CACHE STRING "" FORCE)
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG " CACHE STRING "" FORCE )

# alternative to add_library (but with option)
set(BUILD_SHARED_LIBS ON CACHE BOOL "Build shared libraries.")
mark_as_advanced(BUILD_SHARED_LIBS)

# build type
set(CMAKE_BUILD_TYPE "${CMAKE_BUILD_TYPE}" CACHE STRING
  "Build options: None Debug Release RelWithDebInfo MinSizeRel."
  FORCE )

# Add the cmake folder so the FindSphinx module is found
set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

# include directories
set(ICE_MODEL_INCLUDE_DIRS
  "${CMAKE_CURRENT_BINARY_DIR}/src"
  "${CMAKE_CURRENT_SOURCE_DIR}/src"
  )
include_directories(${ICE_MODEL_INCLUDE_DIRS})

# source directory
add_subdirectory(src)

# examples: option to de/activate compilation of examples
option(ICE_MODEL_EXAMPLES "example" OFF)
if(ICE_MODEL_EXAMPLES)
  add_subdirectory(example)
endif()

#------------------------------------------------------------------------ #
# Doxygen                                                                 #
# ----------------------------------------------------------------------- #
# verbose:
option(ICE_MODEL_DOC "build documentation for the project" OFF)
if(ICE_MODEL_DOC)
  add_subdirectory(docs)  
endif(ICE_MODEL_DOC)
