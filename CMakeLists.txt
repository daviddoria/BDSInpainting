cmake_minimum_required(VERSION 2.6)

PROJECT(BDSInpainting)
SET(CMAKE_INCLUDE_CURRENT_DIR ON)

# Enable C++11
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -std=gnu++0x")

# ITK
if(NOT ITK_FOUND)
  FIND_PACKAGE(ITK REQUIRED ITKCommon ITKIOImageBase ITKDistanceMap ITKIOPNG ITKIOMeta
                            ITKImageIntensity ITKImageFeature ITKMathematicalMorphology
                            ITKBinaryMathematicalMorphology)
  INCLUDE(${USE_ITK_FILE})
endif()

# VTK
if(NOT VTK_FOUND)
  FIND_PACKAGE(VTK REQUIRED)
  INCLUDE(${USE_VTK_FILE})
endif()

# Submodules
if(NOT TARGET PatchMatch)
  add_subdirectory(PatchMatch)
  include_directories(${PatchMatch_includes})
  set(BDSInpainting_includes ${BDSInpainting_includes} ${PatchMatch_includes})
  set(BDSInpainting_libraries ${BDSInpainting_libraries} ${PatchMatch_libraries})
endif()

if(NOT TARGET PoissonEditing)
  add_subdirectory(PoissonEditing)
  include_directories(${PoissonEditing_includes})
  set(BDSInpainting_includes ${BDSInpainting_includes} ${PoissonEditing_includes})
  set(BDSInpainting_libraries ${BDSInpainting_libraries} ${PoissonEditing_libraries})
endif()

if(NOT TARGET ITKVTKHelpers)
  add_subdirectory(ITKVTKHelpers)
  include_directories(${ITKVTKHelpers_includes})
  set(BDSInpainting_includes ${BDSInpainting_includes} ${ITKVTKHelpers_includes})
  set(BDSInpainting_libraries ${BDSInpainting_libraries} ${ITKVTKHelpers_libraries})
endif()

# Allow super-projects to use all of the includes we can use here.
if(NOT ${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_CURRENT_SOURCE_DIR})
  set(BDSInpainting_includes ${BDSInpainting_includes}
      ${CMAKE_CURRENT_SOURCE_DIR} PARENT_SCOPE)
  set(BDSInpainting_libraries ${BDSInpainting_libraries} PARENT_SCOPE)
endif()

# Create a library
add_custom_target(BDSInpainting BDSInpainting.hpp)

SET(BDSInpainting_BuildDrivers ON CACHE BOOL "Build BDSInpainting drivers?")
if(BDSInpainting_BuildDrivers)
 add_subdirectory(Drivers)
endif()

SET(BDSInpainting_BuildTests ON CACHE BOOL "BDSInpainting build tests?")
if(BDSInpainting_BuildTests)
 add_subdirectory(Tests)
endif()

GET_DIRECTORY_PROPERTY(AllBDSInpaintingIncludes INCLUDE_DIRECTORIES)
message("BDSInpainting include directories: ${AllBDSInpaintingIncludes}")
