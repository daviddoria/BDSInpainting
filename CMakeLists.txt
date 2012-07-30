cmake_minimum_required(VERSION 2.6)

PROJECT(BDSInpainting)
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -std=gnu++0x")

# Eigen
#include_directories(/media/portable/src/Eigen)

# ITK
FIND_PACKAGE(ITK REQUIRED ITKCommon ITKIOImageBase ITKDistanceMap ITKIOPNG ITKIOMeta
ITKImageIntensity ITKImageFeature ITKMathematicalMorphology ITKBinaryMathematicalMorphology)
INCLUDE(${USE_ITK_FILE})

# Submodules
if(NOT TARGET PatchMatch)
  add_subdirectory(PatchMatch)
  include_directories(${PatchMatch_includes})
endif()

if(NOT TARGET PoissonEditing)
  add_subdirectory(PoissonEditing)
  include_directories(${PoissonEditing_includes})
endif()

# Allow super-projects to use all of the includes we can use here.
if(NOT ${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_CURRENT_SOURCE_DIR})
  set(BDSInpainting_includes ${PatchMatch_includes}
      ${CMAKE_CURRENT_SOURCE_DIR} PARENT_SCOPE)
endif()

# Create a library
add_custom_target(BDSInpainting BDSInpainting.hpp)
#target_link_libraries(BDSInpainting Mask PatchComparison)

SET(BDSInpainting_BuildDrivers ON CACHE BOOL "Build BDSInpainting drivers?")
if(BDSInpainting_BuildDrivers)
 add_subdirectory(Drivers)
endif()
