cmake_minimum_required(VERSION 2.6)

# clear; ~/bin/cmake-3.2.2-Linux-x86_64/bin/cmake ~/src/BDSInpainting/ -DBoost_DEBUG=ON -DITK_DIR=/home/doriad/build/ITK_devel/ -DVTK_DIR=/home/doriad/build/VTK -DEIGEN3_INCLUDE_DIR=/home/doriad/src/eigen-3.2.1 -DBOOST_ROOT=/home/doriad/src/boost_1_55_0/
PROJECT(BDSInpainting)
SET(CMAKE_INCLUDE_CURRENT_DIR ON)

add_subdirectory(CMakeHelpers)

# Enable C++11
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -std=gnu++11")

# ITK
if(NOT ITK_FOUND)
  FIND_PACKAGE(ITK REQUIRED ITKCommon ITKIOImageBase ITKDistanceMap ITKIOPNG ITKIOMeta
                            ITKImageIntensity ITKImageFeature ITKMathematicalMorphology
                            ITKBinaryMathematicalMorphology)
  INCLUDE(${ITK_USE_FILE})
endif()

# VTK
if(NOT VTK_FOUND)
  FIND_PACKAGE(VTK REQUIRED)
  INCLUDE(${VTK_USE_FILE})
endif()

#### Eigen ####
# Tell CMake to also look in the source directory to find some .cmake files
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_MODULE_PATH})

# Set the CMake variable to the value of the environment variable.
# If the build is absolutely fresh, this will be the case if the CMake variable EIGEN3_INCLUDE_DIR is not set.
if(NOT EIGEN3_INCLUDE_DIR)
  set(EIGEN3_INCLUDE_DIR $ENV{EIGEN3_INCLUDE_DIR})
  message("Set EIGEN3_INCLUDE_DIR to ${EIGEN3_INCLUDE_DIR} from environment variable of the same name.")
else()
  message("EIGEN3_INCLUDE_DIR is ${EIGEN3_INCLUDE_DIR}")
endif()

if(NOT EIGEN3_FOUND)
  FIND_PACKAGE(Eigen3 REQUIRED) #requires FindEigen3.cmake to be in the source directory
  include_directories(${EIGEN3_INCLUDE_DIR})
endif()

# Boost
set(Boost_NO_BOOST_CMAKE TRUE) # seems to avoid "add_library cannot create imported target" errors

FIND_PACKAGE(Boost 1.51 REQUIRED) # The signals2 library is header only, so we don't need to / can't search for it in find_package

IF (Boost_FOUND)
    INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIR})
ENDIF()

# Submodules
UseSubmodule(PatchMatch BDSInpainting)
UseSubmodule(PoissonEditing BDSInpainting)
UseSubmodule(ITKVTKHelpers BDSInpainting)

# Add non-compiled files to the project
add_custom_target(BDSInpainting SOURCES
BDSInpainting.h
BDSInpainting.hpp
BDSInpaintingMultiRes.h
BDSInpaintingMultiRes.hpp
BDSInpaintingRings.h
BDSInpaintingRings.hpp
Compositor.h
Compositor.hpp
InpaintingAlgorithm.h
InpaintingAlgorithm.hpp
PixelCompositors.h)

SET(BDSInpainting_BuildDrivers ON CACHE BOOL "Build BDSInpainting drivers?")
if(BDSInpainting_BuildDrivers)
 add_subdirectory(Drivers)
endif()

#SET(BDSInpainting_BuildTests ON CACHE BOOL "BDSInpainting build tests?")
#if(BDSInpainting_BuildTests)
# add_subdirectory(Tests)
#endif()
