cmake_minimum_required(VERSION 2.6)

PROJECT(BDSInpainting)
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -std=gnu++0x")

# Eigen
include_directories(/media/portable/src/Eigen)

# VTK
FIND_PACKAGE(ITK REQUIRED ITKCommon ITKIOImageBase
ITKImageIntensity ITKImageFeature ITKMathematicalMorphology ITKBinaryMathematicalMorphology
ITKDeprecated)
INCLUDE(${USE_ITK_FILE})

add_subdirectory(PatchMatch)

add_library(BDSInpainting BDSInpainting.cpp)
target_link_libraries(BDSInpainting PatchMatch)

ADD_EXECUTABLE(BDSInpaintingDemo BDSInpaintingDemo.cpp)
TARGET_LINK_LIBRARIES(BDSInpaintingDemo BDSInpainting)
