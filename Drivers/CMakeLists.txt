# Allow headers in tests to be included like
# #include "PatchMatch.h" rather than needing
# #include "PatchMatch/PatchMatch.h"
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/..)

ADD_EXECUTABLE(BDSInpaintingDemo BDSInpaintingDemo.cpp)
TARGET_LINK_LIBRARIES(BDSInpaintingDemo ${PoissonEditingLibs} ${PatchMatchLibs})

#ADD_EXECUTABLE(BDSInpaintingRings BDSInpaintingRings.cpp)
#TARGET_LINK_LIBRARIES(BDSInpaintingRings ${BDSInpainting_libraries} ${PatchMatchLibs})

#message("BDSInpainting_libraries in BDSInpainting/Drivers: ${BDSInpainting_libraries}")

# ADD_EXECUTABLE(BDSInpaintingMultiRes BDSInpaintingMultiRes.cpp)
# TARGET_LINK_LIBRARIES(BDSInpaintingMultiRes ${PoissonEditingLibs} ${PatchMatchLibs})
