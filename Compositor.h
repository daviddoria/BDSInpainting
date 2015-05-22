/*=========================================================================
 *
 *  Copyright David Doria 2012 daviddoria@gmail.com
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/

#ifndef Compositor_H
#define Compositor_H

// ITK
#include "itkCovariantVector.h"
#include "itkImage.h"
#include "itkImageRegion.h"
#include "itkVectorImage.h"

// Submodules
#include <Mask/Mask.h>

#include <PatchMatch/PatchMatchHelpers.h>
#include <PatchMatch/NNField.h>

class CompositorParent
{
  virtual void Composite() = 0;
};

/** This class takes a nearest neighbor field and a target mask and
  * fills the valid pixels in the target mask using one of the specified strategies. */
template <typename TImage, typename TPixelCompositor>
class Compositor : public CompositorParent
{
public:

  /** Constructor. */
  Compositor();

  /** Set the nearest neighbor field to use. */
  void SetNearestNeighborField(NNFieldType* const nnField);

  /** Get the resulting inpainted image. */
  TImage* GetOutput();

  /** Set the patch radius. */
  void SetPatchRadius(const unsigned int patchRadius);

  /** Set the image to fill. */
  void SetImage(TImage* const image);

  /** Set the mask that indicates where to fill the image. Pixels in the Hole region should be filled.*/
  void SetTargetMask(Mask* const mask);

  /** Perform the compositing where the TargetMask is valid.*/
  void Composite();

protected:

  /** The radius of the patches to use for inpainting. */
  unsigned int PatchRadius = 0;

  /** The output image. */
  typename TImage::Pointer Output = nullptr;

  /** The image to fill. */
  typename TImage::Pointer Image = nullptr;

  /** The mask where Hole pixels indicate the pixels to fill. */
  Mask::Pointer TargetMask = nullptr;

  /** The nearest neighbor field to use. */
  NNFieldType* NearestNeighborField = nullptr;
};

#include "Compositor.hpp"

#endif
