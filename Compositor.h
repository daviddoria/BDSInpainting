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

/** This class takes a nearest neighbor field and a target mask and
  * fills the valid pixels in the target mask using one of the specified strategies. */
template <typename TImage>
class Compositor
{
public:

  /** Constructor. */
  Compositor();

  /** Choices of compositing method. */
  enum CompositingMethodEnum {BEST_PATCH, WEIGHTED_AVERAGE, AVERAGE, CLOSEST_TO_AVERAGE};

  /** Set the compositing method to use. */
  void SetCompositingMethod(const CompositingMethodEnum& compositingMethod);

  typedef itk::Image<Match, 2> NNFieldType;

  /** Set the nearest neighbor field to use. */
  void SetNearestNeighborField(NNFieldType* const nnField);

  /** Using the provided NearestNeighborField, compute new values for the hole pixels
    * in the TargetMask. */
  void Compute();

  /** Get the resulting inpainted image. */
  TImage* GetOutput();

  /** Set the patch radius. */
  void SetPatchRadius(const unsigned int patchRadius);

  /** Set the image to fill. */
  void SetImage(TImage* const image);

  /** Set the mask that indicates where to fill the image. Pixels in the Hole region should be filled.*/
  void SetTargetMask(Mask* const mask);

protected:

  /** The radius of the patches to use for inpainting. */
  unsigned int PatchRadius;

  /** The output image. */
  typename TImage::Pointer Output;

  /** The image to fill. */
  typename TImage::Pointer Image;

  /** The mask where Hole pixels indicate the pixels to fill. */
  Mask::Pointer TargetMask;

  /** The compositing method to use. */
  CompositingMethodEnum CompositingMethod;

  /** Composite using the specified CompositingMethod. */
  typename TImage::PixelType Composite(
     const std::vector<typename TImage::PixelType>& contributingPixels,
     const std::vector<float>& contributingScores);

  /** Composite by averaging pixels. */
  typename TImage::PixelType CompositeAverage(
    const std::vector<typename TImage::PixelType>& contributingPixels);

  /** Composite by weighted averaging. */
  typename TImage::PixelType CompositeWeightedAverage(
     const std::vector<typename TImage::PixelType>& contributingPixels,
     const std::vector<float>& contributingScores);

  /** Composite by choosing the pixel closest to the average. */
  typename TImage::PixelType CompositeClosestToAverage(
     const std::vector<typename TImage::PixelType>& contributingPixels,
     const std::vector<float>& contributingScores);

  /** Composite by choosing the pixel from the best patch. */
  typename TImage::PixelType CompositeBestPatch(
     const std::vector<typename TImage::PixelType>& contributingPixels,
     const std::vector<float>& contributingScores);

  /** The nearest neighbor field to use. */
  NNFieldType* NearestNeighborField;
};

#include "Compositor.hpp"

#endif
