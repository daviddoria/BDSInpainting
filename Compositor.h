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
template <typename TImage, typename TPixelCompositor>
class Compositor
{
public:

  /** Constructor. */
  Compositor();

  typedef itk::Image<Match, 2> NNFieldType;

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

protected:

  /** The radius of the patches to use for inpainting. */
  unsigned int PatchRadius;

  /** The output image. */
  typename TImage::Pointer Output;

  /** The image to fill. */
  typename TImage::Pointer Image;

  /** The mask where Hole pixels indicate the pixels to fill. */
  Mask::Pointer TargetMask;

  /** The nearest neighbor field to use. */
  NNFieldType* NearestNeighborField;
};

class PixelCompositorAverage
{
  /** Composite by averaging pixels. */
  static typename TImage::PixelType Composite(
    const std::vector<typename TImage::PixelType>& contributingPixels,
    const std::vector<float>& contributingScores)
  {
    typename TImage::PixelType newValue = Statistics::Average(contributingPixels);
    return newValue;
  }
  
};

class PixelCompositorWeightedAverage
{
  /** Composite by averaging pixels. */
  static typename TImage::PixelType Composite(
    const std::vector<typename TImage::PixelType>& contributingPixels,
    const std::vector<float>& contributingScores)
  {
    // If there is only one element, simply return it.
    if(contributingScores.size() == 1)
    {
      return contributingPixels[0];
    }

    // The weights should be inversely proportional to the patch errors/scores.
    // That is, a patch with a high error should get a low weight. We accomplish this by
    // making the weight equal to 1 - (value - min) / |range|

    float minValue = *std::min_element(contributingScores.begin(), contributingScores.end());
    float maxValue = *std::max_element(contributingScores.begin(), contributingScores.end());
    float range = maxValue - minValue;
    // If the range is zero, all elements are the same so just return the first one
    if(range == 0.0f)
    {
      return contributingPixels[0];
    }

    std::vector<float> weights(contributingScores.size());
    for(unsigned int i = 0; i < contributingScores.size(); ++i)
    {
      weights[i] = 1.0f - (contributingScores[i] - minValue) / range;
    }

    // Make the weights sum to 1
    Helpers::NormalizeVectorInPlace(weights);

    typename TImage::PixelType newValue = Helpers::WeightedAverage(contributingPixels, weights);
    return newValue;
  }
};

class PixelCompositorClosestToAverage
{
  /** Composite by averaging pixels. */
  static typename TImage::PixelType Composite(
    const std::vector<typename TImage::PixelType>& contributingPixels,
    const std::vector<float>& contributingScores)
  {
    // Use the pixel closest to the average pixel
    typename TImage::PixelType averagePixel = Statistics::Average(contributingPixels);
    unsigned int patchId = ITKHelpers::ClosestValueIndex(contributingPixels, averagePixel);
    typename TImage::PixelType newValue = contributingPixels[patchId];
    return newValue;
  }
};

class PixelCompositorBestPatch
{
  /** Composite by averaging pixels. */
  static typename TImage::PixelType Composite(
    const std::vector<typename TImage::PixelType>& contributingPixels,
    const std::vector<float>& contributingScores)
  {
    // Take the pixel from the best matching patch
    unsigned int patchId = Helpers::argmin(contributingScores);
    typename TImage::PixelType newValue = contributingPixels[patchId];
    return newValue;
  }
};

#include "Compositor.hpp"

#endif
