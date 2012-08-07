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

#ifndef Compositor_HPP
#define Compositor_HPP

#include "Compositor.h"

// Submodules
#include <ITKHelpers/ITKHelpers.h>
#include <Mask/MaskOperations.h>

// ITK
#include "itkImageRegionReverseIterator.h"

// STL
#include <ctime>

template <typename TImage>
Compositor<TImage>::Compositor() : PatchRadius(7),
   CompositingMethod(WEIGHTED_AVERAGE)
{
  this->Output = TImage::New();
  this->Image = TImage::New();
  this->TargetMask = Mask::New();
}

template <typename TImage>
TImage* Compositor<TImage>::GetOutput()
{
  return this->Output;
}

template <typename TImage>
void Compositor<TImage>::SetPatchRadius(const unsigned int patchRadius)
{
  this->PatchRadius = patchRadius;
}

template <typename TImage>
void Compositor<TImage>::SetImage(TImage* const image)
{
  ITKHelpers::DeepCopy(image, this->Image.GetPointer());
}

template <typename TImage>
void Compositor<TImage>::SetTargetMask(Mask* const mask)
{
  this->TargetMask->DeepCopyFrom(mask);
}

template <typename TImage>
void Compositor<TImage>::Compute()
{
  // The contribution of each pixel q to the error term (d_cohere) = 1/N_T \sum_{i=1}^m (S(p_i) - T(q))^2
  // To find the best color T(q) (iterative update rule), differentiate with respect to T(q),
  // set to 0, and solve for T(q):
  // T(q) = \frac{1}{m} \sum_{i=1}^m S(p_i)

//   std::cout << "Compositor::Compute()..." << std::endl;
//   ITKHelpers::WriteRGBImage(oldImage, "Compositor_Compute_OldImage.png");
//   ITKHelpers::WriteImage(targetMask, "Compositor_Compute_TargetMask.png");

  // This is done so in the algorithm we can use 'fullRegion', since it refers
  // to the same region for the image and mask.
  // (So there is no confusion such as "why is the mask's region used here instead of the image's?")
  itk::ImageRegion<2> fullRegion = this->Image->GetLargestPossibleRegion();

  // We don't want to change pixels directly on the output image during the iteration,
  // but rather compute them all and then update them all simultaneously.
  typename TImage::Pointer updatedImage = TImage::New();
  ITKHelpers::DeepCopy(this->Image.GetPointer(), updatedImage.GetPointer());

  // We must get a dummy pixel from the image and then fill it with zero to make sure the number
  // of components of the pixel is correct.
  itk::Index<2> zeroIndex = {{0,0}};
  typename TImage::PixelType zeroPixel = this->Image->GetPixel(zeroIndex);
  zeroPixel.Fill(0);

  std::vector<itk::Index<2> > targetPixels = this->TargetMask->GetValidPixels();
  std::cout << "Compositor::Compute(): There are : "
            << targetPixels.size() << " target pixels." << std::endl;

  for(size_t targetPixelId = 0; targetPixelId < targetPixels.size(); ++targetPixelId)
  {
    //std::cout << "Processing " << targetPixelId << " of " << targetPixels.size() << std::endl;
    itk::Index<2> currentPixel = targetPixels[targetPixelId];

    { // debug only
    itk::ImageRegion<2> currentRegion =
          ITKHelpers::GetRegionInRadiusAroundPixel(currentPixel, this->PatchRadius);

    //ITKHelpers::WriteRegion(oldImage, currentRegion, "CurrentRegion.png");
    }

    // Get all patches containing the currentPixel
    std::vector<itk::ImageRegion<2> > patchesContainingPixel =
      ITKHelpers::GetAllPatchesContainingPixel(currentPixel,
                                               this->PatchRadius,
                                               fullRegion);

    // Remove patches from the set if they have invalid NNField values
//     patchesContainingPixel.erase(std::remove_if(patchesContainingPixel.begin(), patchesContainingPixel.end(),
//                                  [NearestNeighborField](const itk::ImageRegion<2>& testRegion)
//                                   {
//                                     itk::Index<2> index = ITKHelpers::GetRegionCenter(testRegion);
//                                     bool valid = NearestNeighborField->GetPixel(index).IsValid();
//                                     return !valid;
//                                   }),
//                                  patchesContainingPixel.end());

    assert(patchesContainingPixel.size() > 0);

    // Compute the list of pixels contributing to this patch and their associated patch scores
    std::vector<typename TImage::PixelType> contributingPixels(patchesContainingPixel.size());
    std::vector<float> contributingScores(patchesContainingPixel.size());

    for(unsigned int containingPatchId = 0;
        containingPatchId < patchesContainingPixel.size(); ++containingPatchId)
    {
      itk::Index<2> containingRegionCenter =
                  ITKHelpers::GetRegionCenter(patchesContainingPixel[containingPatchId]);
      Match bestMatch = this->NearestNeighborField->GetPixel(containingRegionCenter);
      assert(bestMatch.IsValid());

      itk::ImageRegion<2> bestMatchRegion = bestMatch.Region;

      //assert(fullRegion.IsInside(bestMatchRegion));

      itk::Index<2> bestMatchRegionCenter = ITKHelpers::GetRegionCenter(bestMatchRegion);

      // Compute the offset of the pixel in question relative to the center of
      // the current patch that contains the pixel
      itk::Offset<2> offset = currentPixel - containingRegionCenter;

      // Compute the location of the pixel in the best matching patch that is the
      // same position of the pixel in question relative to the containing patch
      itk::Index<2> correspondingPixel = bestMatchRegionCenter + offset;

      contributingPixels[containingPatchId] = this->Image->GetPixel(correspondingPixel);

      contributingScores[containingPatchId] = bestMatch.Score;

    } // end loop over containing patches

    //std::cout << "Collected contributing pixels and scores." << std::endl;
    // Compute new pixel value

    // Select a method to construct new pixel
    //std::cout << "Compositing..." << std::endl;
    typename TImage::PixelType newValue = Composite(contributingPixels, contributingScores);
    //std::cout << "Done compositing." << std::endl;

    updatedImage->SetPixel(currentPixel, newValue);

//         std::cout << "Pixel was " << currentImage->GetPixel(currentPixel)
//                   << " and is now " << updateImage->GetPixel(currentPixel) << std::endl;

  } // end loop over all target pixels

  // Actually update the output (all at the same time)
  ITKHelpers::DeepCopy(updatedImage.GetPointer(), this->Output.GetPointer());
  std::cout << "Finished Compositor::Compute()." << std::endl;
}

template <typename TImage>
void Compositor<TImage>::SetCompositingMethod(const CompositingMethodEnum& compositingMethod)
{
  this->CompositingMethod = compositingMethod;
}

/** Composite using the specified CompositingMethod. */
template <typename TImage>
typename TImage::PixelType Compositor<TImage>::Composite(
    const std::vector<typename TImage::PixelType>& contributingPixels,
    const std::vector<float>& contributingScores)
{
  assert(contributingPixels.size() == contributingScores.size());
  assert(contributingPixels.size() > 0);

  if(this->CompositingMethod == AVERAGE)
  {
    return CompositeAverage(contributingPixels);
  }
  else if(this->CompositingMethod == WEIGHTED_AVERAGE)
  {
    return CompositeWeightedAverage(contributingPixels, contributingScores);
  }
  else if(this->CompositingMethod == CLOSEST_TO_AVERAGE)
  {
    return CompositeClosestToAverage(contributingPixels, contributingScores);
  }
  else if(this->CompositingMethod == BEST_PATCH)
  {
    return CompositeBestPatch(contributingPixels, contributingScores);
  }
  else
  {
    throw std::runtime_error("An invalid CompositingMethod was selected!");
  }
}

/** Composite by averaging pixels. */
template <typename TImage>
typename TImage::PixelType Compositor<TImage>::CompositeAverage(
  const std::vector<typename TImage::PixelType>& contributingPixels)
{
  typename TImage::PixelType newValue = Statistics::Average(contributingPixels);
  return newValue;
}

/** Composite by weighted averaging. */
template <typename TImage>
typename TImage::PixelType Compositor<TImage>::CompositeWeightedAverage(
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

template <typename TImage>
typename TImage::PixelType Compositor<TImage>::CompositeClosestToAverage(
    const std::vector<typename TImage::PixelType>& contributingPixels,
    const std::vector<float>& contributingScores)
{
  // Use the pixel closest to the average pixel
  typename TImage::PixelType averagePixel = Statistics::Average(contributingPixels);
  unsigned int patchId = ITKHelpers::ClosestValueIndex(contributingPixels, averagePixel);
  typename TImage::PixelType newValue = contributingPixels[patchId];
  return newValue;
}

template <typename TImage>
typename TImage::PixelType Compositor<TImage>::CompositeBestPatch(
    const std::vector<typename TImage::PixelType>& contributingPixels,
    const std::vector<float>& contributingScores)
{
  // Take the pixel from the best matching patch
  unsigned int patchId = Helpers::argmin(contributingScores);
  typename TImage::PixelType newValue = contributingPixels[patchId];
  return newValue;
}

template <typename TImage>
void Compositor<TImage>::SetNearestNeighborField(NNFieldType* const nnField)
{
  this->NearestNeighborField = nnField;
}

#endif
