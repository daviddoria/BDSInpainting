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

#ifndef BDSInpainting_HPP
#define BDSInpainting_HPP

#include "BDSInpainting.h"

// Submodules
#include "ITKHelpers/ITKHelpers.h"
#include "Mask/MaskOperations.h"
#include "PatchMatch/PatchMatch.h"
#include "PatchComparison/SSD.h"
#include "PoissonEditing/PoissonEditing.h"

// ITK
#include "itkImageRegionReverseIterator.h"

// STL
#include <ctime>

template <typename TImage>
BDSInpainting<TImage>::BDSInpainting() : ResolutionLevels(3), Iterations(5),
                                         PatchRadius(7), DownsampleFactor(.5),
                                         CompositingMethod(WEIGHTED_AVERAGE)
{
  this->Output = TImage::New();
  this->Image = TImage::New();
  this->SourceMask = Mask::New();
  this->TargetMask = Mask::New();
}

template <typename TImage>
void BDSInpainting<TImage>::Compute()
{

  // The finest scale masks and image are simply the user inputs.
  Mask::Pointer level0sourceMask = Mask::New();
  level0sourceMask->DeepCopyFrom(this->SourceMask);

  Mask::Pointer level0targetMask = Mask::New();
  level0targetMask->DeepCopyFrom(this->TargetMask);

  typename TImage::Pointer level0Image = TImage::New();
  ITKHelpers::DeepCopy(this->Image.GetPointer(), level0Image.GetPointer());

  std::vector<typename TImage::Pointer> imageLevels(this->ResolutionLevels);
  std::vector<Mask::Pointer> sourceMaskLevels(this->ResolutionLevels);
  std::vector<Mask::Pointer> targetMaskLevels(this->ResolutionLevels);

  imageLevels[0] = level0Image;
  sourceMaskLevels[0] = level0sourceMask;
  targetMaskLevels[0] = level0targetMask;

  // Downsample the image and mask to the number of specified resolutions.
  // Start at level 1 because 0 is the full resolution (provided directly by the user)
  for(unsigned int level = 1; level < this->ResolutionLevels; ++level)
  {
    // At each level we want the image to be the previous level size modified by the downsample factor
    itk::Size<2> destinationSize;
    destinationSize[0] = imageLevels[level-1]->GetLargestPossibleRegion().GetSize()[0] *
                         this->DownsampleFactor;
    destinationSize[1] = imageLevels[level-1]->GetLargestPossibleRegion().GetSize()[1] *
                         this->DownsampleFactor;

    // Downsample the image
    typename TImage::Pointer downsampledImage = TImage::New();
    ITKHelpers::ScaleImage(imageLevels[level - 1].GetPointer(), destinationSize,
                           downsampledImage.GetPointer());
    imageLevels[level] = downsampledImage;

    // Downsample the source mask
    Mask::Pointer downsampledSourceMask = Mask::New();
    ITKHelpers::ScaleImage(sourceMaskLevels[level - 1].GetPointer(), destinationSize,
                           downsampledSourceMask.GetPointer());
    downsampledSourceMask->CopyInformationFrom(this->SourceMask);
    sourceMaskLevels[level] = downsampledSourceMask;

    // Downsample the target mask
    Mask::Pointer downsampledTargetMask = Mask::New();
    ITKHelpers::ScaleImage(targetMaskLevels[level - 1].GetPointer(), destinationSize,
                           downsampledTargetMask.GetPointer());
    downsampledTargetMask->CopyInformationFrom(this->TargetMask);
    targetMaskLevels[level] = downsampledTargetMask;
  }

  // Debug only - write the images at every level
  for(unsigned int level = 0; level < this->ResolutionLevels; ++level)
  {
    std::cout << "Level " << level << " image resolution "
              << imageLevels[level]->GetLargestPossibleRegion().GetSize() << std::endl;

    std::cout << "Level " << level << " source mask resolution "
              << sourceMaskLevels[level]->GetLargestPossibleRegion().GetSize() << std::endl;

    std::stringstream ssImage;
    ssImage << "Input_Level_" << Helpers::ZeroPad(level, 2) << ".png";
    ITKHelpers::WriteRGBImage(imageLevels[level].GetPointer(), ssImage.str());

    std::stringstream ssMask;
    ssMask << "SourceMask_Level_" << Helpers::ZeroPad(level, 2) << ".png";
    ITKHelpers::WriteImage(sourceMaskLevels[level].GetPointer(), ssMask.str());
  }

  // Compute the filling, starting at the lowest resolution and working back to the original resolution.
  // At each level, use the output of the previous level as initialization to the next level.
  for(unsigned int level = this->ResolutionLevels - 1; level >= 0; --level)
  {
    std::cout << "BDS level " << level << " (resolution "
              << imageLevels[level]->GetLargestPossibleRegion().GetSize() << ")" << std::endl;
    typename TImage::Pointer output = TImage::New();
    Compute(imageLevels[level].GetPointer(), sourceMaskLevels[level].GetPointer(),
            targetMaskLevels[level].GetPointer(), output);

    { // Debug only
    std::stringstream ss;
    ss << "Output_Level_" << Helpers::ZeroPad(level, 2) << ".png";
    ITKHelpers::WriteRGBImage(output.GetPointer(), ss.str());
    }

    // If this is the finest resolution, we are done
    if(level == 0)
    {
      ITKHelpers::DeepCopy(output.GetPointer(), this->Output.GetPointer());
      break;
    }

    // Upsample result and copy it to the next level.
    typename TImage::Pointer upsampled = TImage::New();
    ITKHelpers::ScaleImage(output.GetPointer(), imageLevels[level-1]->GetLargestPossibleRegion().GetSize(),
                           upsampled.GetPointer());

    { // Debug only
//     std::cout << "Upsampled from " << output->GetLargestPossibleRegion().GetSize() << " to "
//               << upsampled->GetLargestPossibleRegion().GetSize() << std::endl;
//
//     std::cout << "Upsampled size: " << upsampled->GetLargestPossibleRegion().GetSize() << std::endl;
//     std::cout << "Next level size: " << imageLevels[level - 1]->GetLargestPossibleRegion().GetSize()
//               << std::endl;
    }

    // Only keep the computed pixels in the hole - the rest of the pixels are simply from one level up.
    MaskOperations::CopyInHoleRegion(upsampled.GetPointer(),
                                     imageLevels[level - 1].GetPointer(),
                                     targetMaskLevels[level - 1].GetPointer());

    std::stringstream ssUpdated;
    ssUpdated << "UpdatedInput_Level_" << Helpers::ZeroPad(level - 1, 2) << ".png";
    ITKHelpers::WriteRGBImage(imageLevels[level - 1].GetPointer(), ssUpdated.str());
  }

}

template <typename TImage>
void BDSInpainting<TImage>::Compute(TImage* const image, Mask* const sourceMask, Mask* const targetMask,
                                    TImage* const output)
{
  ITKHelpers::WriteRGBImage(image, "ComputeInput.png");

  // Initialize the output with the input
  ITKHelpers::DeepCopy(image, output);

  // Initialize the image to operate on
  typename TImage::Pointer currentImage = TImage::New();
  ITKHelpers::DeepCopy(image, currentImage.GetPointer());

  itk::ImageRegion<2> fullRegion = image->GetLargestPossibleRegion();

  std::cout << "Computing BDS on resolution " << fullRegion.GetSize() << std::endl;

  for(unsigned int iteration = 0; iteration < this->Iterations; ++iteration)
  {
    std::cout << "BDSInpainting Iteration " << iteration << std::endl;

    // Give the PatchMatch functor the data
    this->PatchMatchFunctor.SetImage(currentImage);
    this->PatchMatchFunctor.SetSourceMask(sourceMask);
    this->PatchMatchFunctor.SetTargetMask(targetMask);

    try
    {
      if(iteration == 0)
      {
        this->PatchMatchFunctor.Compute(NULL);
      }
      else
      {
        // For now don't initialize with the previous NN field - though this
        // might work and be a huge speed up.
        this->PatchMatchFunctor.Compute(NULL);
        //this->PatchMatchFunctor.Compute(init);
      }
    }
    catch (std::runtime_error &e)
    {
      std::cout << e.what() << std::endl;
      return;
    }

    typename PatchMatch<TImage>::PMImageType* nnField = this->PatchMatchFunctor.GetOutput();

    ITKHelpers::WriteRGBImage(currentImage.GetPointer(), "BeforeUpdate.png");
    
    // Update the target pixels
    typename TImage::Pointer updatedImage = TImage::New();
    //std::cout << "Updating pixels..." << std::endl;
    UpdatePixels(currentImage, targetMask, nnField, updatedImage);
    //std::cout << "Done updating pixels." << std::endl;
    ITKHelpers::WriteRGBImage(updatedImage.GetPointer(), "Updated.png");
    ITKHelpers::WriteRGBImage(currentImage.GetPointer(), "Current.png");

    MaskOperations::CopyInValidRegion(updatedImage.GetPointer(), currentImage.GetPointer(), targetMask);

    std::stringstream ssPNG;
    ssPNG << "BDS_Iteration_" << Helpers::ZeroPad(iteration, 2) << ".png";
    ITKHelpers::WriteRGBImage(currentImage.GetPointer(), ssPNG.str());
  } // end iterations loop

  ITKHelpers::DeepCopy(currentImage.GetPointer(), output);

  ITKHelpers::WriteRGBImage(output, "ComputeOutput.png");
}

template <typename TImage>
TImage* BDSInpainting<TImage>::GetOutput()
{
  return this->Output;
}

template <typename TImage>
void BDSInpainting<TImage>::SetIterations(const unsigned int iterations)
{
  this->Iterations = iterations;
}

template <typename TImage>
void BDSInpainting<TImage>::SetPatchRadius(const unsigned int patchRadius)
{
  this->PatchRadius = patchRadius;
}

template <typename TImage>
void BDSInpainting<TImage>::SetImage(TImage* const image)
{
  ITKHelpers::DeepCopy(image, this->Image.GetPointer());
}

template <typename TImage>
void BDSInpainting<TImage>::SetSourceMask(Mask* const mask)
{
  this->SourceMask->DeepCopyFrom(mask);
}

template <typename TImage>
void BDSInpainting<TImage>::SetTargetMask(Mask* const mask)
{
  this->TargetMask->DeepCopyFrom(mask);
}

template <typename TImage>
void BDSInpainting<TImage>::SetResolutionLevels(const unsigned int resolutionLevels)
{
  if(resolutionLevels < 1)
  {
    std::cerr << "ResolutionLevels: " << resolutionLevels << std::endl;
    throw std::runtime_error("ResolutionLevels must be >= 1!");
  }
  this->ResolutionLevels = resolutionLevels;
}

template <typename TImage>
void BDSInpainting<TImage>::SetPatchMatchFunctor(const PatchMatch<TImage>& patchMatchFunctor)
{
  this->PatchMatchFunctor = patchMatchFunctor;
}

template <typename TImage>
void BDSInpainting<TImage>::SetDownsampleFactor(const float downsampleFactor)
{
  this->DownsampleFactor = downsampleFactor;
}

template <typename TImage>
void BDSInpainting<TImage>::UpdatePixels(const TImage* const oldImage,
                                         const Mask* const targetMask,
                                         const typename PatchMatch<TImage>::PMImageType* const nnField,
                                         TImage* const updatedImage)
{
  // The contribution of each pixel q to the error term (d_cohere) = 1/N_T \sum_{i=1}^m (S(p_i) - T(q))^2
  // To find the best color T(q) (iterative update rule), differentiate with respect to T(q),
  // set to 0, and solve for T(q):
  // T(q) = \frac{1}{m} \sum_{i=1}^m S(p_i)

  ITKHelpers::WriteRGBImage(oldImage, "UpdatePixels_OldImage.png");

  // This is done so in the algorithm we can use 'fullRegion', since it refers
  // to the same region for the image and mask.
  // (So there is no confusion such as "why is the mask's region used here instead of the image's?")
  itk::ImageRegion<2> fullRegion = oldImage->GetLargestPossibleRegion();

  // Only iterate over this region, as it is the only region where a hole pixel can be found
  itk::ImageRegion<2> targetBoundingBox = MaskOperations::ComputeValidBoundingBox(targetMask);

  // We don't want to change pixels directly on the
  // output image during the iteration, but rather compute them all and then update them all simultaneously.
  ITKHelpers::DeepCopy(oldImage, updatedImage);

  // We must get a dummy pixel from the image and then fill it with zero to make sure the number
  // of components of the pixel is correct.
  itk::Index<2> zeroIndex = {{0,0}};
  typename TImage::PixelType zeroPixel = oldImage->GetPixel(zeroIndex);
  zeroPixel.Fill(0);

  itk::ImageRegionIteratorWithIndex<TImage> imageIterator(updatedImage, targetBoundingBox);

  unsigned int pixelCounter = 0;
  while(!imageIterator.IsAtEnd())
  {
    itk::Index<2> currentPixel = imageIterator.GetIndex();
    if(targetMask->IsValid(currentPixel)) // We have come across a pixel to be filled
    {

      { // debug only
      //std::cout << "Updating " << pixelCounter << " of "
      //          << holeBoundingBox.GetNumberOfPixels() << std::endl;
      itk::ImageRegion<2> currentRegion =
            ITKHelpers::GetRegionInRadiusAroundPixel(currentPixel, this->PatchRadius);

      //ITKHelpers::WriteRegion(oldImage, currentRegion, "CurrentRegion.png");
      }

      std::vector<itk::ImageRegion<2> > patchesContainingPixel =
            ITKHelpers::GetAllPatchesContainingPixel(currentPixel,
                                                     this->PatchRadius,
                                                     fullRegion);

      // Compute the list of pixels contributing to this patch and their associated patch scores
      std::vector<typename TImage::PixelType> contributingPixels(patchesContainingPixel.size());
      std::vector<float> contributingScores(patchesContainingPixel.size());

      for(unsigned int containingPatchId = 0;
          containingPatchId < patchesContainingPixel.size(); ++containingPatchId)
      {
        itk::Index<2> containingRegionCenter =
                    ITKHelpers::GetRegionCenter(patchesContainingPixel[containingPatchId]);
        Match bestMatch = nnField->GetPixel(containingRegionCenter);
        itk::ImageRegion<2> bestMatchRegion = bestMatch.Region;
        itk::Index<2> bestMatchRegionCenter = ITKHelpers::GetRegionCenter(bestMatchRegion);

        { // debug only
//         std::cout << "Containing region center: " << containingRegionCenter << std::endl;
//         std::stringstream ssContainingRegionFile;
//         ssContainingRegionFile << "ContainingRegion_" << containingPatchId << ".png";
//         ITKHelpers::WriteRegion(oldImage, patchesContainingPixel[containingPatchId],
//                                 ssContainingRegionFile.str());

//         std::cout << "Matching region center: " << bestMatchRegionCenter << std::endl;
//         std::stringstream ssMatchingRegionFile;
//         ssMatchingRegionFile << "MatchingRegion_" << containingPatchId << ".png";
        //ITKHelpers::WriteRegion(oldImage, bestMatchRegion, ssMatchingRegionFile.str());
        }

        assert(sourceMask->IsValid(bestMatchRegion));
//         std::cout << "containingRegionCenter: " << containingRegionCenter << std::endl;
//         std::cout << "bestMatchRegionCenter: " << bestMatchRegionCenter << std::endl;

        // Compute the offset of the pixel in question relative to the center of
        // the current patch that contains the pixel
        itk::Offset<2> offset = currentPixel - containingRegionCenter;

        // Compute the location of the pixel in the best matching patch that is the
        // same position of the pixel in question relative to the containing patch
        itk::Index<2> correspondingPixel = bestMatchRegionCenter + offset;

        contributingPixels[containingPatchId] = oldImage->GetPixel(correspondingPixel);
        contributingScores[containingPatchId] = bestMatch.Score;

      } // end loop over containing patches

      // Compute new pixel value

      // Select a method to construct new pixel
      //std::cout << "Compositing..." << std::endl;
      typename TImage::PixelType newValue = Composite(contributingPixels, contributingScores);
      //std::cout << "Done compositing." << std::endl;

      updatedImage->SetPixel(currentPixel, newValue);

//         std::cout << "Pixel was " << currentImage->GetPixel(currentPixel)
//                   << " and is now " << updateImage->GetPixel(currentPixel) << std::endl;
    } // end if is hole

    ++imageIterator;
    pixelCounter++;
  } // end loop over image
}

template <typename TImage>
void BDSInpainting<TImage>::SetCompositingMethod(const CompositingMethodEnum& compositingMethod)
{
  this->CompositingMethod = compositingMethod;
}

/** Composite using the specified CompositingMethod. */
template <typename TImage>
typename TImage::PixelType BDSInpainting<TImage>::Composite(
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
typename TImage::PixelType BDSInpainting<TImage>::CompositeAverage(
  const std::vector<typename TImage::PixelType>& contributingPixels)
{
  typename TImage::PixelType newValue = Statistics::Average(contributingPixels);
  return newValue;
}

/** Composite by weighted averaging. */
template <typename TImage>
typename TImage::PixelType BDSInpainting<TImage>::CompositeWeightedAverage(
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
typename TImage::PixelType BDSInpainting<TImage>::CompositeClosestToAverage(
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
typename TImage::PixelType BDSInpainting<TImage>::CompositeBestPatch(
    const std::vector<typename TImage::PixelType>& contributingPixels,
    const std::vector<float>& contributingScores)
{
  // Take the pixel from the best matching patch
  unsigned int patchId = Helpers::argmin(contributingScores);
  typename TImage::PixelType newValue = contributingPixels[patchId];
  return newValue;
}

#endif
