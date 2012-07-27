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

// ITK
#include "itkImageRegionReverseIterator.h"

// STL
#include <ctime>

template <typename TImage>
BDSInpainting<TImage>::BDSInpainting() : ResolutionLevels(3), Iterations(5), PatchRadius(7), PatchMatchIterations(3), DownsampleFactor(.5)
{
  this->Output = TImage::New();
  this->Image = TImage::New();
  this->SourceMask = Mask::New();
  this->TargetMask = Mask::New();
}

template <typename TImage>
void BDSInpainting<TImage>::Compute()
{
  Mask::Pointer level0sourceMask = Mask::New();
  level0sourceMask->DeepCopyFrom(this->SourceMask);

  Mask::Pointer level0targetMask = Mask::New();
  level0targetMask->DeepCopyFrom(this->TargetMask);

  typename TImage::Pointer level0Image = TImage::New();
  ITKHelpers::DeepCopy(this->Image.GetPointer(), level0Image.GetPointer());

  // Cannot do this! The same image is added as each element of the vector because New() is only evaluated once!
  // std::vector<TImage::Pointer> imageLevels(this->ResolutionLevels, TImage::New());
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
    itk::Size<2> destinationSize;
    destinationSize[0] = imageLevels[level-1]->GetLargestPossibleRegion().GetSize()[0] * this->DownsampleFactor;
    destinationSize[1] = imageLevels[level-1]->GetLargestPossibleRegion().GetSize()[1] * this->DownsampleFactor;

    // Downsample the image
    typename TImage::Pointer downsampledImage = TImage::New();
    ITKHelpers::ScaleImage(imageLevels[level - 1].GetPointer(), destinationSize, downsampledImage.GetPointer());
    imageLevels[level] = downsampledImage;

    // Downsample the source mask
    Mask::Pointer downsampledSourceMask = Mask::New();
    ITKHelpers::ScaleImage(sourceMaskLevels[level - 1].GetPointer(), destinationSize, downsampledSourceMask.GetPointer());
    downsampledSourceMask->CopyInformationFrom(this->SourceMask);
    sourceMaskLevels[level] = downsampledSourceMask;

    // Downsample the target mask
    Mask::Pointer downsampledTargetMask = Mask::New();
    ITKHelpers::ScaleImage(targetMaskLevels[level - 1].GetPointer(), destinationSize, downsampledTargetMask.GetPointer());
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
    Compute(imageLevels[level].GetPointer(), sourceMaskLevels[level].GetPointer(), targetMaskLevels[level].GetPointer(), output);

    std::stringstream ss;
    ss << "Output_Level_" << Helpers::ZeroPad(level, 2) << ".png";
    ITKHelpers::WriteRGBImage(output.GetPointer(), ss.str());

    if(level == 0)
    {
      ITKHelpers::DeepCopy(output.GetPointer(), this->Output.GetPointer());
      break;
    }

    // Upsample result and copy it to the next level. A factor of 2 goes up one level.
    typename TImage::Pointer upsampled = TImage::New();
    ITKHelpers::ScaleImage(output.GetPointer(), imageLevels[level-1]->GetLargestPossibleRegion().GetSize(), upsampled.GetPointer());
//     std::cout << "Upsampled from " << output->GetLargestPossibleRegion().GetSize() << " to "
//               << upsampled->GetLargestPossibleRegion().GetSize() << std::endl;
//
//     std::cout << "Upsampled size: " << upsampled->GetLargestPossibleRegion().GetSize() << std::endl;
//     std::cout << "Next level size: " << imageLevels[level - 1]->GetLargestPossibleRegion().GetSize() << std::endl;

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
void BDSInpainting<TImage>::Compute(TImage* const image, Mask* const sourceMask, Mask* const targetMask, TImage* const output)
{
  ITKHelpers::WriteRGBImage(image, "ComputeInput.png");

  // Initialize the output with the input
  ITKHelpers::DeepCopy(image, output);

  // Initialize the image to operate on
  typename TImage::Pointer currentImage = TImage::New();
  ITKHelpers::DeepCopy(image, currentImage.GetPointer());

  // We must get a dummy pixel from the image and then fill it with zero to make sure the number
  // of components of the pixel is correct.
  itk::Index<2> zeroIndex = {{0,0}};
  typename TImage::PixelType zeroPixel = image->GetPixel(zeroIndex);
  zeroPixel.Fill(0);

  itk::ImageRegion<2> fullRegion = image->GetLargestPossibleRegion();

  std::cout << "Computing BDS on resolution " << fullRegion.GetSize() << std::endl;

  // Setup the PatchMatch object
  PatchMatch<TImage> patchMatch;
  patchMatch.SetPatchRadius(this->PatchRadius);
  patchMatch.SetImage(currentImage);
  SSD<TImage> ssdFunctor;
  ssdFunctor.SetImage(currentImage);
  patchMatch.SetPatchDistanceFunctor(&ssdFunctor);

  for(unsigned int iteration = 0; iteration < this->Iterations; ++iteration)
  {
    std::cout << "BDSInpainting Iteration " << iteration << std::endl;

    // Set the image here even though it was also set outside the loop,
    // because we want to do this at every iteration.
    patchMatch.SetImage(currentImage);
    patchMatch.SetSourceMask(sourceMask);
    patchMatch.SetTargetMask(targetMask);
    patchMatch.SetIterations(this->PatchMatchIterations);

    try
    {
      if(iteration == 0)
      {
        patchMatch.Compute(NULL);
      }
      else
      {
        // For now don't initialize with the previous NN field - though this might work and be a huge speed up.
        patchMatch.Compute(NULL);
        //patchMatch.Compute(init);
      }
    }
    catch (std::runtime_error &e)
    {
      if(e.what() == std::string("PatchMatch: No valid source regions!"))
      {
        //std::cout << "Caught exception." << std::endl;
        return;
      }
    }

    typename PatchMatch<TImage>::PMImageType* nnField = patchMatch.GetOutput();

    // The contribution of each pixel q to the error term (d_cohere) = 1/N_T \sum_{i=1}^m (S(p_i) - T(q))^2
    // To find the best color T(q) (iterative update rule), differentiate with respect to T(q),
    // set to 0, and solve for T(q):
    // T(q) = \frac{1}{m} \sum_{i=1}^m S(p_i)

    // We don't want to change pixels directly on the
    // output image during the iteration, but rather compute them all and then update them all simultaneously.
    typename TImage::Pointer updateImage = TImage::New();
    ITKHelpers::DeepCopy(currentImage.GetPointer(), updateImage.GetPointer());

    // Loop over the whole image (patch centers)
//     itk::ImageRegion<2> internalRegion =
//              ITKHelpers::GetInternalRegion(fullRegion, this->PatchRadius);

    itk::ImageRegion<2> holeBoundingBox = MaskOperations::ComputeHoleBoundingBox(targetMask);

    itk::ImageRegionIteratorWithIndex<TImage> imageIterator(updateImage,
                                                               holeBoundingBox);

    while(!imageIterator.IsAtEnd())
    {
      itk::Index<2> currentPixel = imageIterator.GetIndex();
      if(targetMask->IsHole(currentPixel)) // We have come across a pixel to be filled
      {
        // Zero the pixel - it will be additively updated
        updateImage->SetPixel(currentPixel, zeroPixel);

        itk::ImageRegion<2> currentRegion =
             ITKHelpers::GetRegionInRadiusAroundPixel(currentPixel, this->PatchRadius);

        std::vector<itk::ImageRegion<2> > patchesContainingPixel =
              ITKHelpers::GetAllPatchesContainingPixel(currentPixel,
                                                       this->PatchRadius,
                                                       fullRegion);

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

//           std::cout << "containingRegionCenter: " << containingRegionCenter << std::endl;
//           std::cout << "bestMatchRegionCenter: " << bestMatchRegionCenter << std::endl;

          // Compute the offset of the pixel in question relative to the center of the current patch that contains the pixel
          itk::Offset<2> offset = currentPixel - containingRegionCenter;

          // Compute the location of the pixel in the best matching patch that is the
          // same position of the pixel in question relative to the containing patch
          itk::Index<2> correspondingPixel = bestMatchRegionCenter + offset;

          contributingPixels[containingPatchId] = currentImage->GetPixel(correspondingPixel);
          contributingScores[containingPatchId] = bestMatch.Score;

        } // end loop over containing patches

        // Compute new pixel value

        // Select a method to construct new pixel
        // TImage::PixelType newValue = ITKStatistics::Average(contributingPixels);

        // TImage::PixelType newValue = Helpers::WeightedSum(contributingPixels, contributingScores);

        // Take the pixel from the best matching patch
        unsigned int patchId = Helpers::argmin(contributingScores);
        typename TImage::PixelType newValue = contributingPixels[patchId];

        // Use the pixel closest to the average pixel
//         TImage::PixelType averagePixel = ITKStatistics::Average(contributingPixels);
        //unsigned int patchId = ITKHelpers::ClosestPoint(contributingPixels, averagePixel);
//         TImage::PixelType newValue = contributingPixels[patchId];

        updateImage->SetPixel(currentPixel, newValue);

//         std::cout << "Pixel was " << currentImage->GetPixel(currentPixel)
//                   << " and is now " << updateImage->GetPixel(currentPixel) << std::endl;
      } // end if is hole

      ++imageIterator;
    } // end loop over image

    MaskOperations::CopyInHoleRegion(updateImage.GetPointer(), currentImage.GetPointer(), targetMask);

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
void BDSInpainting<TImage>::SetPatchMatchIterations(const unsigned int patchMatchIterations)
{
  this->PatchMatchIterations = patchMatchIterations;
}

template <typename TImage>
void BDSInpainting<TImage>::SetDownsampleFactor(const float downsampleFactor)
{
  this->DownsampleFactor = downsampleFactor;
}

#endif
