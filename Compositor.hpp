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

template <typename TImage, typename TPixelCompositor>
Compositor<TImage, TPixelCompositor>::Compositor() : PatchRadius(0)
{
  this->Output = TImage::New();
  this->Image = TImage::New();
  this->TargetMask = Mask::New();
}

template <typename TImage, typename TPixelCompositor>
TImage* Compositor<TImage, TPixelCompositor>::GetOutput()
{
  return this->Output;
}

template <typename TImage, typename TPixelCompositor>
void Compositor<TImage, TPixelCompositor>::SetPatchRadius(const unsigned int patchRadius)
{
  this->PatchRadius = patchRadius;
}

template <typename TImage, typename TPixelCompositor>
void Compositor<TImage, TPixelCompositor>::SetImage(TImage* const image)
{
  ITKHelpers::DeepCopy(image, this->Image.GetPointer());
}

template <typename TImage, typename TPixelCompositor>
void Compositor<TImage, TPixelCompositor>::SetTargetMask(Mask* const mask)
{
  this->TargetMask->DeepCopyFrom(mask);
}

template <typename TImage, typename TPixelCompositor>
void Compositor<TImage, TPixelCompositor>::Composite()
{
  // The contribution of each pixel q to the error term (d_cohere) = 1/N_T \sum_{i=1}^m (S(p_i) - T(q))^2
  // To find the best color T(q) (iterative update rule), differentiate with respect to T(q),
  // set to 0, and solve for T(q):
  // T(q) = \frac{1}{m} \sum_{i=1}^m S(p_i)

  assert(this->NearestNeighborField);
  assert(this->TargetMask);
  assert(this->PatchRadius > 0);
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

      itk::ImageRegion<2> bestMatchRegion = bestMatch.GetRegion();

      //assert(fullRegion.IsInside(bestMatchRegion));

      itk::Index<2> bestMatchRegionCenter = ITKHelpers::GetRegionCenter(bestMatchRegion);

      // Compute the offset of the pixel in question relative to the center of
      // the current patch that contains the pixel
      itk::Offset<2> offset = currentPixel - containingRegionCenter;

      // Compute the location of the pixel in the best matching patch that is the
      // same position of the pixel in question relative to the containing patch
      itk::Index<2> correspondingPixel = bestMatchRegionCenter + offset;

      contributingPixels[containingPatchId] = this->Image->GetPixel(correspondingPixel);

      contributingScores[containingPatchId] = bestMatch.GetScore();

    } // end loop over containing patches

    //std::cout << "Collected contributing pixels and scores." << std::endl;
    // Compute new pixel value

    // Select a method to construct new pixel
    //std::cout << "Compositing..." << std::endl;
    typename TImage::PixelType newValue = TPixelCompositor::Composite(contributingPixels, contributingScores);
    //std::cout << "Done compositing." << std::endl;

    updatedImage->SetPixel(currentPixel, newValue);

//         std::cout << "Pixel was " << currentImage->GetPixel(currentPixel)
//                   << " and is now " << updateImage->GetPixel(currentPixel) << std::endl;

  } // end loop over all target pixels

  // Actually update the output (all at the same time)
  ITKHelpers::DeepCopy(updatedImage.GetPointer(), this->Output.GetPointer());
  std::cout << "Finished Compositor::Compute()." << std::endl;
}

template <typename TImage, typename TPixelCompositor>
void Compositor<TImage, TPixelCompositor>::SetNearestNeighborField(NNFieldType* const nnField)
{
  this->NearestNeighborField = nnField;
}

#endif
