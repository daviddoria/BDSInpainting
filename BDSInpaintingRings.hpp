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

#ifndef BDSInpaintingRings_HPP
#define BDSInpaintingRings_HPP

#include "BDSInpaintingRings.h"

template <typename TImage>
BDSInpaintingRings<TImage>::BDSInpaintingRings() : BDSInpainting<TImage>()
{

}

template <typename TImage>
void BDSInpaintingRings<TImage>::Compute(TImage* const image, Mask* const sourceMask, Mask* const targetMask,
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
    this->PatchMatchFunctor->SetImage(currentImage);
    this->PatchMatchFunctor->SetSourceMask(sourceMask);
    this->PatchMatchFunctor->SetTargetMask(targetMask);

    try
    {
      if(iteration == 0)
      {
        this->PatchMatchFunctor->Compute(NULL);
      }
      else
      {
        // For now don't initialize with the previous NN field - though this
        // might work and be a huge speed up.
        this->PatchMatchFunctor->Compute(NULL);
        //this->PatchMatchFunctor.Compute(init);
      }
    }
    catch (std::runtime_error &e)
    {
      std::cout << e.what() << std::endl;
      return;
    }

    typename PatchMatch<TImage>::PMImageType* nnField = this->PatchMatchFunctor->GetOutput();

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

#endif
