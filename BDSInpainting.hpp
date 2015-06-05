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
#include <ITKHelpers/ITKHelpers.h>

#include <Mask/MaskOperations.h>

#include <PatchMatch/PatchMatch.h>
#include <PatchMatch/PatchMatchHelpers.h>
#include <PatchMatch/Propagator.h>

#include <PatchComparison/SSD.h>

// ITK
#include "itkImageRegionReverseIterator.h"

// STL
#include <ctime>

template <typename TImage>
template <typename TPatchMatchFunctor, typename TCompositor>
void BDSInpainting<TImage>::Inpaint(TPatchMatchFunctor* const patchMatchFunctor,
                                    TCompositor* const compositor)
{
  assert(this->Image);
  assert(this->InpaintingMask);

  ConstructValidPatchCentersImage();

  // Initialize the output with the input
  typename TImage::Pointer currentImage = TImage::New();
  ITKHelpers::DeepCopy(this->Image.GetPointer(), currentImage.GetPointer());

  // Allocate the initial NNField
  NNFieldType::Pointer nnField =
    NNFieldType::New();
  nnField->SetRegions(currentImage->GetLargestPossibleRegion());
  nnField->Allocate();

  // Initialize the NNField in the target region
  typedef SSD<TImage> PatchDistanceFunctorType;
  PatchDistanceFunctorType patchDistanceFunctor;
  patchDistanceFunctor.SetImage(currentImage);

  patchMatchFunctor->SetValidPatchCentersImage(this->ValidPatchCentersImage);

  std::vector<itk::Index<2> > pixelsToProcess = this->InpaintingMask->GetHolePixels();
  patchMatchFunctor->SetTargetPixels(pixelsToProcess);
  patchMatchFunctor->SetPatchRadius(this->PatchRadius);

  patchMatchFunctor->GetPropagationFunctor()->SetPatchDistanceFunctor(&patchDistanceFunctor);
  patchMatchFunctor->GetPropagationFunctor()->SetPatchRadius(this->PatchRadius);

  patchMatchFunctor->GetRandomSearchFunctor()->SetImage(this->Image);
  patchMatchFunctor->GetRandomSearchFunctor()->SetPatchRadius(this->PatchRadius);
  patchMatchFunctor->GetRandomSearchFunctor()->SetPatchDistanceFunctor(&patchDistanceFunctor);

  compositor->SetPatchRadius(this->PatchRadius);
  compositor->SetTargetMask(this->InpaintingMask);
  compositor->SetImage(currentImage);
  compositor->SetNearestNeighborField(patchMatchFunctor->GetNNField());

  for(unsigned int iteration = 0; iteration < this->Iterations; ++iteration)
  {
    // Run PatchMatch to compute the NNField
    patchMatchFunctor->Compute();

    std::stringstream ssNNFieldFileName;
    ssNNFieldFileName << "BDS_" << iteration << "_NNField.mha";
    PatchMatchHelpers::WriteNNField(patchMatchFunctor->GetNNField(), ssNNFieldFileName.str());

    // Update the target pixels
    compositor->Composite();
    ITKHelpers::DeepCopy(compositor->GetOutput(), currentImage.GetPointer());
  }

  ITKHelpers::DeepCopy(currentImage.GetPointer(), this->Output.GetPointer());
}

template <typename TImage>
void BDSInpainting<TImage>::ConstructValidPatchCentersImage()
{
    this->ValidPatchCentersImage->SetRegions(this->Image->GetLargestPossibleRegion());
    this->ValidPatchCentersImage->Allocate();

    ITKHelpers::SetImageToConstant(this->ValidPatchCentersImage.GetPointer(), false);

    itk::ImageRegion<2> internalRegion = ITKHelpers::GetInternalRegion(this->InpaintingMask->GetLargestPossibleRegion(),
                                              this->PatchRadius);

//    itk::ImageRegionIteratorWithIndex<Mask> maskIterator(this->InpaintingMask.GetPointer(), internalRegion);
    itk::ImageRegionIteratorWithIndex<BoolImageType> iterator(this->ValidPatchCentersImage.GetPointer(), internalRegion);

    while(!iterator.IsAtEnd())
    {
      itk::ImageRegion<2> targetRegion = ITKHelpers::GetRegionInRadiusAroundPixel(iterator.GetIndex(), this->PatchRadius);

      if(this->InpaintingMask->IsValid(targetRegion))
      {
          iterator.Set(true);
      }

      ++iterator;
    }

    ITKHelpers::WriteBoolImage(this->ValidPatchCentersImage.GetPointer(), "ValidPatchCentersImage.png");
}

#endif
