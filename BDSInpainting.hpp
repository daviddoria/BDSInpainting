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
#include <PatchMatch/InitializerKnownRegion.h>
#include <PatchMatch/InitializerRandom.h>

#include <PatchComparison/SSD.h>

// ITK
#include "itkImageRegionReverseIterator.h"

// STL
#include <ctime>

template <typename TImage, typename TPatchMatchFunctor>
BDSInpainting<TImage, TPatchMatchFunctor>::BDSInpainting() : InpaintingAlgorithm<TImage, TPatchMatchFunctor>()
{
}

template <typename TImage, typename TPatchMatchFunctor>
void BDSInpainting<TImage, TPatchMatchFunctor>::Inpaint()
{
  assert(this->Image->GetLargestPossibleRegion().GetSize()[0] > 0);
  assert(this->SourceMask->GetLargestPossibleRegion().GetSize()[0] > 0);
  assert(this->TargetMask->GetLargestPossibleRegion().GetSize()[0] > 0);

  // Initialize the output with the input
  typename TImage::Pointer currentImage = TImage::New();
  ITKHelpers::DeepCopy(this->Image.GetPointer(), currentImage.GetPointer());

  // This is done so that the "full region" (the entire images) does not have to be
  // referenced using a particular image or mask. That is, 'SourceMask->GetLargestPossibleRegion()'
  // would not raise the question "Why is this region coming from the SourceMask?"
  itk::ImageRegion<2> fullRegion = this->Image->GetLargestPossibleRegion();

  for(unsigned int iteration = 0; iteration < this->Iterations; ++iteration)
  {
    // Allocate the initial NNField
    typename TPatchMatchFunctor::MatchImageType::Pointer matchImage =
      TPatchMatchFunctor::MatchImageType::New();
    matchImage->SetRegions(currentImage->GetLargestPossibleRegion());
    matchImage->Allocate();

    // Initialize the NNField in the known region
    InitializerKnownRegion initializerKnownRegion;
    initializerKnownRegion.SetSourceMask(this->SourceMask);
    initializerKnownRegion.SetPatchRadius(this->PatchRadius);
    initializerKnownRegion.Initialize(matchImage);

    PatchMatchHelpers::WriteNNField(matchImage.GetPointer(),
                                        "InitializedKnownRegionNNField.mha"); // debug only

    // Initialize the NNField in the target region
    InitializerRandom<typename std::remove_pointer<decltype(this->PatchMatchFunctor->GetPatchDistanceFunctor())>::type> randomInitializer;
    randomInitializer.SetTargetMask(this->TargetMask);
    randomInitializer.SetSourceMask(this->SourceMask);
    randomInitializer.SetPatchDistanceFunctor(this->PatchMatchFunctor->GetPatchDistanceFunctor());
    randomInitializer.SetPatchRadius(this->PatchRadius);
    randomInitializer.Initialize(matchImage);

    PatchMatchHelpers::WriteNNField(matchImage.GetPointer(),
                                        "InitializedRandomNNField.mha"); // debug only
    // Give the PatchMatch functor the data
    this->PatchMatchFunctor->SetSourceMask(this->SourceMask);
    this->PatchMatchFunctor->SetTargetMask(this->TargetMask);
    this->PatchMatchFunctor->SetInitialNNField(matchImage);
    this->PatchMatchFunctor->Compute();

    PatchMatchHelpers::WriteNNField(this->PatchMatchFunctor->GetOutput(), "PatchMatchNNField.mha");

    // Update the target pixels
    this->CompositorFunctor->SetTargetMask(this->TargetMask);
    this->CompositorFunctor->SetImage(currentImage);
    this->CompositorFunctor->SetPatchRadius(this->PatchRadius);
    this->CompositorFunctor->SetNearestNeighborField(this->PatchMatchFunctor->GetOutput());
    this->CompositorFunctor->Compute();

    ITKHelpers::DeepCopy(this->CompositorFunctor->GetOutput(), currentImage.GetPointer());
  }

  ITKHelpers::DeepCopy(currentImage.GetPointer(), this->Output.GetPointer());

  ITKHelpers::WriteRGBImage(this->Output.GetPointer(), "ComputeOutput.png");
}

#endif
