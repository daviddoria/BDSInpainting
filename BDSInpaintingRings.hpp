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
                                         typename PatchMatch<TImage>::PMImageType* inputNNField, TImage* const output)
{
  ITKHelpers::WriteImage(this->TargetMask.GetPointer(), "OriginalTargetMask.png");

  // Save the original mask, as we will be modifying the internal masks below
  Mask::Pointer currentTargetMask = Mask::New();
  currentTargetMask->DeepCopyFrom(this->TargetMask);

  Mask::Pointer currentPropagationMask = Mask::New();
  currentPropagationMask->DeepCopyFrom(this->TargetMask);
  currentPropagationMask->Invert();

  typename PatchMatch<TImage>::PMImageType::Pointer previousNNField = NULL;
  if(inputNNField)
  {
    ITKHelpers::DeepCopy(inputNNField, previousNNField.GetPointer());
  }

  unsigned int ringCounter = 0;
  typename TImage::Pointer filledRing = TImage::New();

  while(currentTargetMask->HasValidPixels())
  {
    {
    std::stringstream ssCurrentTargetMask;
    ssCurrentTargetMask << "CurrentTargetMask_" << Helpers::ZeroPad(ringCounter, 4) << ".png";
    ITKHelpers::WriteImage(currentTargetMask.GetPointer(), ssCurrentTargetMask.str());
    }

    // Get the outside boundary of the target region
    Mask::BoundaryImageType::Pointer boundaryImage = Mask::BoundaryImageType::New();
    Mask::BoundaryImageType::PixelType boundaryColor = 255;
    currentTargetMask->FindBoundary(boundaryImage, Mask::HOLE, boundaryColor);

//     {
//     std::stringstream ssBoundary;
//     ssBoundary << "Boundary_" << Helpers::ZeroPad(ringCounter, 4) << ".png";
//     ITKHelpers::WriteImage(boundaryImage.GetPointer(), ssBoundary.str());
//     }

    // Create a mask of just the boundary
    Mask::Pointer boundaryMask = Mask::New();
    Mask::BoundaryImageType::PixelType holeColor = boundaryColor;
    boundaryMask->CreateFromImage(boundaryImage.GetPointer(), holeColor);
    boundaryMask->Invert(); // Make the thin boundary the only valid pixels

//     { // debug only
//     std::stringstream ssBoundaryMask;
//     ssBoundaryMask << "BoundaryMask_" << Helpers::ZeroPad(ringCounter, 4) << ".png";
//     ITKHelpers::WriteImage(boundaryImage.GetPointer(), ssBoundaryMask.str());
//     }

    // Set the mask to use in the PatchMatch algorithm
    this->TargetMask->DeepCopyFrom(boundaryMask);

    this->PatchMatchFunctor->SetAllowedPropagationMask(currentPropagationMask);
    // Compute the filling in the ring
    BDSInpainting<TImage>::Compute(image, this->SourceMask, boundaryMask, previousNNField, filledRing);

    { // debug only
    std::stringstream ssInpaintedRing;
    ssInpaintedRing << "InpaintedRing_" << Helpers::ZeroPad(ringCounter, 4) << ".png";
    ITKHelpers::WriteImage(filledRing.GetPointer(), ssInpaintedRing.str());
    }

    { // debug only
    typename PatchMatch<TImage>::CoordinateImageType::Pointer coordinateImage = PatchMatch<TImage>::CoordinateImageType::New();
    PatchMatch<TImage>::GetPatchCentersImage(this->PatchMatchFunctor->GetOutput(), coordinateImage);
    std::stringstream ssOutput;
    ssOutput << "NNField_" << Helpers::ZeroPad(ringCounter, 3) << ".mha";
    ITKHelpers::WriteImage(coordinateImage.GetPointer(), ssOutput.str());
    }

    // Copy the filled ring into the image for the next iteration
    ITKHelpers::DeepCopy(filledRing.GetPointer(), image);

    // Reduce the size of the target region (we "enlarge the hole", because the "hole" is considered the valid part of the target mask)
    unsigned int kernelRadius = 1;
    currentTargetMask->ExpandHole(kernelRadius);

    currentPropagationMask->DeepCopyFrom(currentTargetMask);
    currentPropagationMask->Invert();

    if(!previousNNField)
    {
      previousNNField = PatchMatch<TImage>::PMImageType::New();
    }

    ITKHelpers::DeepCopy(this->PatchMatchFunctor->GetOutput(), previousNNField.GetPointer());
    this->PatchMatchFunctor->SetAllowedPropagationMask(currentTargetMask);
    ringCounter++;
  }

  ITKHelpers::DeepCopy(filledRing.GetPointer(), output);

}

#endif
