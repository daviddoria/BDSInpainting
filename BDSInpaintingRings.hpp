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

  // The pixels from which information is allowed to propagate
  Mask::Pointer currentPropagationMask = Mask::New();

  // Use a NULL previous field unless one has been provided. We have to
  // create this as a smart pointer because if the inputNNField pointer is NULL,
  // we cannot allocate it later because it is not a smart pointer.
  typename PatchMatch<TImage>::PMImageType::Pointer previousNNField = NULL;
  if(inputNNField)
  {
    ITKHelpers::DeepCopy(inputNNField, previousNNField.GetPointer());
  }

  unsigned int ringCounter = 0;

  // This image will be used to store the filled portion of the image (a ring at a time)
  typename TImage::Pointer filledRing = TImage::New();

  while(currentTargetMask->HasValidPixels())
  {
    { // debug only
    std::stringstream ssCurrentTargetMask;
    ssCurrentTargetMask << "CurrentTargetMask_" << Helpers::ZeroPad(ringCounter, 4) << ".png";
    ITKHelpers::WriteImage(currentTargetMask.GetPointer(), ssCurrentTargetMask.str());
    }

    // We trust the information everywhere except in the hole
    currentPropagationMask->DeepCopyFrom(this->TargetMask);
    currentPropagationMask->Invert();

    ITKHelpers::WriteSequentialImage(currentPropagationMask.GetPointer(), "CurrentPropagationMask", ringCounter, 4, "png");

    // Get the inside boundary of the target region
    Mask::BoundaryImageType::Pointer boundaryImage = Mask::BoundaryImageType::New();
    Mask::BoundaryImageType::PixelType boundaryValue = 255; // In the resulting boundary image, the boundary will be 255.
    currentTargetMask->FindBoundary(boundaryImage, Mask::VALID, boundaryValue);

    ITKHelpers::WriteSequentialImage(boundaryImage.GetPointer(), "Boundary", ringCounter, 4, "png");

    // Create a mask of just the boundary
    Mask::Pointer boundaryMask = Mask::New();
    Mask::BoundaryImageType::PixelType holeValue = boundaryValue;
    Mask::BoundaryImageType::PixelType validValue = 0;
    boundaryMask->CreateFromImage(boundaryImage.GetPointer(), holeValue, validValue);
    //boundaryMask->Invert(); // Make the thin boundary the only valid pixels in the mask

    ITKHelpers::WriteSequentialImage(boundaryMask.GetPointer(), "BoundaryMask", ringCounter, 4, "png");

    // Set the mask to use in the PatchMatch algorithm
    this->TargetMask->DeepCopyFrom(boundaryMask);

    this->PatchMatchFunctor->SetAllowedPropagationMask(currentPropagationMask);
    // Compute the filling in the ring
    BDSInpainting<TImage>::Compute(image, this->SourceMask, boundaryMask, previousNNField, filledRing);

    ITKHelpers::WriteSequentialImage(filledRing.GetPointer(), "InpaintedRing", ringCounter, 4, "png");

    { // debug only
    typename PatchMatch<TImage>::CoordinateImageType::Pointer coordinateImage = PatchMatch<TImage>::CoordinateImageType::New();
    PatchMatch<TImage>::GetPatchCentersImage(this->PatchMatchFunctor->GetOutput(), coordinateImage);
    ITKHelpers::WriteSequentialImage(coordinateImage.GetPointer(), "NNField", ringCounter, 4, "mha");
    }

    // Copy the filled ring into the image for the next iteration
    ITKHelpers::DeepCopy(filledRing.GetPointer(), image);

    // Reduce the size of the target region (we "enlarge the hole", because the "hole" is considered the valid part of the target mask)
    unsigned int kernelRadius = 1;
    currentTargetMask->ExpandHole(kernelRadius);

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