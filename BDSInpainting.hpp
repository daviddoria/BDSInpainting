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
BDSInpainting<TImage>::BDSInpainting() : InpaintingAlgorithm<TImage>()
{
}

template <typename TImage>
template <typename TPatchMatchFunctor, typename TCompositor>
void BDSInpainting<TImage>::Inpaint(TPatchMatchFunctor* const patchMatchFunctor,
                                    TCompositor* const compositor)
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

  // Allocate the initial NNField
  NNFieldType::Pointer nnField =
    NNFieldType::New();
  nnField->SetRegions(currentImage->GetLargestPossibleRegion());
  nnField->Allocate();

  // Initialize the NNField in the target region
  typedef SSD<TImage> PatchDistanceFunctorType;
  PatchDistanceFunctorType patchDistanceFunctor;
  patchDistanceFunctor.SetImage(currentImage);

  //////////////////// The things in this block should probably be passed into this class
  // Set acceptance test to histogram threshold
  // Create the HSV image
  typedef itk::VectorImage<float, 2> HSVImageType;
  HSVImageType::Pointer hsvImage = HSVImageType::New();
  ITKHelpers::ITKImageToHSVImage(currentImage.GetPointer(), hsvImage.GetPointer());
  ITKHelpers::WriteImage(hsvImage.GetPointer(), "HSV.mha");

  typedef Propagator<PatchDistanceFunctorType> PropagatorType;
  PropagatorType propagationFunctor;
  propagationFunctor.SetPatchDistanceFunctor(&patchDistanceFunctor);
  propagationFunctor.SetPatchRadius(this->PatchRadius);

  typedef RandomSearch<TImage, PatchDistanceFunctorType>
    RandomSearchType;
  RandomSearchType randomSearcher;
  randomSearcher.SetImage(this->Image);
  randomSearcher.SetPatchDistanceFunctor(&patchDistanceFunctor);
  //////////////////// The things above this block should probably be passed into this class

  // Setup the PatchMatch functor. Use a generic (parent class) AcceptanceTest.
  patchMatchFunctor->SetIterations(5);

  for(unsigned int iteration = 0; iteration < this->Iterations; ++iteration)
  {
    PatchMatchHelpers::WriteNNField(nnField.GetPointer(),
                                        "InitializedKnownRegionNNField.mha"); // debug only

    PatchMatchHelpers::WriteNNField(nnField.GetPointer(),
                                        "InitializedRandomNNField.mha"); // debug only
    // Give the PatchMatch functor the data
    patchMatchFunctor->Compute();

    PatchMatchHelpers::WriteNNField(nnField.GetPointer(), "PatchMatchNNField.mha");

    // Update the target pixels
    compositor->SetTargetMask(this->TargetMask);
    compositor->SetImage(currentImage);
    compositor->SetPatchRadius(this->PatchRadius);
    compositor->SetNearestNeighborField(nnField);
    compositor->Composite();

    ITKHelpers::DeepCopy(compositor->GetOutput(), currentImage.GetPointer());
    ITKHelpers::ITKImageToHSVImage(currentImage.GetPointer(), hsvImage.GetPointer());
  }

  ITKHelpers::DeepCopy(currentImage.GetPointer(), this->Output.GetPointer());

  ITKHelpers::WriteRGBImage(this->Output.GetPointer(), "ComputeOutput.png");
}

#endif
