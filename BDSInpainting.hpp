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
#include <PatchMatch/InitializerKnownRegion.h>
#include <PatchMatch/InitializerRandom.h>
#include <PatchMatch/AcceptanceTestNeighborHistogram.h>
#include <PatchMatch/PropagatorForwardBackward.h>

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
  PatchMatchHelpers::NNFieldType::Pointer nnField =
    PatchMatchHelpers::NNFieldType::New();
  nnField->SetRegions(currentImage->GetLargestPossibleRegion());
  nnField->Allocate();

  // Initialize the NNField in the target region
  typedef SSD<TImage> PatchDistanceFunctorType;
  PatchDistanceFunctorType patchDistanceFunctor;
  patchDistanceFunctor.SetImage(currentImage);

  // Initialize the NNField in the known region
  InitializerKnownRegion initializerKnownRegion;
  initializerKnownRegion.SetSourceMask(this->SourceMask);
  initializerKnownRegion.SetPatchRadius(this->PatchRadius);
  initializerKnownRegion.Initialize(nnField);


  //////////////////// The things in this block should probably be passed into this class
  // Setup the patch distance functor
  typedef SSD<TImage> SSDFunctorType;
  SSDFunctorType ssdFunctor;
  ssdFunctor.SetImage(this->Image);

  // Set acceptance test to histogram threshold
  // Create the HSV image
  typedef itk::VectorImage<float, 2> HSVImageType;
  HSVImageType::Pointer hsvImage = HSVImageType::New();
  ITKHelpers::ITKImageToHSVImage(currentImage.GetPointer(), hsvImage.GetPointer());
  ITKHelpers::WriteImage(hsvImage.GetPointer(), "HSV.mha");

  typedef AcceptanceTestNeighborHistogram<HSVImageType> AcceptanceTestType;
  AcceptanceTestType acceptanceTest;
  acceptanceTest.SetImage(hsvImage);
  acceptanceTest.SetRangeMin(0.0f);
  acceptanceTest.SetRangeMax(1.0f);
  acceptanceTest.SetPatchRadius(this->PatchRadius);
  acceptanceTest.SetNeighborHistogramMultiplier(2.0f);

  typedef ProcessValidMaskPixels ProcessFunctorType;
  ProcessFunctorType processFunctor(this->TargetMask);

  typedef PropagatorForwardBackward PropagatorType;
  PropagatorType propagationFunctor;
  propagationFunctor.SetPatchRadius(this->PatchRadius);
  propagationFunctor.Propagate(nnField, &ssdFunctor, &processFunctor, &acceptanceTest);

  typedef RandomSearch<TImage> RandomSearchType;
  RandomSearchType randomSearcher;
  randomSearcher.SetImage(this->Image);

  // Setup the PatchMatch functor. Use a generic (parent class) AcceptanceTest.
  patchMatchFunctor->SetPatchRadius(this->PatchRadius);
  patchMatchFunctor->SetIterations(5);

  for(unsigned int iteration = 0; iteration < this->Iterations; ++iteration)
  {
    PatchMatchHelpers::WriteNNField(nnField.GetPointer(),
                                        "InitializedKnownRegionNNField.mha"); // debug only

    InitializerRandom<PatchDistanceFunctorType> randomInitializer;
    randomInitializer.SetTargetMask(this->TargetMask);
    randomInitializer.SetSourceMask(this->SourceMask);
    randomInitializer.SetPatchDistanceFunctor(&patchDistanceFunctor);
    randomInitializer.SetPatchRadius(this->PatchRadius);
    randomInitializer.Initialize(nnField);

    PatchMatchHelpers::WriteNNField(nnField.GetPointer(),
                                        "InitializedRandomNNField.mha"); // debug only
    // Give the PatchMatch functor the data
    patchMatchFunctor->SetSourceMask(this->SourceMask);
    patchMatchFunctor->SetTargetMask(this->TargetMask);
    patchMatchFunctor->Compute(nnField, &propagationFunctor, &randomSearcher);

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
