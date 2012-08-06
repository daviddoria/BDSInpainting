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

#include <PatchComparison/SSD.h>

#include <PoissonEditing/PoissonEditing.h>

// ITK
#include "itkImageRegionReverseIterator.h"

// STL
#include <ctime>

template <typename TImage>
BDSInpainting<TImage>::BDSInpainting() : Iterations(5),
                                         PatchRadius(7),
                                         PatchMatchFunctor(NULL),
                                         CompositorFunctor(NULL)
{
  this->Output = TImage::New();
  this->Image = TImage::New();
  this->SourceMask = Mask::New();
  this->TargetMask = Mask::New();
}

template <typename TImage>
void BDSInpainting<TImage>::Inpaint()
{
  assert(this->Image->GetLargestPossibleRegion().GetSize()[0] > 0);
  assert(this->SourceMask->GetLargestPossibleRegion().GetSize()[0] > 0);
  assert(this->TargetMask->GetLargestPossibleRegion().GetSize()[0] > 0);

  // Initialize the output with the input
  ITKHelpers::DeepCopy(this->Image.GetPointer(), this->Output.GetPointer());

  // This is done so that the "full region" (the entire images) does not have to be
  // referenced using a particular image or mask. That is, 'SourceMask->GetLargestPossibleRegion()'
  // would not raise the question "Why is this region coming from the SourceMask?"
  itk::ImageRegion<2> fullRegion = this->Image->GetLargestPossibleRegion();

  // Allocate the initial NNField
  typename PatchMatchFunctorType::MatchImageType::Pointer matchImage =
    PatchMatchFunctorType::MatchImageType::New();
  matchImage->SetRegions(this->Image->GetLargestPossibleRegion());
  matchImage->Allocate();

  // Initialize the NNField in the known region
  InitializerKnownRegion<TImage> initializerKnownRegion;
  initializerKnownRegion.SetSourceMask(this->SourceMask);
  initializerKnownRegion.SetImage(this->Image);
  initializerKnownRegion.SetPatchRadius(this->PatchRadius);
  initializerKnownRegion.Initialize(matchImage);

  PatchMatchFunctorType::WriteNNField(matchImage.GetPointer(),
                                      "InitializedKnownRegionNNField.mha"); // debug only

  // Initialize the NNField in the target region
  InitializerRandom<TImage> randomInitializer;
  randomInitializer.SetImage(this->Image);
  randomInitializer.SetTargetMask(this->TargetMask);
  randomInitializer.SetSourceMask(this->SourceMask);
  randomInitializer.SetPatchDistanceFunctor(this->PatchMatchFunctor->GetPatchDistanceFunctor());
  randomInitializer.SetPatchRadius(this->PatchRadius);
  randomInitializer.Initialize(matchImage);

  PatchMatchFunctorType::WriteNNField(matchImage.GetPointer(),
                                      "InitializedRandomNNField.mha"); // debug only
  // Give the PatchMatch functor the data
  this->PatchMatchFunctor->SetImage(this->Image);
  this->PatchMatchFunctor->SetSourceMask(this->SourceMask);
  this->PatchMatchFunctor->SetTargetMask(this->TargetMask);
  this->PatchMatchFunctor->SetInitialNNField(matchImage);
  this->PatchMatchFunctor->Compute();

  PatchMatchFunctorType::WriteNNField(this->PatchMatchFunctor->GetOutput(), "PatchMatchNNField.mha");

  // Update the target pixels
  this->CompositorFunctor->SetTargetMask(this->TargetMask);
  this->CompositorFunctor->SetImage(this->Image);
  this->CompositorFunctor->SetPatchRadius(this->PatchRadius);
  this->CompositorFunctor->SetNearestNeighborField(this->PatchMatchFunctor->GetOutput());
  this->CompositorFunctor->Compute();

  ITKHelpers::DeepCopy(this->CompositorFunctor->GetOutput(), this->Output.GetPointer());

  ITKHelpers::WriteRGBImage(this->Output.GetPointer(), "ComputeOutput.png");
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
void BDSInpainting<TImage>::SetPatchMatchFunctor(PatchMatch<TImage>* patchMatchFunctor)
{
  this->PatchMatchFunctor = patchMatchFunctor;
}

template <typename TImage>
void BDSInpainting<TImage>::SetCompositor(Compositor<TImage>* compositor)
{
  this->CompositorFunctor = compositor;
}

#endif
