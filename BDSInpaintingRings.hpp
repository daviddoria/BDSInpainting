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
#include "BDSInpainting.h" // Composition

// Custom
#include "InitializerRandom.h"
#include "InitializerKnownRegion.h"
#include "InitializerNeighborHistogram.h"
#include "AcceptanceTestNeighborHistogram.h"

template <typename TImage>
BDSInpaintingRings<TImage>::BDSInpaintingRings() : InpaintingAlgorithm<TImage>()
{

}

template <typename TImage>
void BDSInpaintingRings<TImage>::Inpaint()
{
  { // Debug only
  ITKHelpers::WriteImage(this->TargetMask.GetPointer(), "BDSInpaintingRings_TargetMask.png");
  ITKHelpers::WriteImage(this->SourceMask.GetPointer(), "BDSInpaintingRings_SourceMask.png");
  }

  // Save the original mask, as we will be modifying the internal masks below
  Mask::Pointer currentTargetMask = Mask::New();
  currentTargetMask->DeepCopyFrom(this->TargetMask);

  // Initialize the image from the original image
  typename TImage::Pointer currentImage = TImage::New();
  ITKHelpers::DeepCopy(this->Image.GetPointer(), currentImage.GetPointer());

  // The pixels from which information is allowed to propagate (everywhere except the target region)
  // (This is computed again at each iteration through the loop)
  Mask::Pointer currentPropagationMask = Mask::New();
  // We trust the information everywhere except in the hole
  currentPropagationMask->DeepCopyFrom(currentTargetMask);
  currentPropagationMask->InvertData();

  ITKHelpers::WriteImage(currentPropagationMask.GetPointer(), "BDSInpaintingRings_InitialPropagationMask.png");

  // Compute the NN-field in the PatchRadius thick region around the target region. This region
  // does not have a trivial (exactly itself) NNField because each patch centered on one of these
  // pixels has some pixels that are in the target region.
  Mask::Pointer expandedTargetMask = Mask::New();
  expandedTargetMask->DeepCopyFrom(this->TargetMask);
  expandedTargetMask->ShrinkHole(this->PatchRadius);
  ITKHelpers::WriteImage(expandedTargetMask.GetPointer(), "BDSInpaintingRings_ExpandedTargetMask.png");

  Mask::Pointer outsideTargetMask = Mask::New();
  ITKHelpers::XORImages(expandedTargetMask.GetPointer(), currentTargetMask.GetPointer(), outsideTargetMask.GetPointer(), this->TargetMask->GetValidValue());
  outsideTargetMask->CopyInformationFrom(this->TargetMask);

  ITKHelpers::WriteImage(outsideTargetMask.GetPointer(), "BDSInpaintingRings_OutsideTargetMask.png");

  // Allocate the initial NNField
  typename PatchMatchFunctorType::MatchImageType::Pointer nnField =
    PatchMatchFunctorType::MatchImageType::New();
  nnField->SetRegions(currentImage->GetLargestPossibleRegion());
  nnField->Allocate();

  // Initialize the NNField in the known region
  InitializerKnownRegion initializerKnownRegion;
  initializerKnownRegion.SetSourceMask(this->SourceMask);
  initializerKnownRegion.SetPatchRadius(this->PatchRadius);
  initializerKnownRegion.Initialize(nnField);

  PatchMatchFunctorType::WriteNNField(nnField.GetPointer(), "BDSInpaintingRings_KnownRegionNNField.mha");

  // Remove the boundary from the source mask, to give the propagation some breathing room.
  // We just remove a 1 pixel thick boundary around the image, then perform an ExpandHole operation.
  // ExpandHole() only operates on the boundary between Valid and Hole, so if we did not first remove the
  // single pixel boundary nothing would happen to the boundary by the morphological filter.
  // Must do this after the InitializerKnownRegion so that the pixels whose patches are fully in the original
  // source region are initialized properly.
  std::vector<itk::Index<2> > boundaryPixels =
    //ITKHelpers::GetBoundaryPixels(this->SourceMask->GetLargestPossibleRegion(), this->PatchRadius);
    ITKHelpers::GetBoundaryPixels(this->SourceMask->GetLargestPossibleRegion(), 1);
  ITKHelpers::SetPixels(this->SourceMask.GetPointer(), boundaryPixels, this->SourceMask->GetHoleValue());

  ITKHelpers::WriteImage(this->SourceMask.GetPointer(), "BDSInpaintingRings_BoundaryRemovedSourceMask.png");

  this->SourceMask->ExpandHole(this->PatchRadius);
  ITKHelpers::WriteImage(this->SourceMask.GetPointer(), "BDSInpaintingRings_FinalSourceMask.png");

  // Create the HSV image
  typedef itk::VectorImage<float, 2> HSVImageType;
  HSVImageType::Pointer hsvImage = HSVImageType::New();
  ITKHelpers::ITKImageToHSVImage(currentImage.GetPointer(), hsvImage.GetPointer());
  ITKHelpers::WriteImage(hsvImage.GetPointer(), "HSV.mha");

//   ITKHelpers::ScaleAllChannelsTo255(hsvImage.GetPointer());
//   typename TImage::Pointer castedHSVImage = TImage::New();
//   ITKHelpers::CastImage(hsvImage.GetPointer(), castedHSVImage.GetPointer());
//   ITKHelpers::WriteImage(castedHSVImage.GetPointer(), "CastedHSV.mha");

  // Initialize the NNField in the PatchRadius thick ring outside of the target region
  //InitializerRandom<TImage> initializer;

  typedef SSD<TImage> SSDFunctorType;
  SSDFunctorType ssdFunctor;
  ssdFunctor.SetImage(this->Image);

  InitializerNeighborHistogram<HSVImageType, SSDFunctorType> initializer;
  initializer.SetNeighborHistogramMultiplier(2.0f);
  initializer.SetPatchDistanceFunctor(&ssdFunctor);
  initializer.SetRangeMin(0.0f);
  initializer.SetRangeMax(1.0f);
  initializer.SetImage(hsvImage);
  initializer.SetTargetMask(outsideTargetMask);
  initializer.SetSourceMask(this->SourceMask);
  initializer.SetPatchRadius(this->PatchRadius);
  initializer.Initialize(nnField);

  PatchMatchFunctorType::WriteNNField(nnField.GetPointer(), "BDSInpaintingRings_InitializedNNField.mha");

  // Setup acceptance test
  AcceptanceTestNeighborHistogram<HSVImageType> acceptanceTest;
  acceptanceTest.SetNeighborHistogramMultiplier(2.0f);
  acceptanceTest.SetImage(hsvImage);
  acceptanceTest.SetRangeMin(0);
  acceptanceTest.SetRangeMax(255);
  this->PatchMatchFunctor->SetAcceptanceTest(&acceptanceTest);

  this->PatchMatchFunctor->SetImage(currentImage);
  this->PatchMatchFunctor->SetAllowedPropagationMask(currentPropagationMask);
  this->PatchMatchFunctor->SetPropagationStrategy(PatchMatchFunctorType::RASTER);
  this->PatchMatchFunctor->SetTargetMask(outsideTargetMask);
  this->PatchMatchFunctor->SetSourceMask(this->SourceMask);
  this->PatchMatchFunctor->SetInitialNNField(nnField);
  this->PatchMatchFunctor->Compute();

  PatchMatchFunctorType::WriteNNField(this->PatchMatchFunctor->GetOutput(), "BDSInpaintingRings_PatchMatchNNField.mha");

  // Keep track of which ring we are on
  unsigned int ringCounter = 0;

  typename PatchMatchFunctorType::MatchImageType::Pointer previousNNField = PatchMatchFunctorType::MatchImageType::New();
  ITKHelpers::DeepCopy(previousNNField.GetPointer(), this->PatchMatchFunctor->GetOutput());

  // Perform ring-at-a-time inpainting
  while(currentTargetMask->HasValidPixels())
  {
    // We trust the information everywhere except in the hole
    currentPropagationMask->DeepCopyFrom(currentTargetMask);
    currentPropagationMask->InvertData();

    //ITKHelpers::WriteSequentialImage(currentPropagationMask.GetPointer(),
//     "BDSRings_CurrentPropagationMask", ringCounter, 4, "png");

    // Get the inside boundary of the target region
    Mask::BoundaryImageType::Pointer boundaryImage = Mask::BoundaryImageType::New();
    // In the resulting boundary image, the boundary will be 255.
    Mask::BoundaryImageType::PixelType boundaryValue = 255; 
    currentTargetMask->FindBoundary(boundaryImage, Mask::VALID, boundaryValue);

    //ITKHelpers::WriteSequentialImage(boundaryImage.GetPointer(), "BDSRings_Boundary", ringCounter, 4, "png");

    // Create a mask of just the boundary
    Mask::Pointer boundaryMask = Mask::New();
    Mask::BoundaryImageType::PixelType holeValue = 0;
    Mask::BoundaryImageType::PixelType validValue = boundaryValue;
    boundaryMask->CreateFromImage(boundaryImage.GetPointer(), holeValue, validValue);
    //boundaryMask->Invert(); // Make the thin boundary the only valid pixels in the mask

    //std::cout << "Boundary mask: " << std::endl; boundaryMask->OutputMembers(); // Debug only
    //ITKHelpers::WriteSequentialImage(boundaryMask.GetPointer(),
//     "BDSRings_BoundaryMask", ringCounter, 4, "png");

    // Set the mask to use in the PatchMatch algorithm
    currentTargetMask->DeepCopyFrom(boundaryMask);

    // We set these properties here, but this object is not used here but rather simply passed along to the composition BDSInpainting object below
    this->PatchMatchFunctor->SetAllowedPropagationMask(currentPropagationMask);
    this->PatchMatchFunctor->SetPropagationStrategy(PatchMatchFunctorType::INWARD);
    //this->PatchMatchFunctor->GetAcceptanceTest()->SetImage(currentImage);
    if(dynamic_cast<AcceptanceTestImage<TImage>*>(this->PatchMatchFunctor->GetAcceptanceTest()))
    {
      dynamic_cast<AcceptanceTestImage<TImage>*>(this->PatchMatchFunctor->GetAcceptanceTest())->SetImage(currentImage);
    }
    else
    {
      throw std::runtime_error("Dynamic cast failed in BDSInpaintingRings::Inpaint()!");
    }
    this->PatchMatchFunctor->SetTargetMask(currentTargetMask);
    this->PatchMatchFunctor->SetAllowedPropagationMask(currentPropagationMask);
    if(previousNNField)
    {
      this->PatchMatchFunctor->SetInitialNNField(previousNNField);
    }

    BDSInpainting<TImage> internalBDSInpaintingFunctor;
    internalBDSInpaintingFunctor.SetImage(this->Image);
    internalBDSInpaintingFunctor.SetPatchRadius(this->PatchRadius);
    internalBDSInpaintingFunctor.SetTargetMask(boundaryMask);
    internalBDSInpaintingFunctor.SetSourceMask(this->SourceMask);
    internalBDSInpaintingFunctor.SetPatchMatchFunctor(this->PatchMatchFunctor);
    internalBDSInpaintingFunctor.Inpaint();

    ITKHelpers::WriteSequentialImage(internalBDSInpaintingFunctor.GetOutput(),
                                     "BDSRings_InpaintedRing", ringCounter, 4, "png");

    // Copy the filled ring into the image for the next iteration
    ITKHelpers::DeepCopy(internalBDSInpaintingFunctor.GetOutput(), currentImage.GetPointer());

    // Reduce the size of the target region (we "enlarge the hole", because
    // the "hole" is considered the valid part of the target mask)
    unsigned int kernelRadius = 1;
    currentTargetMask->ExpandHole(kernelRadius);

    if(!previousNNField)
    {
      previousNNField = PatchMatchFunctorType::MatchImageType::New();
    }

    ITKHelpers::DeepCopy(this->PatchMatchFunctor->GetOutput(), previousNNField.GetPointer());
    this->PatchMatchFunctor->SetAllowedPropagationMask(currentTargetMask);
    ringCounter++;
  }

  ITKHelpers::DeepCopy(currentImage.GetPointer(), this->Output.GetPointer());

}

#endif
