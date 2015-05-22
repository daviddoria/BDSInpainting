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
#include "Slots.h"
#include "PixelCompositors.h"

// Submodules
#include <PatchMatch/PropagatorForwardBackward.h>
#include <PatchMatch/RandomSearch.h>
#include <PatchMatch/AcceptanceTestSSD.h>
#include <PatchMatch/AcceptanceTestSourceRegion.h>
#include <PatchMatch/AcceptanceTestComposite.h>
#include <PatchMatch/Process.h>

#include <Helpers/Helpers.h>

template <typename TImage>
BDSInpaintingRings<TImage>::BDSInpaintingRings() : InpaintingAlgorithm<TImage>()
{

}

template <typename TImage>
void BDSInpaintingRings<TImage>::Inpaint()
{
  { // Debug only
  ITKHelpers::WriteImage(this->TargetMask.GetPointer(), "BDS_TargetMask.png");
  ITKHelpers::WriteImage(this->SourceMask.GetPointer(), "BDS_SourceMask.png");
  }

  // Allocate the initial NNField
  this->NNField = PatchMatchHelpers::NNFieldType::New();
  this->NNField->SetRegions(this->Image->GetLargestPossibleRegion());
  this->NNField->Allocate();

  // Initialize the entire NNfield to be empty matches
  MatchSet emptyMatchSet;
  emptyMatchSet.SetMaximumMatches(10);
  ITKHelpers::SetImageToConstant(this->NNField.GetPointer(), emptyMatchSet);

  PatchMatchHelpers::WriteNNField(this->NNField.GetPointer(), "BDS_OriginalInitialized.mha");

  InitializeKnownRegion();

  PatchMatchHelpers::WriteNNField(this->NNField.GetPointer(), "BDS_InitializeKnownRegion.mha");

  ProducePropagationBuffer();

  ITKHelpers::WriteImage(this->SourceMask.GetPointer(), "BDS_SourceMaskWithBuffer.png");

  // Get the region where we need to compute the NNField but not composite
  Mask::Pointer surroundingRingMask = Mask::New();
  GetSurroundingRingMask(surroundingRingMask);

  ITKHelpers::WriteImage(surroundingRingMask.GetPointer(), "BDS_SurroundingRingMask.png");

  // This is the only step that is separate from the ring-at-a-time filling.
  ComputeNNField(surroundingRingMask);

  PatchMatchHelpers::WriteNNField(this->NNField.GetPointer(), "BDS_SurroundingRing.mha");

  // Initialize the output from the original image
  ITKHelpers::DeepCopy(this->Image.GetPointer(), this->Output.GetPointer());

  PatchRadiusThickRings();
}

template <typename TImage>
void BDSInpaintingRings<TImage>::GetSurroundingRingMask(Mask* surroundingMask)
{
  // Expand the hole (ShrinkHole is called because here we mark the "hole"
  // with "valid" pixels.
  Mask::Pointer expandedTargetMask = Mask::New();
  expandedTargetMask->DeepCopyFrom(this->TargetMask);
  expandedTargetMask->ShrinkHole(this->PatchRadius);

  ITKHelpers::WriteImage(expandedTargetMask.GetPointer(), "BDSInpaintingRings_ExpandedTargetMask.png");

  // Get the difference (XOR) between the original hole and the expanded hole
  ITKHelpers::XORImages(expandedTargetMask.GetPointer(), this->TargetMask.GetPointer(),
                        surroundingMask, this->TargetMask->GetValidValue());
  surroundingMask->CopyInformationFrom(this->TargetMask);

  ITKHelpers::WriteImage(surroundingMask, "BDSInpaintingRings_surroundingMask.png");
}

template <typename TImage>
void BDSInpaintingRings<TImage>::InitializeKnownRegion()
{
  // Initialize the NNField in the known region
  InitializerKnownRegion initializerKnownRegion;
  initializerKnownRegion.SetSourceMask(this->SourceMask);
  initializerKnownRegion.SetPatchRadius(this->PatchRadius);
  initializerKnownRegion.Initialize(this->NNField);

  PatchMatchHelpers::WriteNNField(this->NNField.GetPointer(), "BDSInpaintingRings_KnownRegionNNField.mha");
}

template <typename TImage>
void BDSInpaintingRings<TImage>::ProducePropagationBuffer()
{
  // Remove the boundary from the source mask, to give the propagation some breathing room.
  // We just remove a 1 pixel thick boundary around the image, then perform an ExpandHole operation.
  // ExpandHole() only operates on the boundary between Valid and Hole, so if we did not first remove the
  // single pixel boundary nothing would happen to the boundary by the morphological filter.
  // Must do this after the InitializerKnownRegion so that the pixels whose patches are fully in the original
  // source region are initialized properly.
  std::vector<itk::Index<2> > boundaryPixels =
    ITKHelpers::GetBoundaryPixels(this->SourceMask->GetLargestPossibleRegion(), 1);
  ITKHelpers::SetPixels(this->SourceMask.GetPointer(), boundaryPixels, this->SourceMask->GetHoleValue());

  ITKHelpers::WriteImage(this->SourceMask.GetPointer(), "BDSInpaintingRings_BoundaryRemovedSourceMask.png");

  // Shrink the source region, around the border of the image and around the hole.
  this->SourceMask->ExpandHole(this->PatchRadius);
  ITKHelpers::WriteImage(this->SourceMask.GetPointer(), "BDSInpaintingRings_FinalSourceMask.png");
}

template <typename TImage>
void BDSInpaintingRings<TImage>::ConstrainedPatchMatch(Mask* const targetMask,
                                                       const float histogramRatioStart,
                                                       const float histogramRatioStep, const float maxHistogramRatio)
{
  // Create the HSV image
  typedef itk::VectorImage<float, 2> HSVImageType;
  HSVImageType::Pointer hsvImage = HSVImageType::New();
  ITKHelpers::ITKImageToHSVImage(this->Image.GetPointer(), hsvImage.GetPointer());
  ITKHelpers::WriteImage(hsvImage.GetPointer(), "HSV.mha");

  // Setup the patch distance functor
  typedef SSD<TImage> PatchDistanceFunctorType;
  PatchDistanceFunctorType patchDistanceFunctor;
  patchDistanceFunctor.SetImage(this->Image);

  typedef AcceptanceTestSSD AcceptanceTestSSDType;
  AcceptanceTestSSDType acceptanceTestSSD;
  acceptanceTestSSD.SetIncludeInScore(false);

  typedef AcceptanceTestSourceRegion AcceptanceTestSourceRegionType;
  AcceptanceTestSourceRegion acceptanceTestSourceRegion(this->SourceMask);
  acceptanceTestSourceRegion.SetIncludeInScore(false);

  typedef AcceptanceTestNeighborHistogramRatio<HSVImageType> NeighborHistogramRatioAcceptanceTestType;
  NeighborHistogramRatioAcceptanceTestType neighborHistogramRatioAcceptanceTest;
  neighborHistogramRatioAcceptanceTest.SetImage(hsvImage);
  neighborHistogramRatioAcceptanceTest.SetRangeMin(0.0f); // (0,1) is the range of each channel of the HSV image
  neighborHistogramRatioAcceptanceTest.SetRangeMax(1.0f); // (0,1) is the range of each channel of the HSV image
  neighborHistogramRatioAcceptanceTest.SetPatchRadius(this->PatchRadius);
  //neighborHistogramRatioAcceptanceTest.SetNeighborHistogramMultiplier(2.0f); // this will be set in the loop
  neighborHistogramRatioAcceptanceTest.SetNumberOfBinsPerDimension(30);
  neighborHistogramRatioAcceptanceTest.SetIncludeInScore(true);

  typedef AcceptanceTestComposite AcceptanceTestType;
  AcceptanceTestComposite acceptanceTest;
  acceptanceTest.AddAcceptanceTest(&acceptanceTestSourceRegion);
  acceptanceTest.AddAcceptanceTest(&acceptanceTestSSD);
  acceptanceTest.AddAcceptanceTest(&neighborHistogramRatioAcceptanceTest);

  Process* processFunctor = new ProcessValidMaskPixels(targetMask);

  typedef PropagatorForwardBackward<PatchDistanceFunctorType,
          AcceptanceTestType> PropagatorType;
  PropagatorType propagationFunctor;
  propagationFunctor.SetPatchRadius(this->PatchRadius);
  propagationFunctor.SetAcceptanceTest(&acceptanceTest);
  propagationFunctor.SetPatchDistanceFunctor(&patchDistanceFunctor);
  propagationFunctor.SetProcessFunctor(processFunctor);

  WritePatchPair<TImage> propagatedPatchPairWriter(this->Image, this->PatchRadius, "PropagatedPairs");

  auto propagatedPairWriter = [&propagatedPatchPairWriter](const itk::Index<2>& queryCenter, const itk::Index<2>& matchCenter, const float score)
                    {propagatedPatchPairWriter.Write(queryCenter, matchCenter, score);};
//  propagationFunctor.AcceptedSignal.connect(propagatedPairWriter);

  NeighborTestValidMask validMaskNeighborTest(targetMask);

  Neighbors forwardNeighbors;
  forwardNeighbors.SetRegion(this->NNField->GetLargestPossibleRegion());
  NeighborTestForward forwardNeighborTest;
  forwardNeighbors.AddNeighborTest(&forwardNeighborTest);
  forwardNeighbors.AddNeighborTest(&validMaskNeighborTest);
  propagationFunctor.SetForwardNeighborFunctor(&forwardNeighbors);

  Neighbors backwardNeighbors;
  backwardNeighbors.SetRegion(this->NNField->GetLargestPossibleRegion());
  NeighborTestBackward backwardNeighborTest;
  backwardNeighbors.AddNeighborTest(&backwardNeighborTest);
  backwardNeighbors.AddNeighborTest(&validMaskNeighborTest);
  propagationFunctor.SetBackwardNeighborFunctor(&backwardNeighbors);

  typedef RandomSearch<TImage, PatchDistanceFunctorType, AcceptanceTestType>
    RandomSearchType;
  RandomSearchType randomSearcher;
  randomSearcher.SetImage(this->Image);
  randomSearcher.SetPatchRadius(this->PatchRadius);
  randomSearcher.SetSourceMask(this->SourceMask);
  randomSearcher.SetPatchDistanceFunctor(&patchDistanceFunctor);
  randomSearcher.SetProcessFunctor(processFunctor);
  randomSearcher.SetAcceptanceTest(&acceptanceTest);
  randomSearcher.SetRandom(false);

  WritePatchPair<TImage> patchPairWriter(this->Image, this->PatchRadius, "RandomSearchPairs");

  auto pairWriter = [&patchPairWriter](const itk::Index<2>& queryCenter, const itk::Index<2>& matchCenter, const float score)
                    {patchPairWriter.Write(queryCenter, matchCenter, score);};
//  randomSearcher.AcceptedSignal.connect(pairWriter);

  // Initialize the NNField in the PatchRadius thick ring outside of the target region
  InitializerRandom<PatchDistanceFunctorType> initializer;
  initializer.SetPatchDistanceFunctor(&patchDistanceFunctor);
  initializer.SetTargetMask(targetMask);
  initializer.SetSourceMask(this->SourceMask);
  initializer.SetPatchRadius(this->PatchRadius);
  initializer.Initialize(this->NNField);

  PatchMatchHelpers::WriteNNField(this->NNField.GetPointer(), "BDSInpaintingRings_RandomInit.mha");

  // Setup the PatchMatch functor
  PatchMatch patchMatchFunctor;
  unsigned int patchMatchIterations = 4; // This is 2 forward and 2 backward iterations
  patchMatchFunctor.SetIterations(patchMatchIterations);
  //patchMatchFunctor.SetIterations(1);
  //patchMatchFunctor.Compute(nnField, &propagationFunctor, &randomSearcher, processFunctor); // This will be done in the loop

  PatchMatchHelpers::WriteNNField(this->NNField.GetPointer(), "BDSInpaintingRings_BeforeConstrainedPatchMatch.mha");

  float acceptableHistogramRatio = histogramRatioStart;

  unsigned int iteration = 0;

  auto testNoVerifiedMatch = [](const MatchSet& queryMatchSet)
  {
    return !queryMatchSet.HasVerifiedMatch();
  };

  WriteSlot patchMatchWriter("BDS_Phase1");
  auto pmWriter = [&patchMatchWriter](PatchMatchHelpers::NNFieldType* nnField){patchMatchWriter.Write(nnField);};
//  patchMatchFunctor.UpdatedSignal.connect(pmWriter);

  while((PatchMatchHelpers::CountTestedPixels(this->NNField.GetPointer(),
          targetMask, testNoVerifiedMatch) > 0) &&
        (acceptableHistogramRatio <= maxHistogramRatio))
  {
    std::cout << "There are "
              << PatchMatchHelpers::CountTestedPixels(this->NNField.GetPointer(),
                                                      targetMask, testNoVerifiedMatch)
              << " pixels without a verified match remaining." << std::endl;

    neighborHistogramRatioAcceptanceTest.SetMaxNeighborHistogramRatio(acceptableHistogramRatio);

    patchMatchFunctor.Compute(this->NNField, &propagationFunctor, &randomSearcher,
                              processFunctor);

    PatchMatchHelpers::WriteNNField(this->NNField.GetPointer(),
                                    Helpers::GetSequentialFileName("BDSInpaintingRings_PropagatedNNField",
                                                                   iteration, "mha"));
    acceptableHistogramRatio += histogramRatioStep;
    std::cout << "Increased acceptableHistogramRatio to " << acceptableHistogramRatio << std::endl;

    iteration++;
  }
}

template <typename TImage>
void BDSInpaintingRings<TImage>::ForcePropagation(Mask* const targetMask)
{
  std::cout << "Starting force propagation..." << std::endl;

  ProcessUnverifiedValidMaskPixels unverifiedProcessFunctor(this->NNField, targetMask);

  // Create the HSV image
  typedef itk::VectorImage<float, 2> HSVImageType;
  HSVImageType::Pointer hsvImage = HSVImageType::New();
  ITKHelpers::ITKImageToHSVImage(this->Image.GetPointer(), hsvImage.GetPointer());
  ITKHelpers::WriteImage(hsvImage.GetPointer(), "HSV.mha");

  // Setup the patch distance functor
  typedef SSD<TImage> PatchDistanceFunctorType;
  PatchDistanceFunctorType patchDistanceFunctor;
  patchDistanceFunctor.SetImage(this->Image);

  typedef AcceptanceTestSourceRegion AcceptanceTestSourceRegionType;
  AcceptanceTestSourceRegion acceptanceTestSourceRegion(this->SourceMask);
  acceptanceTestSourceRegion.SetIncludeInScore(false);

  // The only acceptance test we want to apply is to make sure the propagated
  // patch is actually valid (completely in the source region)
  typedef PropagatorForwardBackward<PatchDistanceFunctorType,
          AcceptanceTestSourceRegionType> ForcePropagatorType;
  ForcePropagatorType forcePropagator;
  forcePropagator.SetPatchRadius(this->PatchRadius);
  forcePropagator.SetAcceptanceTest(&acceptanceTestSourceRegion);
  forcePropagator.SetPatchDistanceFunctor(&patchDistanceFunctor);
  forcePropagator.SetProcessFunctor(&unverifiedProcessFunctor);

  WritePatchPair<TImage> forcedPropagatedPatchPairWriter(this->Image, this->PatchRadius, "ForcedPropagatedPairs");

  auto forcedPropagatedPairWriter = [&forcedPropagatedPatchPairWriter](const itk::Index<2>& queryCenter, const itk::Index<2>& matchCenter, const float score)
                    {forcedPropagatedPatchPairWriter.Write(queryCenter, matchCenter, score);};
//  forcePropagator.AcceptedSignal.connect(forcedPropagatedPairWriter);

  OutputPixelSlot outputPixelFunctor;
  auto outputPixelSlot = [&outputPixelFunctor](const itk::Index<2>& index)
                    {outputPixelFunctor.OutputPixel(index);};
//  forcePropagator.ProcessPixelSignal.connect(outputPixelSlot);

  NeighborTestVerified verifiedNeighborTest(this->NNField);

  NeighborTestValidMask validMaskNeighborTest(targetMask);

  Neighbors verifiedForwardNeighbors;
  verifiedForwardNeighbors.SetRegion(this->NNField->GetLargestPossibleRegion());
  verifiedForwardNeighbors.AddNeighborTest(&verifiedNeighborTest);
  verifiedForwardNeighbors.AddNeighborTest(&validMaskNeighborTest);
  NeighborTestForward forwardNeighborTest;
  verifiedForwardNeighbors.AddNeighborTest(&forwardNeighborTest);
  forcePropagator.SetForwardNeighborFunctor(&verifiedForwardNeighbors);

  Neighbors verifiedBackwardNeighbors;
  verifiedBackwardNeighbors.SetRegion(this->NNField->GetLargestPossibleRegion());
  verifiedBackwardNeighbors.AddNeighborTest(&verifiedNeighborTest);
  verifiedBackwardNeighbors.AddNeighborTest(&validMaskNeighborTest);
  NeighborTestBackward backwardNeighborTest;
  verifiedBackwardNeighbors.AddNeighborTest(&backwardNeighborTest);
  forcePropagator.SetBackwardNeighborFunctor(&verifiedBackwardNeighbors);

  auto testNoVerifiedMatch = [](const MatchSet& queryMatchSet)
  {
    return !queryMatchSet.HasVerifiedMatch();
  };

  unsigned int iteration = 0;

  while(PatchMatchHelpers::CountTestedPixels(this->NNField.GetPointer(),
          targetMask, testNoVerifiedMatch) > 0)
  {
    std::cout << "Phase 2: There are "
              << PatchMatchHelpers::CountTestedPixels(this->NNField.GetPointer(),
                                                      targetMask, testNoVerifiedMatch)
              << " pixels without a verified match remaining." << std::endl;

    bool force = true;
    unsigned int numberOfPropagatedPixels = forcePropagator.Propagate(this->NNField, force);

    PatchMatchHelpers::WriteNNField(this->NNField.GetPointer(),
                                    Helpers::GetSequentialFileName("BDSRings_NNField_ForceProp",
                                                                   iteration, "mha"));

    // Mark the pixels that were propagated to in the last iteration so that they can be used in the next iteration.
    // Without this, we would only be able to forcefully propagate to 1 pixel away from verified pixels.
    itk::ImageRegionIteratorWithIndex<PatchMatchHelpers::NNFieldType>
       fieldIterator(this->NNField, this->NNField->GetLargestPossibleRegion());

    while(!fieldIterator.IsAtEnd())
    {
      MatchSet matchSet = fieldIterator.Get();
      if(matchSet.GetNumberOfMatches() > 0 && !matchSet.GetMatch(0).GetAllowPropagation())
      {
        matchSet.GetMatch(0).SetAllowPropagation(true);
      }
      ++fieldIterator;
    }

    if(numberOfPropagatedPixels == 0)
    {
      //throw std::runtime_error("Forced propagated zero pixels!");
      std::cout << "Forced propagated zero pixels, stopping." << std::endl;
      break;
    }

    iteration++;
  }
  PatchMatchHelpers::WriteNNField(this->NNField.GetPointer(),
                                  "BDSInpaintingRings_BoundaryNNField.mha");
}

template <typename TImage>
void BDSInpaintingRings<TImage>::ClearUnverifiedPixels()
{

  // Clear all matches for patches that still do not have a good (verified) match
  // This means it remains from the random initialization.
//    auto noVerifiedMatchFunctor = [NNField](const itk::Index<2>& index)
//                                {
//                                  // Don't process pixels that have a verified match
//                                  if(NNField->GetPixel(index).HasVerifiedMatch())
//                                  {
//                                    return false;
//                                  }
//                                  return true;
//                                };
//    auto clearFunctor = [NNField](const itk::Index<2>& index)
//                                {
//                                  MatchSet matchSet = NNField->GetPixel(index);
//                                  matchSet.Clear();
//                                  NNField->SetPixel(index, matchSet);
//                                };

//   ITKHelpers::ApplyOperationToTestedPixels(nnField.GetPointer(),
//                                            noVerifiedMatchFunctor, clearFunctor);

//   PatchMatchHelpers::WriteVerifiedPixels(nnField.GetPointer(), "VerifiedPixels.mha");

}


template <typename TImage>
void BDSInpaintingRings<TImage>::SinglePixelRings()
{
//  // Keep track of which ring we are on
//  unsigned int ringCounter = 0;

//  Mask::Pointer currentTargetMask = Mask::New();
//  currentTargetMask->DeepCopyFrom(this->TargetMask);

//  // Perform ring-at-a-time inpainting
//  while(currentTargetMask->HasValidPixels())
//  {
//    // Get the inside boundary of the target region
//    Mask::BoundaryImageType::Pointer boundaryImage = Mask::BoundaryImageType::New();
//    // In the resulting boundary image, the boundary will be 255.
//    Mask::BoundaryImageType::PixelType boundaryValue = 255;
//    currentTargetMask->FindBoundary(boundaryImage, Mask::VALID, boundaryValue);

//    //ITKHelpers::WriteSequentialImage(boundaryImage.GetPointer(), "BDSRings_Boundary", ringCounter, 4, "png");

//    // Create a mask of just the boundary
//    Mask::Pointer boundaryMask = Mask::New();
//    Mask::BoundaryImageType::PixelType holeValue = 0;
//    Mask::BoundaryImageType::PixelType validValue = boundaryValue;
//    boundaryMask->CreateFromImage(boundaryImage.GetPointer(), holeValue, validValue);
//    //boundaryMask->Invert(); // Make the thin boundary the only valid pixels in the mask

//    //std::cout << "Boundary mask: " << std::endl; boundaryMask->OutputMembers(); // Debug only
//    //ITKHelpers::WriteSequentialImage(boundaryMask.GetPointer(),
//    //     "BDSRings_BoundaryMask", ringCounter, 4, "png");

//    // Set the mask to use in the PatchMatch algorithm
//    currentTargetMask->DeepCopyFrom(boundaryMask);

//    // We set these properties here, but this object is not used here but rather simply
//    // passed along to the composition BDSInpainting object below
//    BDSInpainting<TImage> internalBDSInpaintingFunctor;
//    internalBDSInpaintingFunctor.SetImage(this->Image);
//    internalBDSInpaintingFunctor.SetPatchRadius(this->PatchRadius);
//    internalBDSInpaintingFunctor.SetTargetMask(boundaryMask);
//    internalBDSInpaintingFunctor.SetSourceMask(this->SourceMask);
//    Compositor<TImage, PixelCompositorAverage> compositor;
//    internalBDSInpaintingFunctor.Inpaint(&patchMatchFunctor, &compositor);

//    ITKHelpers::WriteSequentialImage(internalBDSInpaintingFunctor.GetOutput(),
//                                     "BDSRings_InpaintedRing", ringCounter, 4, "png");

//    // Copy the filled ring into the image for the next iteration
//    ITKHelpers::DeepCopy(internalBDSInpaintingFunctor.GetOutput(), currentImage.GetPointer());

//    // Reduce the size of the target region (we "enlarge the hole", because
//    // the "hole" is considered the valid part of the target mask)
//    unsigned int kernelRadius = 1;
//    currentTargetMask->ExpandHole(kernelRadius);

//    ringCounter++;
//  }

//  ITKHelpers::DeepCopy(currentImage.GetPointer(), this->Output.GetPointer());
}

template <typename TImage>
void BDSInpaintingRings<TImage>::FillHole(Mask* const targetMask)
{
  Compositor<TImage, PixelCompositorAverage> compositor;
  compositor.SetImage(this->Output); // We operate on the current intermediate image
  compositor.SetPatchRadius(this->PatchRadius);
  compositor.SetTargetMask(targetMask);
  compositor.SetNearestNeighborField(this->NNField);
  compositor.Composite();

//  ITKHelpers::WriteSequentialImage(compositor.GetOutput(),
//                                   "BDSRings_InpaintedRing", ringCounter, 4, "png");

  // Copy the filled ring into the image for the next iteration
  ITKHelpers::DeepCopy(compositor.GetOutput(), this->Output.GetPointer());

}

template <typename TImage>
void BDSInpaintingRings<TImage>::PatchRadiusThickRings()
{
  std::cout << "PatchRadiusThickRings()" << std::endl;

  // Keep track of which ring we are on
  unsigned int ringCounter = 0;

  Mask::Pointer remainingHoleMask = Mask::New();
  remainingHoleMask->DeepCopyFrom(this->TargetMask);

  // Perform patch-radius-thick-ring-at-a-time inpainting
  while(remainingHoleMask->HasValidPixels())
  {
    // Get the inside boundary of the target region
    Mask::Pointer shrunkHole = Mask::New();
    shrunkHole->DeepCopyFrom(remainingHoleMask);
    shrunkHole->ExpandHole(this->PatchRadius);

    Mask::Pointer ringMask = Mask::New();
    // Get the difference (XOR) between the original hole and the expanded hole
    ITKHelpers::XORImages(remainingHoleMask.GetPointer(), shrunkHole.GetPointer(),
                          ringMask.GetPointer(), this->TargetMask->GetValidValue());
    ringMask->CopyInformationFrom(this->TargetMask);

    ITKHelpers::WriteSequentialImage(ringMask.GetPointer(), "BDS_RingMask", ringCounter, 3, "png");

    ComputeNNField(ringMask);

    FillHole(ringMask);

    remainingHoleMask->DeepCopyFrom(shrunkHole);

    ITKHelpers::WriteSequentialImage(this->Output.GetPointer(), "BDS_PatchRadiusThickRings", ringCounter, 3, "png");
    PatchMatchHelpers::WriteNNField(this->NNField.GetPointer(), Helpers::GetSequentialFileName("PatchRadiusThickRings", ringCounter, "mha, 3"));
    //ITKHelpers::WriteImage(this->Output.GetPointer(), "BDS_PatchRadiusThickRings.png");
    ringCounter++;
  }
}

template <typename TImage>
void BDSInpaintingRings<TImage>::ComputeNNField(Mask* const targetMaskInput)
{
  // Copy the mask, because we will potentially modify it in this function,
  // but it needs to remain unchanged in the calling function

  Mask::Pointer targetMask = Mask::New();
  targetMask->DeepCopyFrom(targetMaskInput);

  unsigned int numberOfUnverifiedPixels =
      PatchMatchHelpers::CountUnverifiedPixels(this->NNField.GetPointer(), targetMask.GetPointer());

  assert(targetMask->CountValidPixels() == numberOfUnverifiedPixels);

  float histogramRatioStart = 2.0f;
  float histogramRatioStep = 0.2f;
  //float maxHistogramRatio = 5.0f;
  float maxHistogramRatio = 3.2f;

  unsigned int iteration = 0;
  do
  {
    ConstrainedPatchMatch(targetMask, histogramRatioStart, histogramRatioStep, maxHistogramRatio);

    ForcePropagation(targetMask);

    numberOfUnverifiedPixels =
          PatchMatchHelpers::CountUnverifiedPixels(this->NNField.GetPointer(), targetMask.GetPointer());
    std::cout << "After iteration " << iteration << " of ComputeNNField(), there are "
              << numberOfUnverifiedPixels << " numberOfUnverifiedPixels." << std::endl;
    maxHistogramRatio += 1.0f;
    std::cout << "Increased maxHistogramRatio to " << maxHistogramRatio << std::endl;

    // Create a targetMask of only the pixels which still remain to be propagated
    std::vector<itk::Index<2> > remainingPixels = PatchMatchHelpers::GetUnverifiedPixels(this->NNField.GetPointer(), targetMask.GetPointer());
    ITKHelpers::SetImageToConstant(targetMask.GetPointer(), targetMask->GetHoleValue());
    ITKHelpers::SetPixels(targetMask.GetPointer(), remainingPixels, targetMask->GetValidValue());

    PatchMatchHelpers::WriteNNField(this->NNField.GetPointer(), Helpers::GetSequentialFileName("BDS_ComputeNNField", iteration, "mha", 3));

    iteration++;

  } while(numberOfUnverifiedPixels > 0);

}

#endif
