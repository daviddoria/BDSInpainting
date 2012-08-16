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
#include "Verifier.h"
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
  ITKHelpers::WriteImage(this->TargetMask.GetPointer(), "BDSInpaintingRings_TargetMask.png");
  ITKHelpers::WriteImage(this->SourceMask.GetPointer(), "BDSInpaintingRings_SourceMask.png");
  }

  // Save the original mask, as we will be modifying the internal masks below
  Mask::Pointer currentTargetMask = Mask::New();
  currentTargetMask->DeepCopyFrom(this->TargetMask);

  // Initialize the image from the original image
  typename TImage::Pointer currentImage = TImage::New();
  ITKHelpers::DeepCopy(this->Image.GetPointer(), currentImage.GetPointer());

  // We first need to compute the NN-field in the PatchRadius thick region around the target region.
  // This region does not have a trivial (exactly itself) NNField because each patch centered on one of these
  // pixels has some pixels that are in the target region.
  Mask::Pointer expandedTargetMask = Mask::New();
  expandedTargetMask->DeepCopyFrom(this->TargetMask);
  expandedTargetMask->ShrinkHole(this->PatchRadius);
  ITKHelpers::WriteImage(expandedTargetMask.GetPointer(), "BDSInpaintingRings_ExpandedTargetMask.png");

  Mask::Pointer outsideTargetMask = Mask::New();
  ITKHelpers::XORImages(expandedTargetMask.GetPointer(), currentTargetMask.GetPointer(),
                        outsideTargetMask.GetPointer(), this->TargetMask->GetValidValue());
  outsideTargetMask->CopyInformationFrom(this->TargetMask);

  ITKHelpers::WriteImage(outsideTargetMask.GetPointer(), "BDSInpaintingRings_OutsideTargetMask.png");

  // Allocate the initial NNField
  typename PatchMatchHelpers::NNFieldType::Pointer nnField =
    PatchMatchHelpers::NNFieldType::New();
  nnField->SetRegions(currentImage->GetLargestPossibleRegion());
  nnField->Allocate();

  // Initialize the entire NNfield to be empty matches
  MatchSet emptyMatchSet;
  emptyMatchSet.SetMaximumMatches(10);
  ITKHelpers::SetImageToConstant(nnField.GetPointer(), emptyMatchSet);

  PatchMatchHelpers::WriteNNField(nnField.GetPointer(), "BDSInpaintingRings_OriginalInitialized.mha");

  // Initialize the NNField in the known region
  InitializerKnownRegion initializerKnownRegion;
  initializerKnownRegion.SetSourceMask(this->SourceMask);
  initializerKnownRegion.SetPatchRadius(this->PatchRadius);
  initializerKnownRegion.Initialize(nnField);

  PatchMatchHelpers::WriteNNField(nnField.GetPointer(), "BDSInpaintingRings_KnownRegionNNField.mha");

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

  this->SourceMask->ExpandHole(this->PatchRadius);
  ITKHelpers::WriteImage(this->SourceMask.GetPointer(), "BDSInpaintingRings_FinalSourceMask.png");

  // Create the HSV image
  typedef itk::VectorImage<float, 2> HSVImageType;
  HSVImageType::Pointer hsvImage = HSVImageType::New();
  ITKHelpers::ITKImageToHSVImage(currentImage.GetPointer(), hsvImage.GetPointer());
  ITKHelpers::WriteImage(hsvImage.GetPointer(), "HSV.mha");

  //////////////////// The things in this block should probably be passed into this class
  // Setup the patch distance functor
  typedef SSD<TImage> PatchDistanceFunctorType;
  PatchDistanceFunctorType patchDistanceFunctor;
  patchDistanceFunctor.SetImage(this->Image);

  typedef AcceptanceTestSSD AcceptanceTestSSDType;
  AcceptanceTestSSDType acceptanceTestSSD;
  acceptanceTestSSD.SetIncludeInScore(false);

  typedef AcceptanceTestSourceRegion AcceptanceTestSourceRegionType;
  AcceptanceTestSourceRegion acceptanceTestSourceRegion(this->SourceMask);

  typedef AcceptanceTestNeighborHistogramRatio<HSVImageType> NeighborHistogramRatioAcceptanceTestType;
  NeighborHistogramRatioAcceptanceTestType neighborHistogramRatioAcceptanceTest;
  neighborHistogramRatioAcceptanceTest.SetImage(hsvImage);
  neighborHistogramRatioAcceptanceTest.SetRangeMin(0.0f);
  neighborHistogramRatioAcceptanceTest.SetRangeMax(1.0f);
  neighborHistogramRatioAcceptanceTest.SetPatchRadius(this->PatchRadius);
  //neighborHistogramRatioAcceptanceTest.SetNeighborHistogramMultiplier(2.0f); // this will be set in the loop
  neighborHistogramRatioAcceptanceTest.SetNumberOfBinsPerDimension(30);

  typedef AcceptanceTestComposite AcceptanceTestType;
  AcceptanceTestComposite acceptanceTest;
  acceptanceTest.AddAcceptanceTest(&acceptanceTestSourceRegion);
  acceptanceTest.AddAcceptanceTest(&acceptanceTestSSD);
  acceptanceTest.AddAcceptanceTest(&neighborHistogramRatioAcceptanceTest);

  Process* processFunctor = new ProcessValidMaskPixels(outsideTargetMask);

  typedef PropagatorForwardBackward<PatchDistanceFunctorType,
          AcceptanceTestType> PropagatorType;
  PropagatorType propagationFunctor;
  propagationFunctor.SetPatchRadius(this->PatchRadius);
  propagationFunctor.SetAcceptanceTest(&acceptanceTest);
  propagationFunctor.SetPatchDistanceFunctor(&patchDistanceFunctor);
  propagationFunctor.SetProcessFunctor(processFunctor);

  NeighborTestValidMask validMaskNeighborTest(outsideTargetMask);

  Neighbors forwardNeighbors;
  forwardNeighbors.SetRegion(nnField->GetLargestPossibleRegion());
  NeighborTestForward forwardNeighborTest;
  forwardNeighbors.AddNeighborTest(&forwardNeighborTest);
  forwardNeighbors.AddNeighborTest(&validMaskNeighborTest);
  propagationFunctor.SetForwardNeighborFunctor(&forwardNeighbors);

  Neighbors backwardNeighbors;
  backwardNeighbors.SetRegion(nnField->GetLargestPossibleRegion());
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

  // Initialize the NNField in the PatchRadius thick ring outside of the target region
  InitializerRandom<PatchDistanceFunctorType> initializer;
  initializer.SetPatchDistanceFunctor(&patchDistanceFunctor);
  initializer.SetTargetMask(outsideTargetMask);
  initializer.SetSourceMask(this->SourceMask);
  initializer.SetPatchRadius(this->PatchRadius);
  initializer.Initialize(nnField);

  PatchMatchHelpers::WriteNNField(nnField.GetPointer(), "BDSInpaintingRings_RandomInit.mha");

  // Setup the PatchMatch functor
  PatchMatch patchMatchFunctor;
  //patchMatchFunctor.SetIterations(3);
  patchMatchFunctor.SetIterations(1);
  //patchMatchFunctor.Compute(nnField, &propagationFunctor, &randomSearcher, processFunctor); // This will be done in the loop

  PatchMatchHelpers::WriteNNField(nnField.GetPointer(), "BDSInpaintingRings_FirstPatchMatch.mha");

  float histogramMultiplierInitial = 2.0f;
  //float histogramMultiplierStep = 0.2f;
  float histogramMultiplierStep = 10;

//   float histogramMultiplierInitial = 3.0f;
//   float histogramMultiplierStep = 1.0f;
  float histogramMultiplier = histogramMultiplierInitial;

  unsigned int iteration = 0;
  auto testNoVerifiedMatch = [](const MatchSet& queryMatchSet)
  {
    return !queryMatchSet.HasVerifiedMatch();
  };

  while(PatchMatchHelpers::CountTestedPixels(nnField.GetPointer(),
          outsideTargetMask, testNoVerifiedMatch) > 0)
  {
    std::cout << "There are "
              << PatchMatchHelpers::CountTestedPixels(nnField.GetPointer(),
                                                      outsideTargetMask, testNoVerifiedMatch)
              << " pixels without a verified match remaining." << std::endl;

    neighborHistogramRatioAcceptanceTest.SetMaxNeighborHistogramRatio(histogramMultiplier);

    patchMatchFunctor.Compute(nnField, &propagationFunctor, &randomSearcher,
                              processFunctor);

    PatchMatchHelpers::WriteNNField(nnField.GetPointer(),
                                    Helpers::GetSequentialFileName("BDSInpaintingRings_PropagatedNNField",
                                                                   iteration, "mha"));
    histogramMultiplier += histogramMultiplierStep;
    std::cout << "Increased histogramMultiplier to " << histogramMultiplier << std::endl;

    if(histogramMultiplier > 5)
    {
      break;
    }

    iteration++;
  }

  // Clear all matches for patches that still do not have a good (verified) match
  // This means it remains from the random initialization.
  auto noVerifiedMatchFunctor = [nnField](const itk::Index<2>& index)
                              {
                                // Don't process pixels that have a verified match
                                if(nnField->GetPixel(index).HasVerifiedMatch())
                                {
                                  return false;
                                }
                                return true;
                              };
  auto clearFunctor = [nnField](const itk::Index<2>& index)
                              {
                                MatchSet matchSet = nnField->GetPixel(index);
                                matchSet.Clear();
                                nnField->SetPixel(index, matchSet);
                              };

//   ITKHelpers::ApplyOperationToTestedPixels(nnField.GetPointer(),
//                                            noVerifiedMatchFunctor, clearFunctor);

//   PatchMatchHelpers::WriteVerifiedPixels(nnField.GetPointer(), "VerifiedPixels.mha");

  histogramMultiplier = histogramMultiplierInitial;

  auto outputWhichFailed = [](const std::string& whichFailed)
                 {
                   std::cout << whichFailed << std::endl;
                };
  //acceptanceTest.WhichFailedSignal.connect(outputWhichFailed);

  auto outputFailedScore = [](const float failedScore)
                 {
                   std::cout << "Acceptance test failed with score " << failedScore << std::endl;
                };
  //acceptanceTest.FailedScoreSignal.connect(outputFailedScore);

  std::cout << "Starting propagation only phase." << std::endl;

  ProcessUnverifiedValidMaskPixels unverifiedProcessFunctor(nnField, outsideTargetMask);

  // The only acceptance test we want to apply is to make sure the propagated
  // patch is actually valid (completely in the source region)
  typedef PropagatorForwardBackward<PatchDistanceFunctorType,
          AcceptanceTestSourceRegionType> ForcePropagatorType;
  ForcePropagatorType forcePropagator;
  forcePropagator.SetPatchRadius(this->PatchRadius);
  forcePropagator.SetAcceptanceTest(&acceptanceTestSourceRegion);
  forcePropagator.SetPatchDistanceFunctor(&patchDistanceFunctor);
  forcePropagator.SetProcessFunctor(&unverifiedProcessFunctor);

  NeighborTestVerified verifiedNeighborTest(nnField);

  Neighbors verifiedForwardNeighbors;
  verifiedForwardNeighbors.SetRegion(nnField->GetLargestPossibleRegion());
  verifiedForwardNeighbors.AddNeighborTest(&verifiedNeighborTest);
  NeighborTestForward forceForwardNeighborTest; // This is only named 'force*' because 'forwardNeighborTest was already declared for use in the previous (non-forced) propagator
  verifiedForwardNeighbors.AddNeighborTest(&forceForwardNeighborTest);
  forcePropagator.SetForwardNeighborFunctor(&verifiedForwardNeighbors);

  Neighbors verifiedBackwardNeighbors;
  verifiedBackwardNeighbors.SetRegion(nnField->GetLargestPossibleRegion());
  verifiedBackwardNeighbors.AddNeighborTest(&verifiedNeighborTest);
  NeighborTestBackward forceBackwardNeighborTest; // This is only named 'force*' because 'backwardNeighborTest was already declared for use in the previous (non-forced) propagator
  verifiedBackwardNeighbors.AddNeighborTest(&forceBackwardNeighborTest);
  forcePropagator.SetBackwardNeighborFunctor(&verifiedBackwardNeighbors);

  int propCounter = 0;
  auto writePropagatedSlot = [&propCounter](PatchMatchHelpers::NNFieldType* nnField)
                 {
                    PatchMatchHelpers::WriteNNField(nnField,
                                Helpers::GetSequentialFileName("PropOnly",
                                                                propCounter, "mha"));
                    propCounter++;
                 };
  //forcePropagator.PropagatedSignal.connect(writePropagatedSlot);

  iteration = 0;
  while(PatchMatchHelpers::CountTestedPixels(nnField.GetPointer(),
          outsideTargetMask, testNoVerifiedMatch) > 0)
  {
    std::cout << "Phase 2: There are "
              << PatchMatchHelpers::CountTestedPixels(nnField.GetPointer(),
                                                      outsideTargetMask, testNoVerifiedMatch)
              << " pixels without a verified match remaining." << std::endl;

    bool force = true;
    forcePropagator.Propagate(nnField, force);

    PatchMatchHelpers::WriteNNField(nnField.GetPointer(),
                                    Helpers::GetSequentialFileName("BDSRings_NNField_PropOnly",
                                                                   iteration, "mha"));

    // Mark the pixels that were propagated to in the last iteration so that they can be used in the next iteration.
    // Without this, we would only be able to forcefully propagate to 1 pixel away from verified pixels.
    itk::ImageRegionIteratorWithIndex<PatchMatchHelpers::NNFieldType> fieldIterator(nnField, nnField->GetLargestPossibleRegion());
    while(!fieldIterator.IsAtEnd())
    {
      MatchSet matchSet = fieldIterator.Get();
      if(matchSet.GetNumberOfMatches() > 0 && !matchSet.GetMatch(0).GetAllowPropagation())
      {
        matchSet.GetMatch(0).SetAllowPropagation(true);
      }
      ++fieldIterator;
    }

    //break; // TODO: Remove this
    iteration++;
  }
  PatchMatchHelpers::WriteNNField(nnField.GetPointer(),
                                  "BDSInpaintingRings_BoundaryNNField.mha");
  exit(-1); // TODO: remove this 
  // Keep track of which ring we are on
  unsigned int ringCounter = 0;

  // Perform ring-at-a-time inpainting
  while(currentTargetMask->HasValidPixels())
  {
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

    // We set these properties here, but this object is not used here but rather simply
    // passed along to the composition BDSInpainting object below
    BDSInpainting<TImage> internalBDSInpaintingFunctor;
    internalBDSInpaintingFunctor.SetImage(this->Image);
    internalBDSInpaintingFunctor.SetPatchRadius(this->PatchRadius);
    internalBDSInpaintingFunctor.SetTargetMask(boundaryMask);
    internalBDSInpaintingFunctor.SetSourceMask(this->SourceMask);
    Compositor<TImage, PixelCompositorAverage> compositor;
    internalBDSInpaintingFunctor.Inpaint(&patchMatchFunctor, &compositor);

    ITKHelpers::WriteSequentialImage(internalBDSInpaintingFunctor.GetOutput(),
                                     "BDSRings_InpaintedRing", ringCounter, 4, "png");

    // Copy the filled ring into the image for the next iteration
    ITKHelpers::DeepCopy(internalBDSInpaintingFunctor.GetOutput(), currentImage.GetPointer());

    // Reduce the size of the target region (we "enlarge the hole", because
    // the "hole" is considered the valid part of the target mask)
    unsigned int kernelRadius = 1;
    currentTargetMask->ExpandHole(kernelRadius);

    ringCounter++;
  }

  ITKHelpers::DeepCopy(currentImage.GetPointer(), this->Output.GetPointer());
}

#endif
