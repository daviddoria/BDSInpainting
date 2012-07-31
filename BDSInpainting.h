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

#ifndef BDSInpainting_H
#define BDSInpainting_H

// ITK
#include "itkCovariantVector.h"
#include "itkImage.h"
#include "itkImageRegion.h"
#include "itkVectorImage.h"

// Submodules
#include "Mask/Mask.h"
#include "PatchMatch/PatchMatch.h"

/** This class uses PatchMatch to compute the nearest neighbor field, and then the
 *  coherence term from Bidirectional Similarity to perform inpainting. */
template <typename TImage>
class BDSInpainting
{
public:

  /** Constructor. */
  BDSInpainting();

  /** Choices of compositing method. */
  enum CompositingMethodEnum {BEST_PATCH, WEIGHTED_AVERAGE, AVERAGE, CLOSEST_TO_AVERAGE};

  /** Set the compositing method to use. */
  void SetCompositingMethod(const CompositingMethodEnum& compositingMethod);

  /** The main driver. This performs downsampling and inpainting over multiple resolutions. */
  void Inpaint();

  /** This function does the actual work of inpainting a single level.
    * It is called from Compute() at multiple resolutions. */
  virtual void Compute(TImage* const image, Mask* const sourceMask, Mask* const targetMask,
                       typename PatchMatch<TImage>::PMImageType* previousNNField, TImage* const output);

  /** Get the resulting inpainted image. */
  TImage* GetOutput();

  /** Set the number of iterations to run. */
  void SetIterations(const unsigned int iterations);

  /** Set the patch radius. */
  void SetPatchRadius(const unsigned int patchRadius);

  /** Set the number of resolution levels to use. */
  void SetResolutionLevels(const unsigned int resolutionLevels);

  /** Set the PatchMatch functor to use. */
  void SetPatchMatchFunctor(PatchMatch<TImage>* patchMatchFunctor);

  /** Set the amount to downsample the image to construct the different resolutions. */
  void SetDownsampleFactor(const float downsampleFactor);

  /** Set the image to fill. */
  void SetImage(TImage* const image);

  /** Set the mask that indicates where to take source patches from. Source patches
    * are patches that are entirely in the Valid region.*/
  void SetSourceMask(Mask* const mask);

  /** Set the mask that indicates where to fill the image. Pixels in the Hole region should be filled.*/
  void SetTargetMask(Mask* const mask);

protected:

  /** The number of resolutions to use. */
  unsigned int ResolutionLevels;

  /** The number of iterations to run. */
  unsigned int Iterations;

  /** The radius of the patches to use for inpainting. */
  unsigned int PatchRadius;

  /** The PatchMatch functor to use. */
  PatchMatch<TImage>* PatchMatchFunctor;

  /** How much to downsample the image at each level. */
  float DownsampleFactor;

  /** The output image. */
  typename TImage::Pointer Output;

  /** The image to fill. */
  typename TImage::Pointer Image;

  /** The mask where Hole pixels indicate the pixels to fill. */
  Mask::Pointer TargetMask;

  /** The mask where only patches in the Valid region are considered as potential matches. */
  Mask::Pointer SourceMask;

  /** Using the 'nnField', compute new values for the hole pixels in the 'targetMask' and
    * store them in 'updatedImage' */
  void UpdatePixels(const TImage* const oldImage,
                    const Mask* const targetMask,
                    const typename PatchMatch<TImage>::PMImageType* const nnField,
                    TImage* const updatedImage);

  /** The compositing method to use. */
  CompositingMethodEnum CompositingMethod;

  /** Composite using the specified CompositingMethod. */
  typename TImage::PixelType Composite(
     const std::vector<typename TImage::PixelType>& contributingPixels,
     const std::vector<float>& contributingScores);

  /** Composite by averaging pixels. */
  typename TImage::PixelType CompositeAverage(
    const std::vector<typename TImage::PixelType>& contributingPixels);

  /** Composite by weighted averaging. */
  typename TImage::PixelType CompositeWeightedAverage(
     const std::vector<typename TImage::PixelType>& contributingPixels,
     const std::vector<float>& contributingScores);

  /** Composite by choosing the pixel closest to the average. */
  typename TImage::PixelType CompositeClosestToAverage(
     const std::vector<typename TImage::PixelType>& contributingPixels,
     const std::vector<float>& contributingScores);

  /** Composite by choosing the pixel from the best patch. */
  typename TImage::PixelType CompositeBestPatch(
     const std::vector<typename TImage::PixelType>& contributingPixels,
     const std::vector<float>& contributingScores);
};

#include "BDSInpainting.hpp"

#endif
