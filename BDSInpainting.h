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
#include "PatchMatch/Mask/Mask.h"

/** This class uses PatchMatch to compute the nearest neighbor field, and then the
 *  coherence term from Bidirectional Similarity to perform inpainting. */
template <typename TImage>
class BDSInpainting
{
public:

  /** Constructor. */
  BDSInpainting();

  /** The main driver. This performs downsampling and inpainting over multiple resolutions. */
  void Compute();

  /** This function does the actual work of inpainting a single level. It is called from Compute() at multiple resolutions. */
  void Compute(TImage* const image, Mask* const mask, TImage* const output);

  /** Get the resulting inpainted image. */
  TImage* GetOutput();

  /** Set the number of iterations to run. */
  void SetIterations(const unsigned int iterations);

  /** Set the patch radius. */
  void SetPatchRadius(const unsigned int patchRadius);

  /** Set the number of resolution levels to use. */
  void SetResolutionLevels(const unsigned int resolutionLevels);

  /** Set the number of PatchMatch iterations to run. */
  void SetPatchMatchIterations(const unsigned int patchMatchIterations);

  /** Set the amount to downsample the image to construct the different resolutions. */
  void SetDownsampleFactor(const float downsampleFactor);

  /** Set the image to fill. */
  void SetImage(TImage* const image);

  /** Set the mask that indicates where to fill the image. The same mask is used as the source mask and target mask in the PatchMatch
    * algorithm, indicating "use patches outside the hole to fill the pixels inside the hole". This could be generalized to allow user specified
    * regions to take source patches from, for example. */
  void SetMask(Mask* const mask);

private:

  /** The number of resolutions to use. */
  unsigned int ResolutionLevels;

  /** The number of iterations to run. */
  unsigned int Iterations;

  /** The radius of the patches to use for inpainting. */
  unsigned int PatchRadius;

  /** The number of iterations of PatchMatch to run at each iteration. */
  unsigned int PatchMatchIterations;

  /** How much to downsample the image at each level. */
  float DownsampleFactor;

  /** The output image. */
  typename TImage::Pointer Output;

  /** The image to fill. */
  typename TImage::Pointer Image;

  /** The mask indicating where to fill the image. */
  Mask::Pointer MaskImage;
};

#include "BDSInpainting.hpp"

#endif
