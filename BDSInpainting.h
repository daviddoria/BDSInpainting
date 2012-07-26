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

  BDSInpainting();

  /** The main driver. */
  void Compute();

  /** This function does the actual work, and is called from Compute() at multiple resolutions. */
  void Compute(TImage* const image, Mask* const mask, TImage* const output);

  TImage* GetOutput();

  void SetIterations(const unsigned int iterations);

  void SetPatchRadius(const unsigned int patchRadius);

  void SetResolutionLevels(const unsigned int resolutionLevels);

  void SetPatchMatchIterations(const unsigned int patchMatchIterations);

  void SetDownsampleFactor(const float downsampleFactor);

  void SetImage(TImage* const image);

  void SetMask(Mask* const mask);

private:

  unsigned int ResolutionLevels;
  unsigned int Iterations;
  unsigned int PatchRadius;
  unsigned int PatchMatchIterations;
  float DownsampleFactor;

  typename TImage::Pointer Output;

  typename TImage::Pointer Image;

  Mask::Pointer MaskImage;
};

#include "BDSInpainting.hpp"

#endif
