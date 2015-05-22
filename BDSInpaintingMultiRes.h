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

#ifndef BDSInpaintingMultiRes_H
#define BDSInpaintingMultiRes_H

#include "BDSInpainting.h"

/** This class uses a BDSInpainting object to inpaint an image over multiple resolutions. */
template <typename TImage>
class BDSInpaintingMultiRes : public BDSInpainting<TImage>
{
public:

  typedef BDSInpainting<TImage> Superclass;

  // Inherited typedefs
  typedef typename Superclass::PatchMatchFunctorType PatchMatchFunctorType;

  /** This function does the actual work of inpainting a single level.
    * It is called from Compute() at multiple resolutions. */
  void Compute(TImage* const image, Mask* const sourceMask, Mask* const targetMask,
               typename PatchMatch<TImage>::PMImageType* previousNNField, TImage* const output);


  /** Set the number of resolution levels to use. */
  void SetResolutionLevels(const unsigned int resolutionLevels);

  /** Set the amount to downsample the image to construct the different resolutions. */
  void SetDownsampleFactor(const float downsampleFactor);

private:

  /** The number of resolutions to use. */
  unsigned int ResolutionLevels = 3;

  /** How much to downsample the image at each level. */
  float DownsampleFactor = 0.5;

};

#include "BDSInpaintingMultiRes.hpp"

#endif
