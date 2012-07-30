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

#ifndef BDSInpaintingRings_H
#define BDSInpaintingRings_H

#include "BDSInpainting.h"

/** This class uses PatchMatch to compute the nearest neighbor field, and then the
 *  coherence term from Bidirectional Similarity to perform inpainting. */
template <typename TImage>
class BDSInpaintingRings : public BDSInpainting<TImage>
{
public:

  /** Constructor. */
  BDSInpaintingRings();

  /** This function does the actual work of inpainting a single level.
    * It is called from Compute() at multiple resolutions. */
  void Compute(TImage* const image, Mask* const sourceMask, Mask* const targetMask, TImage* const output);

private:

};

#include "BDSInpaintingRings.hpp"

#endif
