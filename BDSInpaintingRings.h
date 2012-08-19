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

/** This class uses composition (uses BDSInpainting objects internally)
 *  to compute the nearest neighbor field one ring at a time, from the outside
 *  in, compositing as it goes along.. */
template <typename TImage>
class BDSInpaintingRings : public InpaintingAlgorithm<TImage>
{
public:
  typedef InpaintingAlgorithm<TImage> Superclass;

  BDSInpaintingRings();

  /** Perform the NNField computation and compositing for the entire hole
    * (and the boundary around it, as prescribed by ExpandMask() ) */
  void Inpaint();

private:

  /** Get the "patch-radius-thick ring" around the original hole. We do not
    * need to composite in this region,
    * but we do need to compute the NNField here (as it is non-trivial
    * (since the trivial NN (itself) is not completely a valid patch),
    * and they will be used in the compositing) */
  void GetSurroundingRingMask(Mask* surroundingMask);

  /** Initialize the known region by setting the NN's of each pixel whose surrounding
    * patch is entirely outside the hole. */
  void InitializeKnownRegion();

  /** Shrink the source region around the image boundary and around the hole.
    * The idea is to allow a patch in this region to be propagated the width of the ring
    * without running into the hole. */
  void ProducePropagationBuffer();

  /** Run several iterations of the PatchMatch algorithm with neighbor-histogram difference
    * verification */
  void ConstrainedPatchMatch(Mask* const targetMask, const float histogramRatioStart,
                             const float histogramRatioStep, const float maxHistogramRatio);

  /** Run forced propagation until either the target region is filled or the
    * propagation has been completely restricted by hole geometry and patch
    * selection location. */
  void ForcePropagation(Mask* const targetMask);

  /** Perform a combination of propagation and random search steps, and composite the result. */
  void FillHole(Mask* const targetMask);

  /** Remove all matches from the MatchSet at pixels which do not have a verified match. */
  void ClearUnverifiedPixels();

  /** Fill the hole one one-pixel-thick ring at a time (from outside in) */
  void SinglePixelRings();

  /** Fill the hole one patch-radius-thick ring at a time (from outside in) */
  void PatchRadiusThickRings();

  /** Compute the NNField using a combination of verified propagation, random search, and forced propagation steps. */
  void ComputeNNField(Mask* const targetMask);

};

#include "BDSInpaintingRings.hpp"

#endif
