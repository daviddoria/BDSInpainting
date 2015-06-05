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

#include "InpaintingAlgorithm.h"

// ITK
#include "itkCovariantVector.h"
#include "itkImage.h"
#include "itkImageRegion.h"
#include "itkVectorImage.h"

// Submodules
#include <Mask/Mask.h>
#include <PatchMatch/PatchMatch.h>
#include <Compositor.h>

/** This class takes a nearest neighbor field and a target mask and
  * uses the coherence term from the paper "Bidirectional Similarity" to perform inpainting.
  * By optionally using more than 1 iteration, the inpainting quality should improve. */
template <typename TImage>
class BDSInpainting : public InpaintingAlgorithm<TImage>
{
public:

  typedef InpaintingAlgorithm<TImage> Superclass;

  /** Compute the nn-field for the target pixels and then composite the patches.*/
  template <typename TPatchMatchFunctor, typename TCompositor>
  void Inpaint(TPatchMatchFunctor* const patchMatchFunctor, TCompositor* const compositor);

protected:
  typedef itk::Image<bool, 2> BoolImageType;
  void ConstructValidPatchCentersImage();
  BoolImageType::Pointer ValidPatchCentersImage = BoolImageType::New();
};

#include "BDSInpainting.hpp"

#endif
