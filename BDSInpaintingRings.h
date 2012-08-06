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
  
  typedef typename Superclass::PatchMatchFunctorType PatchMatchFunctorType;
  
  BDSInpaintingRings();

  void Inpaint();

};

#include "BDSInpaintingRings.hpp"

#endif
