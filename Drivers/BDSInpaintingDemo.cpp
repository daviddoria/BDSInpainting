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

// STL
#include <iostream>

// ITK
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkCovariantVector.h"

// Submodules
#include <Mask/Mask.h>

#include <ITKHelpers/ITKHelpers.h>

#include <PatchComparison/SSD.h>

#include <PatchMatch/PatchMatch.h>
#include <PatchMatch/Propagator.h>
#include <PatchMatch/RandomSearch.h>

#include <PoissonEditing/PoissonEditingWrappers.h>

// Custom
#include "BDSInpainting.h"
#include "Compositor.h"
#include "PixelCompositors.h"

int main(int argc, char*argv[])
{
  // Parse the input
  if(argc < 6)
  {
    std::cerr << "Required arguments: image sourceMask.mask targetMask.mask patchRadius output" << std::endl;
    return EXIT_FAILURE;
  }

  std::stringstream ss;
  for(int i = 1; i < argc; ++i)
  {
    ss << argv[i] << " ";
  }

  std::string imageFilename;
  std::string sourceMaskFilename;
  std::string targetMaskFilename;
  unsigned int patchRadius;
  std::string outputFilename;

  ss >> imageFilename >> sourceMaskFilename >> targetMaskFilename >> patchRadius >> outputFilename;

  // Output the parsed values
  std::cout << "imageFilename: " << imageFilename << std::endl
            << "sourceMaskFilename: " << sourceMaskFilename << std::endl
            << "targetMaskFilename: " << targetMaskFilename << std::endl
            << "patchRadius: " << patchRadius << std::endl
            << "outputFilename: " << outputFilename << std::endl;

  typedef itk::Image<itk::CovariantVector<unsigned char, 3>, 2> ImageType;

  // Read the image and the masks
  typedef itk::ImageFileReader<ImageType> ImageReaderType;
  ImageReaderType::Pointer imageReader = ImageReaderType::New();
  imageReader->SetFileName(imageFilename);
  imageReader->Update();

  ImageType* image = imageReader->GetOutput();

  Mask::Pointer sourceMask = Mask::New();
  sourceMask->Read(sourceMaskFilename);

  Mask::Pointer targetMask = Mask::New();
  targetMask->Read(targetMaskFilename);

  //std::cout << "target mask has " << targetMask->CountHolePixels() << " hole pixels." << std::endl;

  // Poisson fill the input image


  typename PoissonEditingParent::GuidanceFieldType::Pointer zeroGuidanceField =
            PoissonEditingParent::GuidanceFieldType::New();
  zeroGuidanceField->SetRegions(image->GetLargestPossibleRegion());
  zeroGuidanceField->Allocate();
  typename PoissonEditingParent::GuidanceFieldType::PixelType zeroPixel;
  zeroPixel.Fill(0);
  ITKHelpers::SetImageToConstant(zeroGuidanceField.GetPointer(), zeroPixel);

  ImageType::Pointer filledImage = ImageType::New();

  FillImage(image, targetMask, zeroGuidanceField, filledImage.GetPointer(),
            image->GetLargestPossibleRegion());

  ITKHelpers::WriteRGBImage(filledImage.GetPointer(), "PoissonFilled.png");

  // PatchMatch requires that the target region be specified by valid pixels
  targetMask->InvertData();

  // Setup the patch distance functor
  typedef SSD<ImageType> PatchDistanceFunctorType;
  PatchDistanceFunctorType* patchDistanceFunctor = new PatchDistanceFunctorType;
  patchDistanceFunctor->SetImage(image);

  typedef Propagator<PatchDistanceFunctorType> PropagatorType;
  PropagatorType* propagator = new PropagatorType;
  propagator->SetPatchDistanceFunctor(patchDistanceFunctor);

  typedef RandomSearch<ImageType, PatchDistanceFunctorType> RandomSearchType;
  RandomSearchType* randomSearchFunctor = new RandomSearchType;
  randomSearchFunctor->SetPatchDistanceFunctor(patchDistanceFunctor);
  randomSearchFunctor->SetImage(image);
  randomSearchFunctor->SetPatchRadius(patchRadius);


  // Setup the PatchMatch functor
  PatchMatch<ImageType, PropagatorType, RandomSearchType> patchMatchFunctor;
  patchMatchFunctor.SetPatchRadius(patchRadius);
  patchMatchFunctor.SetIterations(1);
  patchMatchFunctor.SetPropagationFunctor(propagator);
  patchMatchFunctor.SetRandomSearchFunctor(randomSearchFunctor);

  // Test the result of PatchMatch here
  //patchMatchFunctor.SetRandom(false);

  // Here, the source match and target match are the same, specifying the classicial
  // "use pixels outside the hole to fill the pixels inside the hole".
  // In an interactive algorith, the user could manually specify a source region,
  // improving the resulting inpainting.
  BDSInpainting<ImageType> bdsInpainting;
  bdsInpainting.SetPatchRadius(patchRadius);
  bdsInpainting.SetImage(image);
  bdsInpainting.SetSourceMask(sourceMask);
  bdsInpainting.SetTargetMask(targetMask);

  bdsInpainting.SetIterations(1);
  //bdsInpainting.SetIterations(4);

  Compositor<ImageType, PixelCompositorAverage> compositor;
  bdsInpainting.Inpaint(&patchMatchFunctor, &compositor);

  ITKHelpers::WriteRGBImage(bdsInpainting.GetOutput(), outputFilename);

  return EXIT_SUCCESS;
}
