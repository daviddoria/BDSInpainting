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

// Custom
#include "BDSInpaintingRings.h"

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
  typedef PoissonEditing<typename TypeTraits<ImageType::PixelType>::ComponentType> PoissonEditingType;

  typename PoissonEditingType::GuidanceFieldType::Pointer zeroGuidanceField =
            PoissonEditingType::GuidanceFieldType::New();
  zeroGuidanceField->SetRegions(image->GetLargestPossibleRegion());
  zeroGuidanceField->Allocate();
  typename PoissonEditingType::GuidanceFieldType::PixelType zeroPixel;
  zeroPixel.Fill(0);
  ITKHelpers::SetImageToConstant(zeroGuidanceField.GetPointer(), zeroPixel);

  PoissonEditingType::FillImage(image, targetMask,
                                zeroGuidanceField.GetPointer(), image);

  ITKHelpers::WriteRGBImage(image, "PoissonFilled.png");

  // PatchMatch requires that the target region be specified by valid pixels
  targetMask->InvertData();

  // Setup the patch distance functor
  SSD<ImageType> ssdFunctor;
  ssdFunctor.SetImage(image);

  // Setup the PatchMatch functor
  PatchMatch<ImageType> patchMatchFunctor;
  patchMatchFunctor.SetPatchRadius(patchRadius);
  patchMatchFunctor.SetPatchDistanceFunctor(&ssdFunctor);
  patchMatchFunctor.SetIterations(20);
  patchMatchFunctor.SetInitializationStrategy(PatchMatch<ImageType>::RANDOM);

  // Test the result of PatchMatch here
//   patchMatchFunctor.SetImage(image);
//   patchMatchFunctor.SetTargetMask(targetMask);
//   patchMatchFunctor.SetSourceMask(sourceMask);
   // patchMatchFunctor.SetRandom(false);
//   patchMatchFunctor.Compute(NULL);

  // Here, the source match and target match are the same, specifying the classicial
  // "use pixels outside the hole to fill the pixels inside the hole".
  // In an interactive algorith, the user could manually specify a source region,
  // improving the resulting inpainting.
  BDSInpaintingRings<ImageType> bdsInpainting;
  bdsInpainting.SetPatchRadius(patchRadius);
  bdsInpainting.SetImage(image);
  bdsInpainting.SetSourceMask(sourceMask);
  bdsInpainting.SetTargetMask(targetMask);

  // This means simply fill the image at the original resolution without any downsampling
  bdsInpainting.SetResolutionLevels(1);

  //bdsInpainting.SetResolutionLevels(2);

  bdsInpainting.SetIterations(1);
  //bdsInpainting.SetIterations(4);

  //bdsInpainting.SetCompositingMethod(BDSInpainting<ImageType>::WEIGHTED_AVERAGE);
  bdsInpainting.SetCompositingMethod(BDSInpainting<ImageType>::AVERAGE);
  //bdsInpainting.SetCompositingMethod(BDSInpainting<ImageType>::BEST_PATCH);
  //bdsInpainting.SetCompositingMethod(BDSInpainting<ImageType>::CLOSEST_TO_AVERAGE);

  bdsInpainting.SetDownsampleFactor(.5);

  bdsInpainting.SetPatchMatchFunctor(&patchMatchFunctor);
  bdsInpainting.Inpaint();

  ITKHelpers::WriteRGBImage(bdsInpainting.GetOutput(), outputFilename);

  return EXIT_SUCCESS;
}
