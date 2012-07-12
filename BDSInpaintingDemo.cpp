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
#include "PatchMatch/Mask/Mask.h"
#include "PatchMatch/Mask/ITKHelpers/ITKHelpers.h"

#include "BDSInpainting.h"

typedef itk::Image<itk::CovariantVector<float, 3>, 2> ImageType;

int main(int argc, char*argv[])
{
  if(argc < 4)
  {
    std::cerr << "Required arguments: image mask output" << std::endl;
    return EXIT_FAILURE;
  }

  std::string imageFilename = argv[1];
  std::string maskFilename = argv[2];
  std::string outputFilename = argv[3];

  std::cout << "imageFilename: " << imageFilename << std::endl;
  std::cout << "maskFilename: " << maskFilename << std::endl;
  std::cout << "outputFilename: " << outputFilename << std::endl;

  typedef itk::ImageFileReader<ImageType> ImageReaderType;
  ImageReaderType::Pointer imageReader = ImageReaderType::New();
  imageReader->SetFileName(imageFilename);
  imageReader->Update();

  Mask::Pointer mask = Mask::New();
  mask->SetHoleValue(0);
  mask->SetValidValue(255);
  mask->Read(maskFilename);

  BDSInpainting bdsInpainting;
  bdsInpainting.SetPatchRadius(7);
  bdsInpainting.SetImage(imageReader->GetOutput());
  bdsInpainting.SetMask(mask);
  bdsInpainting.SetResolutionLevels(1);
  //bdsInpainting.SetResolutionLevels(2);

  bdsInpainting.SetIterations(4);
  
  bdsInpainting.SetDownsampleFactor(.5);
  
  bdsInpainting.SetPatchMatchIterations(3);
  bdsInpainting.Compute();

  ITKHelpers::WriteRGBImage(bdsInpainting.GetOutput(), outputFilename);

  return EXIT_SUCCESS;
}
