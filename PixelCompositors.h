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

#ifndef PixelCompositors_H
#define PixelCompositors_H

struct PixelCompositorAverage
{
  /** Composite by averaging pixels. */
  template <typename TPixel>
  static TPixel Composite(
    const std::vector<TPixel>& contributingPixels,
    const std::vector<float>& contributingScores)
  {
    TPixel newValue = Statistics::Average(contributingPixels);
    return newValue;
  }
};

struct PixelCompositorWeightedAverage
{
  /** Composite by averaging pixels. */
  template <typename TPixel>
  static TPixel Composite(
    const std::vector<TPixel>& contributingPixels,
    const std::vector<float>& contributingScores)
  {
    // If there is only one element, simply return it.
    if(contributingScores.size() == 1)
    {
      return contributingPixels[0];
    }

    // The weights should be inversely proportional to the patch errors/scores.
    // That is, a patch with a high error should get a low weight. We accomplish this by
    // making the weight equal to 1 - (value - min) / |range|

    float minValue = *std::min_element(contributingScores.begin(), contributingScores.end());
    float maxValue = *std::max_element(contributingScores.begin(), contributingScores.end());
    float range = maxValue - minValue;
    // If the range is zero, all elements are the same so just return the first one
    if(range == 0.0f)
    {
      return contributingPixels[0];
    }

    std::vector<float> weights(contributingScores.size());
    for(unsigned int i = 0; i < contributingScores.size(); ++i)
    {
      weights[i] = 1.0f - (contributingScores[i] - minValue) / range;
    }

    // Make the weights sum to 1
    Helpers::NormalizeVectorInPlace(weights);

    TPixel newValue = Helpers::WeightedAverage(contributingPixels, weights);
    return newValue;
  }
};

struct PixelCompositorClosestToAverage
{
  /** Composite by averaging pixels. */
  template <typename TPixel>
  static TPixel Composite(
    const std::vector<TPixel>& contributingPixels,
    const std::vector<float>& contributingScores)
  {
    // Use the pixel closest to the average pixel
    TPixel averagePixel = Statistics::Average(contributingPixels);
    unsigned int patchId = ITKHelpers::ClosestValueIndex(contributingPixels, averagePixel);
    TPixel newValue = contributingPixels[patchId];
    return newValue;
  }
};

struct PixelCompositorBestPatch
{
  /** Composite by averaging pixels. */
  template <typename TPixel>
  static TPixel Composite(
    const std::vector<TPixel>& contributingPixels,
    const std::vector<float>& contributingScores)
  {
    // Take the pixel from the best matching patch
    unsigned int patchId = Helpers::Argmin(contributingScores);
    TPixel newValue = contributingPixels[patchId];
    return newValue;
  }
};

#endif
