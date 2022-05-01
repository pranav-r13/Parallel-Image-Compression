#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES // Pi
#endif

#include "math.h"
#include "image.h"

// Forward DCT operation for NxN block
// If <all> set, then YCbCr each have DCT performed on them.
// Else, only Y has DCT performed on it.
std::vector<std::shared_ptr<PixelYcbcr>> DCT(std::vector<std::shared_ptr<PixelYcbcr>> pixels, int block_size, bool all);

// Inverse DCT operation for NxN block
// If <all> set, then YCbCr each have IDCT performed on them.
// Else, only Y has IDCT performed on it.
std::vector<std::shared_ptr<PixelYcbcr>> IDCT(std::vector<std::shared_ptr<PixelYcbcr>> pixels, int block_size, bool all);

