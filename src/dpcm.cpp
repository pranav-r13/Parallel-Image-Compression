#include "image.h"
#include "dpcm.h"

// Take in set of macroblocks post DCT and quantization to compress.
// Updates macroblocks to encode DC values i.e. (0,0) in each block
// as a series of deltas, where first block value is an actual value.
//
// Updates values in-place.
void DPCM(std::vector<std::vector<std::shared_ptr<PixelYcbcr>>> blocks) {
    for (unsigned int i = blocks.size() - 1; i > 0; i--) {
        std::shared_ptr<PixelYcbcr> prev = blocks[i-1][0];
        std::shared_ptr<PixelYcbcr> curr = blocks[i][0];
        blocks[i][0]->y = curr->y - prev->y;
        blocks[i][0]->cr = curr->cr - prev->cr;
        blocks[i][0]->cb = curr->cb - prev->cb;
    }
}

void unDPCM(std::vector<std::vector<std::shared_ptr<PixelYcbcr>>> blocks) {
    for (unsigned int i = 1; i < blocks.size(); i++) {
        std::shared_ptr<PixelYcbcr> prev = blocks[i-1][0];
        std::shared_ptr<PixelYcbcr> curr = blocks[i][0];
        blocks[i][0]->y = curr->y + prev->y;
        blocks[i][0]->cr = curr->cr + prev->cr;
        blocks[i][0]->cb = curr->cb + prev->cb;
    }
}
