#include "quantize.h"

// Quantization per channel for NxN block
// If <all> set, then YCbCr each have Quantize operation performed on them.
// Else, only Y has Quantize operation performed on it.
std::vector<std::shared_ptr<PixelYcbcr>> quantize(std::vector<std::shared_ptr<PixelYcbcr>> pixels, int block_size, bool all) {

    if (block_size != QUANTIZEBLOCK_SIZE) {
        fprintf(stderr, "Quantization only supports %dx%d blocks!\n", QUANTIZEBLOCK_SIZE, QUANTIZEBLOCK_SIZE);
        exit(1);
    }

    for (int i = 0; i < block_size*block_size; i++) {
        pixels[i]->y = pixels[i]->y / quant_matrix[i];
        if (all) {
            pixels[i]->cr = pixels[i]->cr / quant_matrix[i];
            pixels[i]->cb = pixels[i]->cb / quant_matrix[i];
        }
    }
    return pixels;
}

// Undo quantization per channel for NxN block
// If <all> set, then YCbCr each have undo quantize operation performed on them.
// Else, only Y has undo quantize operation performed on it.
std::vector<std::shared_ptr<PixelYcbcr>> unquantize(std::vector<std::shared_ptr<PixelYcbcr>> pixels, int block_size, bool all) {

    if (block_size != QUANTIZEBLOCK_SIZE) {
        fprintf(stderr, "Undo quantization only supports %dx%d blocks!\n", QUANTIZEBLOCK_SIZE, QUANTIZEBLOCK_SIZE);
        exit(1);
    }

    for (int i = 0; i < block_size*block_size; i++) {
        pixels[i]->y = pixels[i]->y * quant_matrix[i];
        if (all) {
            pixels[i]->cr = pixels[i]->cr * quant_matrix[i];
            pixels[i]->cb = pixels[i]->cb * quant_matrix[i];
        }
    }
    return pixels;

}
