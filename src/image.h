#include <vector>
#include <memory>

#ifndef IMAGE_H
#define IMAGE_H

struct PixelRgba {
    unsigned char r;
    unsigned char g;
    unsigned char b;
    unsigned char a;
};

struct PixelYcbcr {
    double y;
    double cb;
    double cr;
};

struct ImageRgb {
    std::vector<std::shared_ptr<PixelRgba>> pixels;
    int numPixels;
    int width;
    int height;
};

struct ImageYcbcr {
    std::vector<std::shared_ptr<PixelYcbcr>> pixels;
    int numPixels;
    int width;
    int height;
};

struct ImageBlocks {
    std::vector<std::vector<std::shared_ptr<PixelYcbcr>>> blocks;
    int numBlocks;
    int width;
    int height;
};

std::shared_ptr<ImageRgb> convertBytesToImage(std::vector<unsigned char> bytes, unsigned int width, unsigned int height, int start = 0, int end = -1);
std::vector<unsigned char> convertImageToBytes(std::shared_ptr<ImageRgb> image);

std::shared_ptr<ImageYcbcr> convertRgbToYcbcr(std::shared_ptr<ImageRgb> input);
std::shared_ptr<ImageRgb> convertYcbcrToRgb(std::shared_ptr<ImageYcbcr> input);

std::shared_ptr<ImageBlocks> convertYcbcrToBlocks(std::shared_ptr<ImageYcbcr> input, int block_size);
std::shared_ptr<ImageYcbcr> convertBlocksToYcbcr(std::shared_ptr<ImageBlocks> input, int block_size);

void downsampleCbcr(std::shared_ptr<ImageBlocks> image, int block_size);
void upsampleCbcr(std::shared_ptr<ImageBlocks> image, int block_size);

// image utils

struct Coord {
    int col;
    int row;
};

bool pixel_in_bounds(int row, int col, int width, int height);

// Convert from (x,y) given size NxN array to vectorized idx
int sub2ind(int width, int col, int row);
int sub2ind(int width, Coord coord);

// Convert from vectorized idx to (x,y) given size NxN array
Coord ind2sub(int width, int idx);

#endif
