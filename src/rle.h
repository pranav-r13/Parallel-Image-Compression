#include <vector>
#include <map>
#include <memory>
#include "image.h"

#define MACROBLOCK_SIZE 8
#define COLOR_Y  0
#define COLOR_CR 1
#define COLOR_CB 2

// Store (encoded, count) structs. Since macroblocks are 8x8 the char datatype
// (-128 to 127) is sufficient for storing this information.
struct RleTuple {
    char encoded;
    char count;
};

struct EncodedBlockColor {
    double dc_val;
    std::shared_ptr<std::vector<RleTuple>> encoded;
    std::shared_ptr<std::map<char, double>> decode_table;
    std::shared_ptr<std::map<double, char>> encode_table;
};

struct EncodedBlock {
    std::shared_ptr<EncodedBlockColor> y;
    std::shared_ptr<EncodedBlockColor> cr;
    std::shared_ptr<EncodedBlockColor> cb;
};

struct JpegEncoded { // for helper function ease of use
    unsigned int width;
    unsigned int height;
    std::vector<std::shared_ptr<EncodedBlock>> encodedBlocks;
};

// MIRROR STRUCTURES FOR STACK ALLOC MEMORY MATH IN MPI
struct EncodedBlockColorNoPtr {
    double dc_val;
    char encoded_len;
    RleTuple encoded[MACROBLOCK_SIZE * MACROBLOCK_SIZE];
    char table_size;
    char char_vals[MACROBLOCK_SIZE * MACROBLOCK_SIZE];
    double double_vals[MACROBLOCK_SIZE * MACROBLOCK_SIZE];
};

// MIRROR STRUCTURES FOR STACK ALLOC MEMORY MATH IN MPI
struct EncodedBlockNoPtr {
    EncodedBlockColorNoPtr y;
    EncodedBlockColorNoPtr cr;
    EncodedBlockColorNoPtr cb;
};

std::shared_ptr<EncodedBlock> RLE(
    std::vector<std::shared_ptr<PixelYcbcr>> block,
    int block_size
);

std::vector<double> extractChannel(
    std::vector<std::shared_ptr<PixelYcbcr>> block,
    int chan
);

std::shared_ptr<EncodedBlockColor> buildTable(
    std::vector<std::shared_ptr<PixelYcbcr>> block,
    int chan,
    int block_size
);

void encodeValues(
    std::vector<std::shared_ptr<PixelYcbcr>> block,
    std::shared_ptr<EncodedBlockColor> color,
    int chan
);

std::vector<std::shared_ptr<PixelYcbcr>> decodeRLE(
    std::shared_ptr<EncodedBlock> encoded,
    int block_size
);

void writeToBuffer(
    std::shared_ptr<EncodedBlockNoPtr> encodedBlockBuffer,
    std::vector<std::shared_ptr<EncodedBlock>> encodedBlocks,
    int idx, int chan
);

std::vector<std::shared_ptr<EncodedBlock>> convertBufferToEncodedBlocks(
    std::shared_ptr<EncodedBlockNoPtr> encodedBlocksBuffer, int numEncodedBlocks
);
