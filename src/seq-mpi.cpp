#include <iostream>
#include <fstream>
#include <cstdarg>
#include <string>
#include "CycleTimer.h"
#include "getopt.h"
#include "stdio.h"
#include "lodepng/lodepng.h"
#include "dct.h"
#include "quantize.h"
#include "dpcm.h"
#include "rle.h"
#include "mpi.h"

#define MACROBLOCK_SIZE 8

#ifndef LOGLEVEL
#define LOGLEVEL 0 // set to 1 for logging
#endif

void log(int rank, const char* format, ...) {
    if (rank != 0) {
        return;
    }
    va_list args;
    va_start(args, format);
    if (LOGLEVEL) {
        fprintf(stdout, format, args);
    }
    va_end(args);
}

std::shared_ptr<JpegEncoded> jpegSeq(const char* infile, const char* compressedFile) {
    fprintf(stdout, "running sequential version\n");

    double startTime = CycleTimer::currentSeconds();

    std::vector<unsigned char> bytes; // The raw pixels
    unsigned int width, height;

    double loadImageStartTime = CycleTimer::currentSeconds();
    unsigned int error = lodepng::decode(bytes, width, height, infile);
    double loadImageStopTime = CycleTimer::currentSeconds();

    if(error) {
      std::cout << "decoder error " << error << ": " << lodepng_error_text(error) << std::endl;
    } else {
      log(0, "success decoding %s!\n", infile);
    }

    // 4 bytes per pixel, ordered RGBARGBA
    log(0, "convertBytesToImage()...\n");
    double convertBytesToImageStartTime = CycleTimer::currentSeconds();
    std::shared_ptr<ImageRgb> imageRgb = convertBytesToImage(bytes, width, height);
    double convertBytesToImageEndTime = CycleTimer::currentSeconds();

    log(0, "convertRgbToYcbcr()...\n");
    double convertRgbToYcbcrStartTime = CycleTimer::currentSeconds();
    std::shared_ptr<ImageYcbcr> imageYcbcr = convertRgbToYcbcr(imageRgb);
    double convertRgbToYcbcrEndTime = CycleTimer::currentSeconds();

    log(0, "convertYcbcrToBlocks()...\n");
    double convertYcbcrToBlocksStartTime = CycleTimer::currentSeconds();
    std::shared_ptr<ImageBlocks> imageBlocks = convertYcbcrToBlocks(imageYcbcr, MACROBLOCK_SIZE);
    double convertYcbcrToBlocksEndTime = CycleTimer::currentSeconds();

    log(0, "DCT()...\n");
    double dctStartTime = CycleTimer::currentSeconds();
    std::vector<std::vector<std::shared_ptr<PixelYcbcr>>> dcts;
    for (auto block : imageBlocks->blocks) {
        dcts.push_back(DCT(block, MACROBLOCK_SIZE, true));
    }
    double dctEndTime = CycleTimer::currentSeconds();

    log(0, "quantize()...\n");
    double quantizeStartTime = CycleTimer::currentSeconds();
    std::vector<std::vector<std::shared_ptr<PixelYcbcr>>> quantizedBlocks;
    for (auto dct : dcts) {
        quantizedBlocks.push_back(quantize(dct, MACROBLOCK_SIZE, true));
    }
    double quantizeEndTime = CycleTimer::currentSeconds();

    log(0, "DPCM()...\n");
    double dpcmStartTime = CycleTimer::currentSeconds();
    DPCM(quantizedBlocks);
    double dpcmEndTime = CycleTimer::currentSeconds();

    log(0, "RLE()...\n");
    double rleStartTime = CycleTimer::currentSeconds();
    std::vector<std::shared_ptr<EncodedBlock>> encodedBlocks;
    for (auto quantizedBlock : quantizedBlocks) {
        encodedBlocks.push_back(RLE(quantizedBlock, MACROBLOCK_SIZE));
    }
    double rleEndTime = CycleTimer::currentSeconds();

    log(0, "done encoding!\n");
    log(0, "writing to file...\n");
    double writeCompressedImageStartTime = CycleTimer::currentSeconds();
    std::ofstream jpegFile(compressedFile);
    for (const auto &block : encodedBlocks) {
        jpegFile << block;
    }
    double writeCompressedImageEndTime = CycleTimer::currentSeconds();
    log(0, "jpeg stored!\n");

    std::shared_ptr<JpegEncoded> result = std::make_shared<JpegEncoded>();
    result->encodedBlocks = encodedBlocks;
    result->width = width;
    result->height = height;

    double endTime = CycleTimer::currentSeconds();

    fprintf(stdout,
    "=======================================\n"
    "= Sequential encoding performance: \n"
    "=======================================\n"
    "Load Image: %.3fs\n"
    "Convert Bytes to Image: %.3fs\n"
    "Convert RGB to YCbCr: %.3fs\n"
    "Convert YCbCr to Blocks: %.3fs\n"
    "DCT: %.3fs\n"
    "Quantize: %.3fs\n"
    "DPCM: %.3fs\n"
    "RLE: %.3fs\n"
    "Encode Compressed Image: %.3fs\n"
    "Total time: %.3fs\n",
    loadImageStopTime - loadImageStartTime,
    convertBytesToImageEndTime - convertBytesToImageStartTime,
    convertRgbToYcbcrEndTime - convertRgbToYcbcrStartTime,
    convertYcbcrToBlocksEndTime - convertYcbcrToBlocksStartTime,
    dctEndTime - dctStartTime,
    quantizeEndTime - quantizeStartTime,
    dpcmEndTime - dpcmStartTime,
    rleEndTime - rleStartTime,
    writeCompressedImageEndTime - writeCompressedImageStartTime,
    endTime - startTime);

    return result;
}

std::vector<unsigned char> jpegDecodeSeq(std::shared_ptr<JpegEncoded> jpegEncoded, const char* outfile) {

    unsigned int width = jpegEncoded->width;
    unsigned int height = jpegEncoded->height;
    std::vector<std::shared_ptr<EncodedBlock>> encodedBlocks = jpegEncoded->encodedBlocks;

    log(0, "==============\n");
    log(0, "now let's undo the process...\n");

    log(0, "undoing RLE()...\n");
    std::vector<std::vector<std::shared_ptr<PixelYcbcr>>> decodedQuantizedBlocks;
    for (auto encodedBlock : encodedBlocks) {
        decodedQuantizedBlocks.push_back(decodeRLE(encodedBlock, MACROBLOCK_SIZE));
    }

    log(0, "undoing DPCM()...\n");
    unDPCM(decodedQuantizedBlocks);

    log(0, "undoing quantize()...\n");
    std::vector<std::vector<std::shared_ptr<PixelYcbcr>>> unquantizedBlocks;
    for (auto decodedQuantizedBlock : decodedQuantizedBlocks) {
        unquantizedBlocks.push_back(unquantize(decodedQuantizedBlock, MACROBLOCK_SIZE, true));
    }

    log(0, "undoing DCT()...\n");
    std::vector<std::vector<std::shared_ptr<PixelYcbcr>>> idcts;
    for (auto unquantized : unquantizedBlocks) {
        idcts.push_back(IDCT(unquantized, MACROBLOCK_SIZE, true));
    }

    log(0, "undoing convertYcbcrToBlocks()...\n");
    std::shared_ptr<ImageBlocks> imageBlocksIdct(new ImageBlocks);
    imageBlocksIdct->blocks = idcts;
    imageBlocksIdct->width = width;
    imageBlocksIdct->height = height;
    std::shared_ptr<ImageYcbcr> imgFromBlocks = convertBlocksToYcbcr(imageBlocksIdct, MACROBLOCK_SIZE);

    log(0, "undoing convertRgbToYcbcr()...\n");
    std::shared_ptr<ImageRgb> imageRgbRecovered = convertYcbcrToRgb(imgFromBlocks);

    log(0, "undoing convertBytesToImage()...\n");
    std::vector<unsigned char> imgRecovered = convertImageToBytes(imageRgbRecovered);

    unsigned int error = lodepng::encode(outfile, imgRecovered, width, height);

    if(error) {
        std::cout << "encoder error " << error << ": "<< lodepng_error_text(error) << std::endl;
    } else {
        fprintf(stdout, "success encoding to %s!\n", outfile);
    }

    return imgRecovered;
}

void encodeSeq(const char* infile, const char* outfile, const char* compressedFile) {

    std::shared_ptr<JpegEncoded> jpegEncoded = jpegSeq(infile, compressedFile);
    std::vector<std::shared_ptr<EncodedBlock>> encodedBlocks = jpegEncoded->encodedBlocks;
    std::vector<unsigned char> imgRecovered = jpegDecodeSeq(jpegEncoded, outfile);

}

void encodeMpi(const char* infile, const char* outfile, const char* compressedFile) {

    // Start parallel area
    MPI_Status mpiStatus;
    int numTasks, rank;
    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &numTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    double endTime = CycleTimer::currentSeconds();

    log(rank, "running MPI version\n");
    double startTime = CycleTimer::currentSeconds();
    std::vector<unsigned char> bytes; // The raw pixels
    unsigned int width, height;

    double loadImageStartTime = CycleTimer::currentSeconds();
    unsigned int error = lodepng::decode(bytes, width, height, infile);
    double loadImageEndTime = CycleTimer::currentSeconds();

    if(error) {
        std::cout << "decoder error " << error << ": " << lodepng_error_text(error) << std::endl;
    } else {
        if (rank == 0) {
            log(0, "success decoding %s!\n", infile);
        }
    }

    // Begin setup MPI structs
    double mpiSetupStartTime = CycleTimer::currentSeconds();
    // set up mpi pixelycbcr
    // float y, cb, cr
    MPI_Datatype MPI_PixelYcbcr;
    MPI_Type_contiguous(3, MPI_DOUBLE, &MPI_PixelYcbcr);
    MPI_Type_commit(&MPI_PixelYcbcr);

    // Set up RleTuple datatype
    MPI_Datatype MPI_RleTuple;
    MPI_Type_contiguous(2, MPI_CHAR, &MPI_RleTuple);
    MPI_Type_commit(&MPI_RleTuple);

    // Set up EncodedBlockColor
    // Change to data structure: convert std::vector to {char size, [RleTuple]}
    // Change to data stucture: convert std::map to
    //      {char num_entries, [char], [double]}
    // New struct:
    //      - double dc_val
    //      - char encoded_len
    //      - rleTupleVector encoded
    //      - char table_size
    //      - charVector char_vals
    //      - doubleVector double_vals

    // 1. Set up encoded[RleTuple]
    MPI_Datatype MPI_RleTupleVector;
    MPI_Type_contiguous(MACROBLOCK_SIZE * MACROBLOCK_SIZE, MPI_RleTuple, &MPI_RleTupleVector);
    MPI_Type_commit(&MPI_RleTupleVector);

    // 2. Set up vector for map: char values
    MPI_Datatype MPI_CharVector;
    MPI_Type_contiguous(MACROBLOCK_SIZE * MACROBLOCK_SIZE, MPI_CHAR, &MPI_CharVector);
    MPI_Type_commit(&MPI_CharVector);

    // 3. Set up vector for map: double values
    MPI_Datatype MPI_DoubleVector;
    MPI_Type_contiguous(MACROBLOCK_SIZE * MACROBLOCK_SIZE, MPI_DOUBLE, &MPI_DoubleVector);
    MPI_Type_commit(&MPI_DoubleVector);

    // 4. Set up structure for EncodedBlockColor
    int encodedBlockColorLen = 6;
    MPI_Datatype MPI_EncodedBlockColor, encodedBlockColorTypes[encodedBlockColorLen];
    int encodedBlockColorBlocks[encodedBlockColorLen];
    MPI_Aint encodedBlockColorOffsets[encodedBlockColorLen], lb, extent;
    // double dc_val
    encodedBlockColorOffsets[0] = 0;
    encodedBlockColorTypes[0] = MPI_DOUBLE;
    encodedBlockColorBlocks[0] = 1;
    // char encoded_len
    MPI_Type_get_extent(MPI_DOUBLE, &lb, &extent);
    encodedBlockColorOffsets[1] = encodedBlockColorOffsets[0] + extent;
    encodedBlockColorTypes[1] = MPI_CHAR;
    encodedBlockColorBlocks[1] = 1;
    // rleTupleVector encoded
    MPI_Type_get_extent(MPI_CHAR, &lb, &extent);
    encodedBlockColorOffsets[2] = encodedBlockColorOffsets[1] + extent;
    encodedBlockColorTypes[2] = MPI_RleTupleVector;
    encodedBlockColorBlocks[2] = 1;
    // char table_size
    MPI_Type_get_extent(MPI_RleTupleVector, &lb, &extent);
    encodedBlockColorOffsets[3] = encodedBlockColorOffsets[2] + extent;
    encodedBlockColorTypes[3] = MPI_CHAR;
    encodedBlockColorBlocks[3] = 1;
    // charVector char_vals
    MPI_Type_get_extent(MPI_CHAR, &lb, &extent);
    encodedBlockColorOffsets[4] = encodedBlockColorOffsets[3] + extent;
    encodedBlockColorTypes[4] = MPI_CharVector;
    encodedBlockColorBlocks[4] = 1;
    // doubleVector double_vals
    MPI_Type_get_extent(MPI_CharVector, &lb, &extent);
    encodedBlockColorOffsets[5] = encodedBlockColorOffsets[4] + extent;
    encodedBlockColorTypes[5] = MPI_DoubleVector;
    encodedBlockColorBlocks[5] = 1;
    MPI_Type_create_struct(encodedBlockColorLen, encodedBlockColorBlocks,
        encodedBlockColorOffsets, encodedBlockColorTypes, &MPI_EncodedBlockColor);
    MPI_Type_commit(&MPI_EncodedBlockColor);

    // Set up EncodedBlock
    // Change to struct: change it to array of EncodedBlockColor[3] instead
    //      of y, cr, cb fields
    MPI_Datatype MPI_EncodedBlock;
    MPI_Type_contiguous(3, MPI_EncodedBlockColor, &MPI_EncodedBlock);
    MPI_Type_commit(&MPI_EncodedBlock);
    double mpiSetupEndTime = CycleTimer::currentSeconds();
    // End setup for MPI Structs

    int tag = 1;

    /*
     * BEGIN SCATTER
     * Purpose: give each thread a section of bytes to convert to an image
     */

    log(rank, "SCATTER\n");
    log(rank, "convertBytesToImage()...\n");
    double convertBytesToImageStartTime = CycleTimer::currentSeconds();
    int numPixels = bytes.size() / 4;
    int len = numPixels / numTasks;
    int start = rank * len * 4;
    int end = (rank + 1) * len * 4;
    // handle rounding issues
    if (rank == numTasks - 1) {
        end = bytes.size();
    }
    std::shared_ptr<ImageRgb> imageRgb = convertBytesToImage(bytes, width, height, start, end);
    double convertBytesToImageEndTime = CycleTimer::currentSeconds();

    // each thread will continue and convert its image to ycbcr
    log(rank, "convertRgbToYcbcr()...\n");
    double convertRgbToYcbcrStartTime = CycleTimer::currentSeconds();
    std::shared_ptr<ImageYcbcr> imageYcbcr = convertRgbToYcbcr(imageRgb);
    double convertRgbToYcbcrEndTime = CycleTimer::currentSeconds();

    /*
     * BEGIN GATHER
     * Purpose: collect all images into the full image in master
     */
    log(rank, "GATHER\n");
    double gatherYcbcrPixelsStartTime = CycleTimer::currentSeconds();
    if (rank == 0) {
        for (int i = 1; i < numTasks; i++) {
            start = i * len;
            end = (i + 1) * len;
            // handle rounding issues
            if (i == numTasks - 1) {
                end = bytes.size();
            }
            // receive image metadata
            MPI_Recv(&numPixels, 1, MPI_INT, i, MPI_ANY_TAG, MPI_COMM_WORLD, &mpiStatus);
            // receive image pixels
            std::shared_ptr<PixelYcbcr> buffer(new PixelYcbcr[numPixels]);
            MPI_Recv(buffer.get(), numPixels, MPI_PixelYcbcr, i, MPI_ANY_TAG, MPI_COMM_WORLD, &mpiStatus);
            // add pixels to full image
            for (int j = 0; j < numPixels; j++) {
                std::shared_ptr<PixelYcbcr> pixel(new PixelYcbcr());
                pixel->y = (buffer.get())[j].y;
                pixel->cb = (buffer.get())[j].cb;
                pixel->cr = (buffer.get())[j].cr;
                imageYcbcr->pixels.push_back(pixel);
                imageYcbcr->numPixels++;
            }
        }
    } else {
        MPI_Send(&(imageYcbcr->numPixels), 1, MPI_INT, 0, tag, MPI_COMM_WORLD);
        std::shared_ptr<PixelYcbcr> buffer(new PixelYcbcr[imageYcbcr->numPixels]);
        for (int i = 0; i < imageYcbcr->numPixels; i++) {
            (buffer.get())[i].y = imageYcbcr->pixels[i]->y;
            (buffer.get())[i].cb = imageYcbcr->pixels[i]->cb;
            (buffer.get())[i].cr = imageYcbcr->pixels[i]->cr;
        }
        MPI_Send(buffer.get(), imageYcbcr->numPixels, MPI_PixelYcbcr, 0, tag, MPI_COMM_WORLD);
    }
    double gatherYcbcrPixelsEndTime = CycleTimer::currentSeconds();

    /*
     * END GATHER
     */

    std::vector<std::vector<std::shared_ptr<PixelYcbcr>>> imageBlocksWorker;

    double convertYcbcrToBlocksStartTime = CycleTimer::currentSeconds();
    if (rank == 0) {
        log(rank, "convertYcbcrToBlocks()...\n");

        std::shared_ptr<ImageBlocks> imageBlocks = convertYcbcrToBlocks(imageYcbcr, MACROBLOCK_SIZE);

        /*
         * BEGIN SCATTER
         * Purpose: send blocks out to threads
         */
        log(rank, "SCATTER\n");

        // pixels in the blocks for each thread
        std::vector<std::shared_ptr<PixelYcbcr>> workerPixels[numTasks];

        // choose blocks for each thread
        for (unsigned int i = 0; i < imageBlocks->blocks.size(); i++) {
            // append the pixels in the block to the running list for the thread
            int index = i % numTasks;
            std::vector<std::shared_ptr<PixelYcbcr>> block = imageBlocks->blocks[i];
            workerPixels[index].insert(workerPixels[index].end(), block.begin(), block.end());
        }

        // send blocks to each thread
        for (int i = 1; i < numTasks; i++) {
            numPixels = workerPixels[i].size();
            std::shared_ptr<PixelYcbcr> buffer(new PixelYcbcr[numPixels]);
            for (int j = 0; j < numPixels; j++) {
                (buffer.get())[j].y = workerPixels[i][j]->y;
                (buffer.get())[j].cb = workerPixels[i][j]->cb;
                (buffer.get())[j].cr = workerPixels[i][j]->cr;
            }
            MPI_Send(&numPixels, 1, MPI_INT, i, tag, MPI_COMM_WORLD);
            MPI_Send(buffer.get(), numPixels, MPI_PixelYcbcr, i, tag, MPI_COMM_WORLD);
        }

        int pixelsPerBlock = MACROBLOCK_SIZE * MACROBLOCK_SIZE;
        // blocks for master thread
        std::vector<std::vector<std::shared_ptr<PixelYcbcr>>> blocks(workerPixels[0].size() / pixelsPerBlock);
        for (unsigned int i = 0; i < workerPixels[0].size(); i++) {
            int index = i / pixelsPerBlock;
            blocks[index].push_back(workerPixels[0][i]);
        }

        imageBlocksWorker = blocks;
    } else {
        // Get number of pixels to recv from master
        MPI_Recv(&numPixels, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &mpiStatus);
        std::shared_ptr<PixelYcbcr> buffer(new PixelYcbcr[numPixels]);
        MPI_Recv(buffer.get(), numPixels, MPI_PixelYcbcr, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &mpiStatus);

        // pixels per block
        int pixelsPerBlock = MACROBLOCK_SIZE * MACROBLOCK_SIZE;
        // blocks for master thread
        std::vector<std::vector<std::shared_ptr<PixelYcbcr>>> blocks(numPixels / pixelsPerBlock);
        for (int i = 0; i < numPixels; i++) {
            int index = i / pixelsPerBlock;
            std::shared_ptr<PixelYcbcr> pixel(new PixelYcbcr());
            pixel->y = (buffer.get())[i].y;
            pixel->cb = (buffer.get())[i].cb;
            pixel->cr = (buffer.get())[i].cr;
            blocks[index].push_back(pixel);
        }

        imageBlocksWorker = blocks;
    }
    double convertYcbcrToBlocksEndTime = CycleTimer::currentSeconds();

    /*
     * END SCATTER
     * Purpose: send blocks out to threads
     */

    // blocks stay in their threads
    log(rank, "DCT()...\n");
    double dctStartTime = CycleTimer::currentSeconds();
    std::vector<std::vector<std::shared_ptr<PixelYcbcr>>> dcts;
    for (auto block : imageBlocksWorker) {
        dcts.push_back(DCT(block, MACROBLOCK_SIZE, true));
    }
    double dctEndTime = CycleTimer::currentSeconds();

    // blocks stay in their threads
    log(rank, "quantize()...\n");
    double quantizeStartTime = CycleTimer::currentSeconds();
    std::vector<std::vector<std::shared_ptr<PixelYcbcr>>> quantizedBlocks;
    for (auto dct : dcts) {
        quantizedBlocks.push_back(quantize(dct, MACROBLOCK_SIZE, true));
    }
    double quantizeEndTime = CycleTimer::currentSeconds();

    // blocks stay in their threads
    log(rank, "DPCM()...\n");
    double dpcmStartTime = CycleTimer::currentSeconds();
    DPCM(quantizedBlocks);
    double dpcmEndTime = CycleTimer::currentSeconds();

    // blocks stay in their threads
    log(rank, "RLE()...\n");
    double encodedBlocksStartTime = CycleTimer::currentSeconds();
    std::vector<std::shared_ptr<EncodedBlock>> encodedBlocks;
    for (auto quantizedBlock : quantizedBlocks) {
        encodedBlocks.push_back(RLE(quantizedBlock, MACROBLOCK_SIZE));
    }
    double encodedBlocksEndTime = CycleTimer::currentSeconds();

    /*
     * BEGIN GATHER
     * Purpose: collect all blocks back into a vector of encoded blocks in master
     */
    log(rank, "GATHER\n");
    double gatherEncodedBlocksStartTime = CycleTimer::currentSeconds();
    // For master thread, collect encoded blocks in order from workers
    int numEncodedBlocks;
    // encoded blocks for all threads
    std::vector<std::shared_ptr<EncodedBlock>> allEncodedBlocks[numTasks];
    // encoded blocks in original order
    std::vector<std::shared_ptr<EncodedBlock>> finalEncodedBlocks;
    if (rank == 0) {
        allEncodedBlocks[0] = encodedBlocks;
        // grab all encoded blocks from workers
        for (int i = 1; i < numTasks; i++) {
            // recv the number of encoded blocks
            MPI_Recv(&numEncodedBlocks, 1, MPI_INT, i, MPI_ANY_TAG, MPI_COMM_WORLD, &mpiStatus);
            // recv the encoded blocks
            std::shared_ptr<EncodedBlockNoPtr> encodedBlocksBuffer(new EncodedBlockNoPtr[numEncodedBlocks]);
            MPI_Recv(encodedBlocksBuffer.get(), numEncodedBlocks, MPI_EncodedBlock, i, MPI_ANY_TAG, MPI_COMM_WORLD, &mpiStatus);
            std::vector<std::shared_ptr<EncodedBlock>> encodedBlocksRecv = convertBufferToEncodedBlocks(encodedBlocksBuffer, numEncodedBlocks);
            std::vector<std::shared_ptr<EncodedBlock>> threadEncodedBlocks;
            for (int j = 0; j < numEncodedBlocks; j++) {
                threadEncodedBlocks.push_back(encodedBlocksRecv[j]);
            }
            allEncodedBlocks[i] = threadEncodedBlocks;
        }
        // find number of total encoded blocks
        int totalEncodedBlocks = 0;
        int indices[numTasks];
        for (int i = 0; i < numTasks; i++) {
            totalEncodedBlocks += allEncodedBlocks[i].size();
            indices[i] = 0;
        }
        // reorder encoded blocks by original order
        for (int i = 0; i < totalEncodedBlocks; i++) {
            int index = i % numTasks;
            finalEncodedBlocks.push_back(allEncodedBlocks[index][indices[index]]);
            indices[index]++;
        }
    } else {
        // For each worker, send its encoded blocks back to master
        numEncodedBlocks = encodedBlocks.size();
        // Send number of encoded blocks
        MPI_Send(&numEncodedBlocks, 1, MPI_INT, 0, tag, MPI_COMM_WORLD);
        std::shared_ptr<EncodedBlockNoPtr> encodedBlockBuffer(new EncodedBlockNoPtr[numEncodedBlocks]);
        // Build the encoded blocks structure
        for (unsigned int i = 0; i < encodedBlocks.size(); i++) {
            writeToBuffer(encodedBlockBuffer, encodedBlocks, i, COLOR_Y);
            writeToBuffer(encodedBlockBuffer, encodedBlocks, i, COLOR_CR);
            writeToBuffer(encodedBlockBuffer, encodedBlocks, i, COLOR_CB);
        }
        // Send the encoded blocks
        MPI_Send(encodedBlockBuffer.get(), numEncodedBlocks, MPI_EncodedBlock, 0, tag, MPI_COMM_WORLD);
        MPI_Finalize();
        return;
    }
    double gatherEncodedBlocksEndTime = CycleTimer::currentSeconds();

    /*
     * END GATHER
     * Purpose: collect all blocks back into a vector of encoded blocks in master
     */

    double encodeCompressedStartTime = CycleTimer::currentSeconds();
    if (rank == 0) {
        log(0, "done encoding!\n");
        log(0, "writing to file...\n");
        std::ofstream jpegFile(compressedFile);
        for (const auto &block : finalEncodedBlocks) {
            jpegFile << block;
        }
        log(0, "jpeg stored!\n");
        endTime = CycleTimer::currentSeconds();
        log(0, "Time Elapsed: %.3fs\n", (endTime - startTime));
        log(0, "==============\n");
        log(0, "now let's undo the process...\n");
    }

    log(rank, "undoing RLE()...\n");
    std::vector<std::vector<std::shared_ptr<PixelYcbcr>>> decodedQuantizedBlocks;
    for (auto encodedBlock : finalEncodedBlocks) {
        decodedQuantizedBlocks.push_back(decodeRLE(encodedBlock, MACROBLOCK_SIZE));
    }

    log(rank, "undoing DPCM()...\n");
    unDPCM(decodedQuantizedBlocks);

    log(rank, "undoing quantize()...\n");
    std::vector<std::vector<std::shared_ptr<PixelYcbcr>>> unquantizedBlocks;
    for (auto decodedQuantizedBlock : decodedQuantizedBlocks) {
        unquantizedBlocks.push_back(unquantize(decodedQuantizedBlock, MACROBLOCK_SIZE, true));
    }

    log(rank, "undoing DCT()...\n");
    std::vector<std::vector<std::shared_ptr<PixelYcbcr>>> idcts;
    for (auto unquantized : unquantizedBlocks) {
        idcts.push_back(IDCT(unquantized, MACROBLOCK_SIZE, true));
    }

    log(rank, "undoing convertYcbcrToBlocks()...\n");
    std::shared_ptr<ImageBlocks> imageBlocksIdct(new ImageBlocks);
    imageBlocksIdct->blocks = idcts;
    imageBlocksIdct->width = width;
    imageBlocksIdct->height = height;
    std::shared_ptr<ImageYcbcr> imgFromBlocks = convertBlocksToYcbcr(imageBlocksIdct, MACROBLOCK_SIZE);
    log(rank, "undoing convertRgbToYcbcr()...\n");
    std::shared_ptr<ImageRgb> imageRgbRecovered = convertYcbcrToRgb(imgFromBlocks);

    log(rank, "undoing convertImageToBytes()...\n");
    std::vector<unsigned char> imgRecovered = convertImageToBytes(imageRgbRecovered);

    if (rank == 0) {
        error = lodepng::encode(outfile, imgRecovered, width, height);

        if(error) {
            std::cout << "encoder error " << error << ": "<< lodepng_error_text(error) << std::endl;
        } else {
            fprintf(stdout, "success encoding to %s!\n", outfile);
        }
    }

    // Free derived types
    MPI_Type_free(&MPI_PixelYcbcr);
    MPI_Type_free(&MPI_RleTuple);
    MPI_Type_free(&MPI_RleTupleVector);
    MPI_Type_free(&MPI_CharVector);
    MPI_Type_free(&MPI_DoubleVector);
    MPI_Type_free(&MPI_EncodedBlockColor);
    MPI_Type_free(&MPI_EncodedBlock);

    // End parallel area
    MPI_Finalize();

    // Print statistics
    if (rank == 0) {
        fprintf(stdout,
        "=======================================\n"
        "= MPI encoding performance: \n"
        "=======================================\n"
        "Load image: %.3fs\n"
        "Setup MPI: %.3fs\n"
        "Convert Bytes to Image: %.3fs\n"
        "Convert RGB to YCbCr: %.3fs\n"
        "Gather YCbCr Pixels: %.3fs\n"
        "Convert YCbCr to Blocks: %.3fs\n"
        "DCT: %.3fs\n"
        "Quantize: %.3fs\n"
        "DPCM: %.3fs\n"
        "RLE: %.3fs\n"
        "Gather Encoded Blocks: %.3fs\n"
        "Encode Compressed Image: %.3fs\n"
        "Total time: %.3fs\n",
        loadImageEndTime - loadImageStartTime,
        mpiSetupEndTime - mpiSetupStartTime,
        convertBytesToImageEndTime - convertBytesToImageStartTime,
        convertRgbToYcbcrEndTime - convertRgbToYcbcrStartTime,
        gatherYcbcrPixelsEndTime - gatherYcbcrPixelsStartTime,
        convertYcbcrToBlocksEndTime - convertYcbcrToBlocksStartTime,
        dctEndTime - dctStartTime,
        quantizeEndTime - quantizeStartTime,
        dpcmEndTime - dpcmStartTime,
        encodedBlocksEndTime - encodedBlocksStartTime,
        gatherEncodedBlocksEndTime - gatherEncodedBlocksStartTime,
        endTime - encodeCompressedStartTime,
        endTime - startTime);
    }
}

int main(int argc, char** argv) {
    std::string filename = argv[1];
    int opt;
    int mpi = 0;
    while ((opt = getopt(argc, argv, ":p")) != -1) {
        switch (opt) {
            case 'p':
                mpi = 1;
                break;
            default:
                fprintf(stderr, "Usage: %s [image] [-p]\n", argv[0]);
                exit(EXIT_FAILURE);
        }
    }

    std::string raw_image = std::string("raw_images/") + filename + std::string(".png");
    std::string image = std::string("images/") + filename + std::string(".png");
    std::string compressed = std::string("compressed/") + filename + std::string(".jpeg");

    if (mpi) {
        encodeMpi(raw_image.c_str(), image.c_str(), compressed.c_str());
    } else {
        encodeSeq(raw_image.c_str(), image.c_str(), compressed.c_str());
    }

    exit(EXIT_SUCCESS);
}
