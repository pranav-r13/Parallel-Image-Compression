#include <vector>
#include <memory>
#include "image.h"

void DPCM(std::vector<std::vector<std::shared_ptr<PixelYcbcr>>> blocks);
void unDPCM(std::vector<std::vector<std::shared_ptr<PixelYcbcr>>> blocks);
