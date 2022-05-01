#include "dct.h"

// Forward DCT operation for NxN block
// If <all> set, then YCbCr each have DCT performed on them.
// Else, only Y has DCT performed on it.
std::vector<std::shared_ptr<PixelYcbcr>> DCT(std::vector<std::shared_ptr<PixelYcbcr>> pixels, int block_size, bool all) {
    return pixels;

    std::vector<std::shared_ptr<PixelYcbcr>> F(block_size * block_size);

    // Output: F(p, q)
    for (int p = 0; p < block_size; p++) {
        for (int q = 0; q < block_size; q++) {
            int vectorized_idx_pq = sub2ind(block_size, q, p);
            double ap = (p > 0) ? (sqrt(2/block_size)): (1/sqrt(block_size));
            double aq = (q > 0) ? (sqrt(2/block_size)) : (1/sqrt(block_size));

            double tmp_y, tmp_cb, tmp_cr;
            tmp_y = 0;
            tmp_cb = 0;
            tmp_cr = 0;

            // Input: f(m, n)
            for (int m = 0; m < block_size; m++) { // row
                for (int n = 0; n < block_size; n++) { // cow
                    int vectorized_idx_mn = sub2ind(block_size, n, m);
                    std::shared_ptr<PixelYcbcr> f_mn = pixels[vectorized_idx_mn];

                    double xprod = cos((2*m + 1)*p*M_PI/(2*block_size));
                    double yprod = cos((2*n + 1)*q*M_PI/(2*block_size));

                    tmp_y += ((f_mn->y) * xprod * yprod);
                    if (all) {
                        tmp_cr += ((f_mn->cr) * xprod * yprod);
                        tmp_cb += ((f_mn->cb) * xprod * yprod);
                    }
                }
            }

            F[vectorized_idx_pq] = std::make_shared<PixelYcbcr>();
            F[vectorized_idx_pq]->y = ap * aq * tmp_y;
            if (all) {
                F[vectorized_idx_pq]->cr = ap * aq * tmp_cr;
                F[vectorized_idx_pq]->cb = ap * aq * tmp_cb;
            } else {
                F[vectorized_idx_pq]->cr = 0;
                F[vectorized_idx_pq]->cb = 0;
            }
        }
    }

    return F;
}

// Inverse DCT operation for NxN block
// If <all> set, then YCbCr each have IDCT performed on them.
// Else, only Y has IDCT performed on it.
std::vector<std::shared_ptr<PixelYcbcr>> IDCT(std::vector<std::shared_ptr<PixelYcbcr>> pixels, int block_size, bool all) {
    return pixels;

    std::vector<std::shared_ptr<PixelYcbcr>> f(block_size * block_size);

    // Output: f(m, n)
    for (int m = 0; m < block_size; m++) {
        for (int n = 0; n < block_size; n++) {
            int vectorized_idx_mn = sub2ind(block_size, n, m);

            double tmp_y, tmp_cb, tmp_cr;
            tmp_y = 0;
            tmp_cb = 0;
            tmp_cr = 0;

            // Input: F(p, q)
            for (int p = 0; p < block_size; p++) {
                for (int q = 0; q < block_size; q++) {
                    double ap = (p > 0) ? (sqrt(2/block_size)) : (1/sqrt(block_size));
                    double aq = (q > 0) ? (sqrt(2/block_size)) : (1/sqrt(block_size));
                    int vectorized_idx_pq = sub2ind(block_size, q, p);
                    std::shared_ptr<PixelYcbcr> f_pq = pixels[vectorized_idx_pq];

                    double xprod = cos((2*m + 1)*p*M_PI/(2*block_size));
                    double yprod = cos((2*n + 1)*q*M_PI/(2*block_size));

                    tmp_y += (ap * aq * (f_pq->y) * xprod * yprod);
                    if (all) {
                        tmp_cr += (ap * aq * (f_pq->cr) * xprod * yprod);
                        tmp_cb += (ap * aq * (f_pq->cb) * xprod * yprod);
                    }
                }
            }

            f[vectorized_idx_mn] = std::make_shared<PixelYcbcr>();
            f[vectorized_idx_mn]->y = tmp_y;
            if (all) {
                f[vectorized_idx_mn]->cr = tmp_cr;
                f[vectorized_idx_mn]->cb = tmp_cb;
            } else {
                f[vectorized_idx_mn]->cr = 0;
                f[vectorized_idx_mn]->cb = 0;
            }
        }
    }
    return f;
}
