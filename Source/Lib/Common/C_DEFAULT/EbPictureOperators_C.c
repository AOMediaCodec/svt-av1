/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#include "EbPictureOperators_C.h"
#include "EbUtility.h"

/*********************************
* Picture Average
*********************************/
void PictureAverageKernel(
    EbByte   src0,
    uint32_t   src0Stride,
    EbByte   src1,
    uint32_t   src1Stride,
    EbByte   dst,
    uint32_t   dstStride,
    uint32_t   areaWidth,
    uint32_t   areaHeight)
{
    uint32_t x, y;

    for (y = 0; y < areaHeight; y++) {
        for (x = 0; x < areaWidth; x++) {
            dst[x] = (src0[x] + src1[x] + 1) >> 1;
        }
        src0 += src0Stride;
        src1 += src1Stride;
        dst += dstStride;
    }
}

/*********************************
* Picture Copy Kernel
*********************************/
void picture_copy_kernel(
    EbByte                  src,
    uint32_t                   src_stride,
    EbByte                  dst,
    uint32_t                   dst_stride,
    uint32_t                   area_width,
    uint32_t                   area_height,
    uint32_t                   bytes_per_sample)  //=1 always)
{
    uint32_t sampleCount = 0;
    const uint32_t sampleTotalCount = area_width * area_height;
    const uint32_t copyLength = area_width * bytes_per_sample;

    src_stride *= bytes_per_sample;
    dst_stride *= bytes_per_sample;

    while (sampleCount < sampleTotalCount) {
        EB_MEMCPY(dst, src, copyLength);
        src += src_stride;
        dst += dst_stride;
        sampleCount += area_width;
    }

    return;
}

// C equivalents

uint64_t spatial_full_distortion_kernel_c(
    uint8_t   *input,
    uint32_t   input_offset,
    uint32_t   input_stride,
    uint8_t   *recon,
    uint32_t   recon_offset,
    uint32_t   recon_stride,
    uint32_t   area_width,
    uint32_t   area_height)
{
    uint32_t  columnIndex;
    uint32_t  row_index = 0;
    uint64_t  spatialDistortion = 0;
    input += input_offset;
    recon += recon_offset;

    while (row_index < area_height) {
        columnIndex = 0;
        while (columnIndex < area_width) {
            spatialDistortion += (int64_t)SQR((int64_t)(input[columnIndex]) - (recon[columnIndex]));
            ++columnIndex;
        }

        input += input_stride;
        recon += recon_stride;
        ++row_index;
    }
    return spatialDistortion;
}
