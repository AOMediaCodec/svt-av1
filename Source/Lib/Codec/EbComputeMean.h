/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbComputeMean_h
#define EbComputeMean_h

#include "EbComputeMean_SSE2.h"
#include "EbCombinedAveragingSAD_Intrinsic_AVX2.h"
#include "EbDefinitions.h"

#ifdef __cplusplus
extern "C" {
#endif


    typedef uint64_t(*EB_COMPUTE_MEAN_FUNC)(
        uint8_t *input_samples,
        uint32_t input_stride,
        uint32_t input_area_width,
        uint32_t input_area_height);

    static const EB_COMPUTE_MEAN_FUNC ComputeMeanFunc[2][ASM_TYPE_TOTAL] = {
        {
            // NON_AVX2
            ComputeMean8x8_SSE2_INTRIN,
            // AVX2
            compute_mean8x8_avx2_intrin
        },
        {
            // NON_AVX2
            ComputeMeanOfSquaredValues8x8_SSE2_INTRIN,
            // AVX2
            ComputeMeanOfSquaredValues8x8_SSE2_INTRIN
        }
    };

#ifdef __cplusplus
}
#endif

#endif