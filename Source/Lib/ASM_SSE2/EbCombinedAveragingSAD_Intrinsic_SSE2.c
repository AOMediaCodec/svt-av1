/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#include "stdint.h"
#include "emmintrin.h"
#include "EbComputeSAD_SSE2.h"

uint32_t Compute4xMSadSub_SSE2_INTRIN(
    uint8_t  *src,                            // input parameter, source samples Ptr
    uint32_t  src_stride,                      // input parameter, source stride
    uint8_t  *ref,                            // input parameter, reference samples Ptr
    uint32_t  ref_stride,                      // input parameter, reference stride
    uint32_t  height,                         // input parameter, block height (M)
    uint32_t  width)                          // input parameter, block width (N)
{
    uint32_t y;
    (void)width;
    __m128i xmm_sad = _mm_setzero_si128();

    // for (y = 0; y < height; y+=4) {
    for (y = 0; y < height; y++) {

        xmm_sad = _mm_add_epi16(xmm_sad, _mm_sad_epu8(_mm_cvtsi32_si128(*(uint32_t*)src), _mm_cvtsi32_si128(*(uint32_t*)ref)));
        //  xmm_sad = _mm_add_epi16(xmm_sad, _mm_sad_epu8(_mm_cvtsi32_si128(*(uint32_t*)(src+src_stride)), _mm_cvtsi32_si128(*(uint32_t*)(ref+ref_stride))));
        //  xmm_sad = _mm_add_epi16(xmm_sad, _mm_sad_epu8(_mm_cvtsi32_si128(*(uint32_t*)(src+(src_stride << 1))), _mm_cvtsi32_si128(*(uint32_t*)(ref+(ref_stride << 1)))));
        //  xmm_sad = _mm_add_epi16(xmm_sad, _mm_sad_epu8(_mm_cvtsi32_si128(*(uint32_t*)(src+3*src_stride)), _mm_cvtsi32_si128(*(uint32_t*)(ref+ 3*ref_stride))));

       //   src += (src_stride << 2);
        //  ref += (ref_stride << 2);
        src += src_stride;
        ref += ref_stride;
    }
    return _mm_cvtsi128_si32(xmm_sad);
}

uint32_t CombinedAveraging4xMSAD_SSE2_INTRIN(
    uint8_t  *src,
    uint32_t  src_stride,
    uint8_t  *ref1,
    uint32_t  ref1_stride,
    uint8_t  *ref2,
    uint32_t  ref2_stride,
    uint32_t  height,
    uint32_t  width)
{
    __m128i sad0, sad1;
    uint32_t y;
    (void)width;
    sad0 = sad1 = _mm_setzero_si128();

    for (y = 0; y < height; y += 2) {

        sad0 = _mm_add_epi32(sad0, _mm_sad_epu8(_mm_cvtsi32_si128(*(uint32_t *)src), _mm_avg_epu8(_mm_cvtsi32_si128(*(uint32_t *)ref1), _mm_cvtsi32_si128(*(uint32_t *)ref2))));

        sad1 = _mm_add_epi32(sad1, _mm_sad_epu8(_mm_cvtsi32_si128(*(uint32_t *)(src + src_stride)), _mm_avg_epu8(_mm_cvtsi32_si128(*(uint32_t *)(ref1 + ref1_stride)), _mm_cvtsi32_si128(*(uint32_t *)(ref2 + ref2_stride)))));
        src += src_stride << 1;
        ref1 += ref1_stride << 1;
        ref2 += ref2_stride << 1;
    }
    return _mm_cvtsi128_si32(_mm_add_epi32(sad0, sad1));
}
