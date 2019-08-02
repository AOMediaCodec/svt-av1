/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbNoiseExtractAVX2_h
#define EbNoiseExtractAVX2_h

#include "immintrin.h"
#include "EbDefinitions.h"
#include "EbPictureBufferDesc.h"
#ifdef __cplusplus
extern "C" {
#endif

    /*******************************************
    * noise_extract_luma_weak
    *  weak filter Luma and store noise.
    *******************************************/
       void chroma_strong_avx2_intrin(
        __m256i   top,
        __m256i   curr,
        __m256i   bottom,
        __m256i   curr_prev,
        __m256i   curr_next,
        __m256i   top_prev,
        __m256i   top_next,
        __m256i   bottom_prev,
        __m256i   bottom_next,
        uint8_t  *ptr_denoised);

    void luma_weak_filter_avx2_intrin(
        __m256i  top,
        __m256i  curr,
        __m256i  bottom,
        __m256i  curr_prev,
        __m256i  curr_next,
        uint8_t *ptr_denoised,
        uint8_t *ptr_noise);

    void chroma_weak_luma_strong_filter_avx2_intrin(
        __m256i  top,
        __m256i  curr,
        __m256i  bottom,
        __m256i  curr_prev,
        __m256i  curr_next,
        __m256i  top_prev,
        __m256i  top_next,
        __m256i  bottom_prev,
        __m256i  bottom_next,
        uint8_t *ptr_denoised);

#ifdef __cplusplus
}
#endif
#endif // EbNoiseExtractAVX2_h
