/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbComputeSAD_h
#define EbComputeSAD_h

#include "EbDefinitions.h"

#include "EbCombinedAveragingSAD_Intrinsic_AVX2.h"
#include "EbComputeSAD_C.h"
#include "EbComputeSAD_SSE2.h"
#include "EbComputeSAD_SSE4_1.h"
#include "EbComputeSAD_AVX2.h"
#include "EbUtility.h"
#ifdef __cplusplus
extern "C" {
#endif

uint32_t compute4x_m_sad_avx2_intrin(const uint8_t  *src, uint32_t  src_stride, const uint8_t  *ref, uint32_t  ref_stride, uint32_t  height, uint32_t  width);

uint32_t sad_16b_kernel(
    uint16_t  *src,                           // input parameter, source samples Ptr
    uint32_t  src_stride,                     // input parameter, source stride
    uint16_t  *ref,                           // input parameter, reference samples Ptr
    uint32_t  ref_stride,                     // input parameter, reference stride
    uint32_t  height,                         // input parameter, block height (M)
    uint32_t  width);                         // input parameter, block width (N)

#ifdef __cplusplus
}
#endif
#endif // EbComputeSAD_h
