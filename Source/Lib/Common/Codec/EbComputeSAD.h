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

    /***************************************
    * Function Ptr Types
    ***************************************/

    static void nxm_sad_kernel_void_func() {}

    typedef void(*EbSadLoopKernelNxMType)(
        uint8_t  *src,                            // input parameter, source samples Ptr
        uint32_t  src_stride,                      // input parameter, source stride
        uint8_t  *ref,                            // input parameter, reference samples Ptr
        uint32_t  ref_stride,                      // input parameter, reference stride
        uint32_t  height,                         // input parameter, block height (M)
        uint32_t  width,                          // input parameter, block width (N)
        uint64_t *best_sad,
        int16_t *x_search_center,
        int16_t *y_search_center,
        uint32_t  src_stride_raw,                   // input parameter, source stride (no line skipping)
        int16_t search_area_width,
        int16_t search_area_height);

    typedef uint32_t(*EbSadAvgKernelNxMType)(
        uint8_t  *src,
        uint32_t  src_stride,
        uint8_t  *ref1,
        uint32_t  ref1_stride,
        uint8_t  *ref2,
        uint32_t  ref2_stride,
        uint32_t  height,
        uint32_t  width);

    typedef uint32_t(*EbCompute8x4SadType)(
        uint8_t  *src,                            // input parameter, source samples Ptr
        uint32_t  src_stride,                      // input parameter, source stride
        uint8_t  *ref,                            // input parameter, reference samples Ptr
        uint32_t  ref_stride);                     // input parameter, reference stride

    typedef uint32_t(*EB_COMPUTE8X8SAD_TYPE)(
        uint8_t  *src,                            // input parameter, source samples Ptr
        uint32_t  src_stride,                      // input parameter, source stride
        uint8_t  *ref,                            // input parameter, reference samples Ptr
        uint32_t  ref_stride);                     // input parameter, reference stride
    typedef void(*EbGetEightSad8x8)(
        uint8_t   *src,
        uint32_t   src_stride,
        uint8_t   *ref,
        uint32_t   ref_stride,
        uint32_t  *p_best_sad8x8,
        uint32_t  *p_best_mv8x8,
        uint32_t  *p_best_sad16x16,
        uint32_t  *p_best_mv16x16,
        uint32_t   mv,
        uint16_t  *p_sad16x16,
        EbBool     sub_sad);

    typedef void(*EbGetEightSad32x32)(
        uint16_t  *p_sad16x16,
        uint32_t  *p_best_sad32x32,
        uint32_t  *p_best_sad64x64,
        uint32_t  *p_best_mv32x32,
        uint32_t  *p_best_mv64x64,
        uint32_t   mv);

    typedef uint32_t(*CombinedAveragingSsd)(
        uint8_t   *src,
        ptrdiff_t  src_stride,
        uint8_t   *ref1,
        ptrdiff_t  ref1_stride,
        uint8_t   *ref2,
        ptrdiff_t  ref2_stride,
        uint32_t   height,
        uint32_t   width
        );

    static EbSadAvgKernelNxMType FUNC_TABLE nxm_sad_averaging_kernel_func_ptr_array[ASM_TYPE_TOTAL][9] =   // [asm_type][SAD - block height]
    {
        // NON_AVX2
        {
            /*0 4xM  */     combined_averaging_sad,
            /*1 8xM  */     combined_averaging_sad,
            /*2 16xM */     combined_averaging_sad,
            /*3 24xM */     combined_averaging_sad,
            /*4 32xM */     combined_averaging_sad,
            /*5      */     (EbSadAvgKernelNxMType)nxm_sad_kernel_void_func,
            /*6 48xM */     combined_averaging_sad,
            /*7      */     (EbSadAvgKernelNxMType)nxm_sad_kernel_void_func,
            /*8 64xM */     combined_averaging_sad
        },
        // AVX2
        {
            /*0 4xM  */     combined_averaging4x_msad_sse2_intrin,
            /*1 8xM  */     combined_averaging8x_msad_avx2_intrin,
            /*2 16xM */     combined_averaging16x_msad_avx2_intrin,
            /*3 24xM */     combined_averaging24x_msad_avx2_intrin,
            /*4 32xM */     combined_averaging32x_msad_avx2_intrin,
            /*5      */     (EbSadAvgKernelNxMType)nxm_sad_kernel_void_func,
            /*6 48xM */     combined_averaging48x_msad_avx2_intrin,
            /*7      */     (EbSadAvgKernelNxMType)nxm_sad_kernel_void_func,
            /*8 64xM */     combined_averaging64x_msad_avx2_intrin
        },
    };

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
