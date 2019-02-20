/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EBAVCSTYLEMCP_H
#define EBAVCSTYLEMCP_H

#include "EbAvcStyleMcp_SSE2.h"
#include "EbAvcStyleMcp_SSSE3.h"

#include "EbPictureOperators.h"

#include "EbMcp.h"

#include "EbDefinitions.h"
#include "EbPictureBufferDesc.h"
#include "EbPictureControlSet.h"

#ifdef __cplusplus
extern "C" {
#endif

    void estimate_bi_pred_interpolation_avc_luma(
        EbPictureBufferDesc_t *ref_pic_list0,
        EbPictureBufferDesc_t *ref_pic_list1,
        uint32_t                 refList0PosX,
        uint32_t                 refList0PosY,
        uint32_t                 refList1PosX,
        uint32_t                 refList1PosY,
        uint32_t                 pu_width,
        uint32_t                 pu_height,
        EbPictureBufferDesc_t *biDst,
        uint32_t                 dstLumaIndex,
        uint32_t                 dstChromaIndex,
        uint32_t                 component_mask,
        EbByte                refList0TempDst,
        EbByte                refList1TempDst,
        EbByte                firstPassIFTempDst,
        EbBool                subSamplePredFlag,
        EbAsm                 asm_type);

    void estimate_uni_pred_interpolation_avc_luma(
        EbPictureBufferDesc_t *ref_pic,
        uint32_t                 posX,
        uint32_t                 posY,
        uint32_t                 pu_width,
        uint32_t                 pu_height,
        EbPictureBufferDesc_t *dst,
        uint32_t                 dstLumaIndex,
        uint32_t                 dstChromaIndex,
        uint32_t                 component_mask,
        EbByte                temp_buf,
        EbBool                subSamplePredFlag,
        EbAsm                 asm_type);

    void estimate_bi_pred_interpolation_unpacked_avc_style(
        EbPictureBufferDesc_t *ref_pic_list0,
        EbPictureBufferDesc_t *ref_pic_list1,
        uint32_t                 refList0PosX,
        uint32_t                 refList0PosY,
        uint32_t                 refList1PosX,
        uint32_t                 refList1PosY,
        uint32_t                 pu_width,
        uint32_t                 pu_height,
        EbPictureBufferDesc_t *biDst,
        uint32_t                 dstLumaIndex,
        uint32_t                 dstChromaIndex,
        uint32_t                 component_mask,
        EbByte                refList0TempDst,
        EbByte                refList1TempDst,
        EbByte                firstPassIFTempDst,
        EbBool                subSamplePredFlag,
        EbAsm                 asm_type);

    void estimate_uni_pred_interpolation_unpacked_avc_style(
        EbPictureBufferDesc_t *ref_pic,
        uint32_t                 posX,
        uint32_t                 posY,
        uint32_t                 pu_width,
        uint32_t                 pu_height,
        EbPictureBufferDesc_t *dst,
        uint32_t                 dstLumaIndex,
        uint32_t                 dstChromaIndex,
        uint32_t                 component_mask,
        EbByte                temp_buf,
        EbBool                subSamplePredFlag,
        EbAsm                 asm_type);

    void estimate_uni_pred_interpolation_avc_chroma_ref10_bit(
        EbPictureBufferDesc_t *refFramePicList0,
        uint32_t                 posX,
        uint32_t                 posY,
        uint32_t                 pu_width,
        uint32_t                 pu_height,
        EbPictureBufferDesc_t *dst,
        uint32_t                 dstLumaIndex,
        uint32_t                 dstChromaIndex,          //input parameter, please refer to the detailed explanation above.
        uint32_t                 component_mask,
        EbByte                temp_buf,
        EbBool                sub_pred,
        EbAsm                 asm_type);
    void estimate_bi_pred_interpolation_avc_chroma_ref10_bit(
        EbPictureBufferDesc_t *refFramePicList0,
        EbPictureBufferDesc_t *refFramePicList1,
        uint32_t                 refList0PosX,
        uint32_t                 refList0PosY,
        uint32_t                 refList1PosX,
        uint32_t                 refList1PosY,
        uint32_t                 pu_width,
        uint32_t                 pu_height,
        EbPictureBufferDesc_t *biDst,
        uint32_t                 dstLumaIndex,
        uint32_t                 dstChromaIndex,
        uint32_t                 component_mask,
        EbByte                refList0TempDst,
        EbByte                refList1TempDst,
        EbByte                firstPassIFTempDst,
        EbBool                sub_pred,
        EbAsm                 asm_type);

    void estimate_uni_pred_interpolation_avc_lumaRef10Bit(
        EbPictureBufferDesc_t *refFramePicList0,
        uint32_t                 posX,
        uint32_t                 posY,
        uint32_t                 pu_width,
        uint32_t                 pu_height,
        EbPictureBufferDesc_t *dst,
        uint32_t                 dstLumaIndex,
        uint32_t                 dstChromaIndex,          //input parameter, please refer to the detailed explanation above.
        uint32_t                 component_mask,
        EbByte                temp_buf,
        EbBool                sub_pred,
        EbBool                subPredChroma,
        EbAsm                 asm_type);

    void estimate_bi_pred_interpolation_avc_luma_ref10_bit(
        EbPictureBufferDesc_t *refFramePicList0,
        EbPictureBufferDesc_t *refFramePicList1,
        uint32_t                 refList0PosX,
        uint32_t                 refList0PosY,
        uint32_t                 refList1PosX,
        uint32_t                 refList1PosY,
        uint32_t                 pu_width,
        uint32_t                 pu_height,
        EbPictureBufferDesc_t *biDst,
        uint32_t                 dstLumaIndex,
        uint32_t                 dstChromaIndex,
        uint32_t                 component_mask,
        EbByte                refList0TempDst,
        EbByte                refList1TempDst,
        EbByte                firstPassIFTempDst,
        EbBool                sub_pred,
        EbBool                subPredChroma,
        EbAsm                 asm_type);

    void uni_pred_i_free_ref8_bit(
        EbPictureBufferDesc_t *ref_pic,
        uint32_t                 posX,
        uint32_t                 posY,
        uint32_t                 pu_width,
        uint32_t                 pu_height,
        EbPictureBufferDesc_t *dst,
        uint32_t                 dstLumaIndex,
        uint32_t                 dstChromaIndex,          //input parameter, please refer to the detailed explanation above.
        uint32_t                 component_mask,
        EbByte                temp_buf,
        EbBool                subSamplePredFlag,
        EbBool                subSamplePredFlagChroma,
        EbAsm                 asm_type);
    void bi_pred_i_free_ref8_bit(
        EbPictureBufferDesc_t *ref_pic_list0,
        EbPictureBufferDesc_t *ref_pic_list1,
        uint32_t                 refList0PosX,
        uint32_t                 refList0PosY,
        uint32_t                 refList1PosX,
        uint32_t                 refList1PosY,
        uint32_t                 pu_width,
        uint32_t                 pu_height,
        EbPictureBufferDesc_t *biDst,
        uint32_t                 dstLumaIndex,
        uint32_t                 dstChromaIndex,
        uint32_t                 component_mask,
        EbByte                refList0TempDst,
        EbByte                refList1TempDst,
        EbByte                firstPassIFTempDst,
        EbBool                subSamplePredFlag,
        EbBool                subSamplePredFlagChroma,
        EbAsm                 asm_type);


    typedef void(*AvcStyleInterpolationFilterNew)(
        EbByte               ref_pic,
        uint32_t                src_stride,
        EbByte               dst,
        uint32_t                dst_stride,
        uint32_t                pu_width,
        uint32_t                pu_height,
        EbByte               temp_buf,
        EbBool               skip,
        uint32_t                frac_pos);

    typedef void(*AvcStyleChromaInterpolationFilterNew)(
        EbByte               ref_pic,
        uint32_t                src_stride,
        EbByte               dst,
        uint32_t                dst_stride,
        uint32_t                pu_width,
        uint32_t                pu_height,
        EbByte               temp_buf,
        EbBool               skip,
        uint32_t                frac_pos_x,
        uint32_t                frac_pos_y);

    typedef void(*PictureAverage)(
        EbByte                  src0,
        uint32_t                   src0_stride,
        EbByte                  src1,
        uint32_t                   src1_stride,
        EbByte                  dst,
        uint32_t                   dst_stride,
        uint32_t                   area_width,
        uint32_t                   area_height);


    /***************************************
    * Function Tables
    ***************************************/
    static const AvcStyleInterpolationFilterNew FUNC_TABLE avc_style_uni_pred_luma_if_function_ptr_array[ASM_TYPE_TOTAL][16] = {
        // NON_AVX2
        {
            avc_style_copy_sse2,                                    //A
            AvcStyleLumaInterpolationFilterHorizontal_SSSE3_INTRIN,       //a
            AvcStyleLumaInterpolationFilterHorizontal_SSSE3_INTRIN,       //b
            AvcStyleLumaInterpolationFilterHorizontal_SSSE3_INTRIN,       //c
            AvcStyleLumaInterpolationFilterVertical_SSSE3_INTRIN,         //d
            AvcStyleLumaInterpolationFilterPose_SSSE3,             //e
            AvcStyleLumaInterpolationFilterPosf_SSSE3,             //f
            AvcStyleLumaInterpolationFilterPosg_SSSE3,             //g
            AvcStyleLumaInterpolationFilterVertical_SSSE3_INTRIN,         //h
            AvcStyleLumaInterpolationFilterPosi_SSSE3,             //i
            AvcStyleLumaInterpolationFilterPosj_SSSE3,             //j
            AvcStyleLumaInterpolationFilterPosk_SSSE3,             //k
            AvcStyleLumaInterpolationFilterVertical_SSSE3_INTRIN,         //n
            AvcStyleLumaInterpolationFilterPosp_SSSE3,             //p
            AvcStyleLumaInterpolationFilterPosq_SSSE3,             //q
            AvcStyleLumaInterpolationFilterPosr_SSSE3,             //r
        },
        // AVX2
        {
            avc_style_copy_sse2,                                    //A
            AvcStyleLumaInterpolationFilterHorizontal_SSSE3_INTRIN,       //a
            AvcStyleLumaInterpolationFilterHorizontal_SSSE3_INTRIN,       //b
            AvcStyleLumaInterpolationFilterHorizontal_SSSE3_INTRIN,       //c
            AvcStyleLumaInterpolationFilterVertical_SSSE3_INTRIN,         //d
            AvcStyleLumaInterpolationFilterPose_SSSE3,             //e
            AvcStyleLumaInterpolationFilterPosf_SSSE3,             //f
            AvcStyleLumaInterpolationFilterPosg_SSSE3,             //g
            AvcStyleLumaInterpolationFilterVertical_SSSE3_INTRIN,         //h
            AvcStyleLumaInterpolationFilterPosi_SSSE3,             //i
            AvcStyleLumaInterpolationFilterPosj_SSSE3,             //j
            AvcStyleLumaInterpolationFilterPosk_SSSE3,             //k
            AvcStyleLumaInterpolationFilterVertical_SSSE3_INTRIN,         //n
            AvcStyleLumaInterpolationFilterPosp_SSSE3,             //p
            AvcStyleLumaInterpolationFilterPosq_SSSE3,             //q
            AvcStyleLumaInterpolationFilterPosr_SSSE3,             //r
        },
    };

    static const PictureAverage FUNC_TABLE picture_average_array[ASM_TYPE_TOTAL] = {
        // NON_AVX2
        picture_average_kernel_sse2_intrin,
        // AVX2
        picture_average_kernel_sse2_intrin,
    };

    typedef void(*PictureAverage1Line)(
        EbByte                  src0,
        EbByte                  src1,
        EbByte                  dst,
        uint32_t                   area_width);

    static const PictureAverage1Line FUNC_TABLE picture_average1_line_array[ASM_TYPE_TOTAL] = {
        // NON_AVX2
        picture_average_kernel1_line_sse2_intrin,
        // AVX2
        picture_average_kernel1_line_sse2_intrin,
    };

#ifdef __cplusplus
}
#endif
#endif //EBAVCSTYLEMCP_H