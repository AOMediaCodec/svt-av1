/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

/*
* Copyright (c) 2016, Alliance for Open Media. All rights reserved
*
* This source code is subject to the terms of the BSD 2 Clause License and
* the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
* was not distributed with this source code in the LICENSE file, you can
* obtain it at www.aomedia.org/license/software. If the Alliance for Open
* Media Patent License 1.0 was not distributed with this source code in the
* PATENTS file, you can obtain it at www.aomedia.org/license/patent.
*/

// This file is generated. Do not edit.
#ifndef AOM_DSP_RTCD_H_
#define AOM_DSP_RTCD_H_

#include "EbDefinitions.h"
#include "EbPictureBufferDesc.h"

#ifdef RTCD_C
#define RTCD_EXTERN                //CHKN RTCD call in effect. declare the function pointers in  encHandle.
#else
#define RTCD_EXTERN extern         //CHKN run time externing the fucntion pointers.
#endif

/*
 * DSP
 */

#define HAS_MMX 0x01
#define HAS_SSE 0x02
#define HAS_SSE2 0x04
#define HAS_SSE3 0x08
#define HAS_SSSE3 0x10
#define HAS_SSE4_1 0x20
#define HAS_AVX 0x40
#define HAS_AVX2 0x80
#define HAS_SSE4_2 0x100

#ifdef __cplusplus
extern "C" {
#endif

    //to not include convolve.h, just forward declare what's needed.
    struct ConvolveParams;
    struct InterpFilterParams;

    void picture_average_kernel_sse2_intrin(EbByte   src0, uint32_t src0_stride, EbByte   src1, uint32_t src1_stride, EbByte   dst, uint32_t dst_stride, uint32_t area_width, uint32_t area_height);
    RTCD_EXTERN void(*picture_average)(EbByte   src0, uint32_t src0_stride, EbByte   src1, uint32_t src1_stride, EbByte   dst, uint32_t dst_stride, uint32_t area_width, uint32_t area_height);
    
    void picture_average_kernel1_line_sse2_intrin(EbByte src0,EbByte src1,EbByte dst,uint32_t area_width);
    RTCD_EXTERN void(*picture_average1_line)(EbByte src0,EbByte src1,EbByte dst,uint32_t area_width);

    void sad_loop_kernel_sse4_1_intrin(uint8_t  *src, uint32_t  src_stride, uint8_t  *ref, uint32_t  ref_stride, uint32_t  height, uint32_t  width, uint64_t *best_sad, int16_t *x_search_center, int16_t *y_search_center, uint32_t  src_stride_raw, int16_t search_area_width, int16_t search_area_height);
    void sad_loop_kernel_avx2_intrin(uint8_t  *src, uint32_t  src_stride, uint8_t  *ref, uint32_t  ref_stride, uint32_t  height, uint32_t  width, uint64_t *best_sad, int16_t *x_search_center, int16_t *y_search_center, uint32_t  src_stride_raw, int16_t search_area_width, int16_t search_area_height);
    RTCD_EXTERN void(*nxm_sad_loop_kernel)(uint8_t  *src, uint32_t  src_stride, uint8_t  *ref, uint32_t  ref_stride, uint32_t  height, uint32_t  width, uint64_t *best_sad, int16_t *x_search_center, int16_t *y_search_center, uint32_t  src_stride_raw, int16_t search_area_width, int16_t search_area_height);

    void sad_loop_kernel_sparse_sse4_1_intrin(uint8_t  *src, uint32_t  src_stride, uint8_t  *ref, uint32_t  ref_stride, uint32_t  height, uint32_t  width, uint64_t *best_sad, int16_t *x_search_center, int16_t *y_search_center, uint32_t  src_stride_raw, int16_t search_area_width, int16_t search_area_height);
    void sad_loop_kernel_sparse_avx2_intrin(uint8_t  *src, uint32_t  src_stride, uint8_t  *ref, uint32_t  ref_stride, uint32_t  height, uint32_t  width, uint64_t *best_sad, int16_t *x_search_center, int16_t *y_search_center, uint32_t  src_stride_raw, int16_t search_area_width, int16_t search_area_height);
    RTCD_EXTERN void(*nxm_sad_loop_kernel_sparse)(uint8_t  *src, uint32_t  src_stride, uint8_t  *ref, uint32_t  ref_stride, uint32_t  height, uint32_t  width, uint64_t *best_sad, int16_t *x_search_center, int16_t *y_search_center, uint32_t  src_stride_raw, int16_t search_area_width, int16_t search_area_height);

    void get_eight_horizontal_search_point_results_8x8_16x16_pu_sse41_intrin(uint8_t *src, uint32_t src_stride, uint8_t *ref, uint32_t ref_stride, uint32_t  *p_best_sad8x8, uint32_t *p_best_mv8x8, uint32_t *p_best_sad16x16, uint32_t *p_best_mv16x16, uint32_t mv, uint16_t *p_sad16x16, EbBool sub_sad);
    void get_eight_horizontal_search_point_results_8x8_16x16_pu_avx2_intrin(uint8_t *src, uint32_t src_stride, uint8_t *ref, uint32_t ref_stride, uint32_t  *p_best_sad8x8, uint32_t *p_best_mv8x8, uint32_t *p_best_sad16x16, uint32_t *p_best_mv16x16, uint32_t mv, uint16_t *p_sad16x16, EbBool sub_sad);
    RTCD_EXTERN void(*get_eight_horizontal_search_point_results_8x8_16x16)(uint8_t *src, uint32_t src_stride, uint8_t *ref, uint32_t ref_stride, uint32_t  *p_best_sad8x8, uint32_t *p_best_mv8x8, uint32_t *p_best_sad16x16, uint32_t *p_best_mv16x16, uint32_t mv, uint16_t *p_sad16x16, EbBool sub_sad);

    void get_eight_horizontal_search_point_results_32x32_64x64_pu_sse41_intrin(uint16_t  *p_sad16x16, uint32_t  *p_best_sad32x32, uint32_t  *p_best_sad64x64, uint32_t  *p_best_mv32x32, uint32_t  *p_best_mv64x64, uint32_t   mv);
    void get_eight_horizontal_search_point_results_32x32_64x64_pu_avx2_intrin(uint16_t  *p_sad16x16, uint32_t  *p_best_sad32x32, uint32_t  *p_best_sad64x64, uint32_t  *p_best_mv32x32, uint32_t  *p_best_mv64x64, uint32_t   mv);
    RTCD_EXTERN void(*get_eight_horizontal_search_point_results_32x32_64x64)(uint16_t  *p_sad16x16,uint32_t  *p_best_sad32x32,uint32_t  *p_best_sad64x64,uint32_t  *p_best_mv32x32,uint32_t  *p_best_mv64x64,uint32_t   mv);

    uint32_t combined_averaging_ssd_c(uint8_t *src, ptrdiff_t src_stride, uint8_t *ref1, ptrdiff_t ref1_stride, uint8_t *ref2, ptrdiff_t ref2_stride, uint32_t height, uint32_t width);
    uint32_t combined_averaging_ssd_avx2(uint8_t *src, ptrdiff_t src_stride, uint8_t *ref1, ptrdiff_t ref1_stride, uint8_t *ref2, ptrdiff_t ref2_stride, uint32_t height, uint32_t width);
    RTCD_EXTERN uint32_t(*combined_averaging_ssd)(uint8_t *src, ptrdiff_t src_stride,uint8_t *ref1, ptrdiff_t ref1_stride,uint8_t *ref2, ptrdiff_t ref2_stride,uint32_t height, uint32_t width);
    
    void ebav1_smooth_h_predictor_c(const uint32_t size, uint8_t *ref_samples, uint8_t *dst, const uint32_t prediction_buffer_stride, const EbBool skip);
    RTCD_EXTERN void(*ebav1_smooth_h_predictor)(const uint32_t size, uint8_t *ref_samples, uint8_t *dst, const uint32_t prediction_buffer_stride, const EbBool skip);

    void ebav1_smooth_v_predictor_c(const uint32_t   size, uint8_t  *ref_samples, uint8_t *dst, const uint32_t prediction_buffer_stride, const EbBool skip);
    RTCD_EXTERN void(*ebav1_smooth_v_predictor)(const uint32_t   size, uint8_t  *ref_samples, uint8_t *dst, const uint32_t prediction_buffer_stride, const EbBool skip);

    void noise_extract_luma_weak_avx2_intrin(EbPictureBufferDesc *input_picture_ptr, EbPictureBufferDesc *denoised_picture_ptr, EbPictureBufferDesc *noise_picture_ptr, uint32_t sb_origin_y, uint32_t sb_origin_x);
    void noise_extract_luma_weak_c(EbPictureBufferDesc *input_picture_ptr, EbPictureBufferDesc *denoised_picture_ptr, EbPictureBufferDesc *noise_picture_ptr, uint32_t sb_origin_y, uint32_t sb_origin_x);
    RTCD_EXTERN void(*noise_extract_luma_weak)(EbPictureBufferDesc *input_picture_ptr, EbPictureBufferDesc *denoised_picture_ptr, EbPictureBufferDesc *noise_picture_ptr, uint32_t sb_origin_y, uint32_t sb_origin_x);

    void noise_extract_luma_weak_lcu_c(EbPictureBufferDesc *input_picture_ptr, EbPictureBufferDesc *denoised_picture_ptr, EbPictureBufferDesc *noise_picture_ptr, uint32_t sb_origin_y, uint32_t sb_origin_x);
    void noise_extract_luma_weak_lcu_avx2_intrin(EbPictureBufferDesc *input_picture_ptr, EbPictureBufferDesc *denoised_picture_ptr, EbPictureBufferDesc *noise_picture_ptr, uint32_t sb_origin_y, uint32_t sb_origin_x);
    RTCD_EXTERN void(*noise_extract_luma_weak_lcu)(EbPictureBufferDesc *input_picture_ptr, EbPictureBufferDesc *denoised_picture_ptr, EbPictureBufferDesc *noise_picture_ptr, uint32_t sb_origin_y, uint32_t sb_origin_x);

    void noise_extract_luma_strong_c(EbPictureBufferDesc *input_picture_ptr, EbPictureBufferDesc *denoised_picture_ptr, uint32_t sb_origin_y, uint32_t sb_origin_x);
    void noise_extract_luma_strong_avx2_intrin(EbPictureBufferDesc *input_picture_ptr, EbPictureBufferDesc *denoised_picture_ptr, uint32_t sb_origin_y, uint32_t sb_origin_x);
    RTCD_EXTERN void(*strong_luma_filter)(EbPictureBufferDesc *input_picture_ptr, EbPictureBufferDesc *denoised_picture_ptr, uint32_t sb_origin_y,uint32_t sb_origin_x);

    void noise_extract_chroma_strong_c(EbPictureBufferDesc *input_picture_ptr, EbPictureBufferDesc *denoised_picture_ptr, uint32_t sb_origin_y, uint32_t sb_origin_x);
    void noise_extract_chroma_strong_avx2_intrin(EbPictureBufferDesc *input_picture_ptr, EbPictureBufferDesc *denoised_picture_ptr, uint32_t sb_origin_y, uint32_t sb_origin_x);
    RTCD_EXTERN void(*noise_extract_chroma_strong)(EbPictureBufferDesc *input_picture_ptr, EbPictureBufferDesc *denoised_picture_ptr, uint32_t sb_origin_y, uint32_t sb_origin_x);

    void noise_extract_chroma_weak_c(EbPictureBufferDesc *input_picture_ptr, EbPictureBufferDesc *denoised_picture_ptr, uint32_t sb_origin_y, uint32_t sb_origin_x);
    void noise_extract_chroma_weak_avx2_intrin(EbPictureBufferDesc *input_picture_ptr, EbPictureBufferDesc *denoised_picture_ptr, uint32_t sb_origin_y, uint32_t sb_origin_x);
    RTCD_EXTERN void(*noise_extract_chroma_weak)(EbPictureBufferDesc *input_picture_ptr, EbPictureBufferDesc *denoised_picture_ptr, uint32_t sb_origin_y, uint32_t sb_origin_x);

    int32_t sum_residual_c(int16_t * in_ptr, uint32_t size, uint32_t stride_in);
    int32_t sum_residual8bit_avx2_intrin(int16_t * in_ptr, uint32_t size, uint32_t stride_in);
    RTCD_EXTERN int32_t(*sum_residual)(int16_t * in_ptr, uint32_t size, uint32_t stride_in);

    void memset16bit_block_c(int16_t * in_ptr, uint32_t stride_in, uint32_t size, int16_t value);
    void memset16bit_block_avx2_intrin(int16_t * in_ptr, uint32_t stride_in, uint32_t size, int16_t value);
    RTCD_EXTERN void(*memset16bit_block)(int16_t * in_ptr, uint32_t stride_in, uint32_t size, int16_t value);

    void eb_apply_selfguided_restoration_c(const uint8_t *dat, int32_t width, int32_t height, int32_t stride, int32_t eps, const int32_t *xqd, uint8_t *dst, int32_t dst_stride, int32_t *tmpbuf, int32_t bit_depth, int32_t highbd);
    void eb_apply_selfguided_restoration_avx2(const uint8_t *dat, int32_t width, int32_t height, int32_t stride, int32_t eps, const int32_t *xqd, uint8_t *dst, int32_t dst_stride, int32_t *tmpbuf, int32_t bit_depth, int32_t highbd);
    RTCD_EXTERN void(*eb_apply_selfguided_restoration)(const uint8_t *dat, int32_t width, int32_t height, int32_t stride, int32_t eps, const int32_t *xqd, uint8_t *dst, int32_t dst_stride, int32_t *tmpbuf, int32_t bit_depth, int32_t highbd);

    void eb_av1_wiener_convolve_add_src_c(const uint8_t *src, ptrdiff_t src_stride, uint8_t *dst, ptrdiff_t dst_stride, const int16_t *filter_x, int32_t x_step_q4, const int16_t *filter_y, int32_t y_step_q4, int32_t w, int32_t h, const ConvolveParams *conv_params);
    void eb_av1_wiener_convolve_add_src_avx2(const uint8_t *src, ptrdiff_t src_stride, uint8_t *dst, ptrdiff_t dst_stride, const int16_t *filter_x, int32_t x_step_q4, const int16_t *filter_y, int32_t y_step_q4, int32_t w, int32_t h, const ConvolveParams *conv_params);
    RTCD_EXTERN void(*eb_av1_wiener_convolve_add_src)(const uint8_t *src, ptrdiff_t src_stride, uint8_t *dst, ptrdiff_t dst_stride, const int16_t *filter_x, int32_t x_step_q4, const int16_t *filter_y, int32_t y_step_q4, int32_t w, int32_t h, const ConvolveParams *conv_params);

    void eb_av1_highbd_wiener_convolve_add_src_c(const uint8_t *src, ptrdiff_t src_stride, uint8_t *dst, ptrdiff_t dst_stride, const int16_t *filter_x, int32_t x_step_q4, const int16_t *filter_y, int32_t y_step_q4, int32_t w, int32_t h, const ConvolveParams *conv_params, int32_t bps);
    void eb_av1_highbd_wiener_convolve_add_src_avx2(const uint8_t *src, ptrdiff_t src_stride, uint8_t *dst, ptrdiff_t dst_stride, const int16_t *filter_x, int32_t x_step_q4, const int16_t *filter_y, int32_t y_step_q4, int32_t w, int32_t h, const ConvolveParams *conv_params, int32_t bps);
    RTCD_EXTERN void(*eb_av1_highbd_wiener_convolve_add_src)(const uint8_t *src, ptrdiff_t src_stride, uint8_t *dst, ptrdiff_t dst_stride, const int16_t *filter_x, int32_t x_step_q4, const int16_t *filter_y, int32_t y_step_q4, int32_t w, int32_t h, const ConvolveParams *conv_params, int32_t bps);

    void apply_temp_filtering_32x32_c(const uint8_t *y_frame1, int y_stride, const uint8_t *y_pred, int y_buf_stride, const uint8_t *u_frame1, const uint8_t *v_frame1, int uv_stride, const uint8_t *u_pred, const uint8_t *v_pred, int uv_buf_stride, unsigned int block_width, unsigned int block_height, int ss_x, int ss_y, int strength, const int *blk_fw, int use_32x32, uint32_t *y_accumulator, uint16_t *y_count, uint32_t *u_accumulator, uint16_t *u_count, uint32_t *v_accumulator, uint16_t *v_count);
    void apply_temp_filtering_32x32_sse4_1(const uint8_t *y_frame1, int y_stride, const uint8_t *y_pred, int y_buf_stride, const uint8_t *u_frame1, const uint8_t *v_frame1, int uv_stride, const uint8_t *u_pred, const uint8_t *v_pred, int uv_buf_stride, unsigned int block_width, unsigned int block_height, int ss_x, int ss_y, int strength, const int *blk_fw, int use_32x32, uint32_t *y_accumulator, uint16_t *y_count, uint32_t *u_accumulator, uint16_t *u_count, uint32_t *v_accumulator, uint16_t *v_count);
    RTCD_EXTERN void(*apply_temp_filtering_32x32)(const uint8_t *y_frame1, int y_stride, const uint8_t *y_pred, int y_buf_stride, const uint8_t *u_frame1, const uint8_t *v_frame1, int uv_stride, const uint8_t *u_pred, const uint8_t *v_pred, int uv_buf_stride, unsigned int block_width, unsigned int block_height, int ss_x, int ss_y, int strength, const int *blk_fw, int use_32x32, uint32_t *y_accumulator, uint16_t *y_count, uint32_t *u_accumulator, uint16_t *u_count, uint32_t *v_accumulator, uint16_t *v_count);

    uint32_t compute4x4SAD_Kernel_c(const uint8_t  *src, uint32_t  src_stride, const uint8_t  *ref, uint32_t  ref_stride, uint32_t  height, uint32_t  width);
    uint32_t compute4x_m_sad_avx2_intrin(const uint8_t  *src, uint32_t  src_stride, const uint8_t  *ref, uint32_t  ref_stride, uint32_t  height, uint32_t  width);
    RTCD_EXTERN uint32_t(*compute4x4_SAD)(const uint8_t  *src, uint32_t  src_stride, const uint8_t  *ref, uint32_t  ref_stride, uint32_t  height, uint32_t  width);

    void full_distortion_kernel32_bits_c(int32_t *coeff, uint32_t coeff_stride, int32_t *recon_coeff, uint32_t recon_coeff_stride, uint64_t distortion_result[DIST_CALC_TOTAL], uint32_t area_width, uint32_t area_height);
    void full_distortion_kernel32_bits_avx2(int32_t *coeff, uint32_t coeff_stride, int32_t *recon_coeff, uint32_t recon_coeff_stride, uint64_t distortion_result[DIST_CALC_TOTAL], uint32_t area_width, uint32_t area_height);
    RTCD_EXTERN void(*full_distortion_kernel32_bits)(int32_t *coeff, uint32_t coeff_stride, int32_t *recon_coeff, uint32_t recon_coeff_stride, uint64_t distortion_result[DIST_CALC_TOTAL], uint32_t area_width, uint32_t area_height);

    void full_distortion_kernel_cbf_zero32_bits_c(int32_t *coeff,uint32_t coeff_stride,int32_t *recon_coeff,uint32_t recon_coeff_stride,uint64_t distortion_result[DIST_CALC_TOTAL],uint32_t area_width,uint32_t area_height);
    void full_distortion_kernel_cbf_zero32_bits_avx2(int32_t *coeff, uint32_t coeff_stride, int32_t *recon_coeff, uint32_t recon_coeff_stride, uint64_t distortion_result[DIST_CALC_TOTAL], uint32_t area_width, uint32_t area_height);
    RTCD_EXTERN void (*full_distortion_kernel_cbf_zero32_bits)(int32_t *coeff, uint32_t coeff_stride, int32_t *recon_coeff, uint32_t recon_coeff_stride, uint64_t distortion_result[DIST_CALC_TOTAL], uint32_t area_width, uint32_t area_height);
      
    uint64_t spatial_full_distortion_kernel_c(uint8_t *input,uint32_t input_offset,uint32_t input_stride,uint8_t *recon,uint32_t recon_offset,uint32_t recon_stride,uint32_t area_width,uint32_t area_height);
#ifndef NON_AVX512_SUPPORT
    uint64_t spatial_full_distortion_kernel_avx512(uint8_t *input, uint32_t input_offset, uint32_t input_stride, uint8_t *recon, uint32_t recon_offset, uint32_t recon_stride, uint32_t area_width, uint32_t area_height);
#endif
    uint64_t spatial_full_distortion_kernel_avx2(uint8_t *input, uint32_t input_offset, uint32_t input_stride, uint8_t *recon, uint32_t recon_offset, uint32_t recon_stride, uint32_t area_width, uint32_t area_height);
    RTCD_EXTERN uint64_t(*spatial_full_distortion_kernel)(uint8_t *input, uint32_t input_offset, uint32_t input_stride, uint8_t *recon, uint32_t recon_offset, uint32_t recon_stride, uint32_t area_width, uint32_t area_height);

    void eb_av1_selfguided_restoration_c(const uint8_t *dgd8, int32_t width, int32_t height, int32_t dgd_stride, int32_t *flt0, int32_t *flt1, int32_t flt_stride, int32_t sgr_params_idx, int32_t bit_depth, int32_t highbd);
    void eb_av1_selfguided_restoration_avx2(const uint8_t *dgd8, int32_t width, int32_t height, int32_t dgd_stride, int32_t *flt0, int32_t *flt1, int32_t flt_stride, int32_t sgr_params_idx, int32_t bit_depth, int32_t highbd);
    RTCD_EXTERN void(*eb_av1_selfguided_restoration)(const uint8_t *dgd8, int32_t width, int32_t height, int32_t dgd_stride, int32_t *flt0, int32_t *flt1, int32_t flt_stride, int32_t sgr_params_idx, int32_t bit_depth, int32_t highbd);

    uint32_t compute8x4_sad_kernel_c(uint8_t *src, uint32_t src_stride, uint8_t *ref, uint32_t ref_stride);
    RTCD_EXTERN uint32_t (*compute8x4_sad_kernel)(uint8_t *src, uint32_t src_stride, uint8_t *ref, uint32_t ref_stride);

    uint64_t compute8x8_satd_u8_sse4(uint8_t  *src, uint64_t *dc_value, uint32_t  src_stride);
    RTCD_EXTERN uint64_t (*compute8x8_satd_u8)(uint8_t  *src, uint64_t *dc_value, uint32_t  src_stride);
    
    void picture_addition_kernel16bit_sse2_intrin(uint16_t  *pred_ptr,uint32_t  pred_stride,int16_t *residual_ptr,uint32_t  residual_stride,uint16_t  *recon_ptr,uint32_t  recon_stride,uint32_t  width,uint32_t  height);
    RTCD_EXTERN void (*picture_addition_kernel16bit)(uint16_t  *pred_ptr, uint32_t  pred_stride, int16_t *residual_ptr, uint32_t  residual_stride, uint16_t  *recon_ptr, uint32_t  recon_stride, uint32_t  width, uint32_t  height);
    
    void avc_style_luma_interpolation_filter_ssse3_helper(EbByte ref_pic,uint32_t src_stride,EbByte dst,uint32_t dst_stride,uint32_t pu_width,uint32_t pu_height,EbByte temp_buf,EbBool skip,uint32_t frac_pos,uint8_t choice);
    RTCD_EXTERN void (*avc_style_luma_interpolation_filter)(EbByte ref_pic, uint32_t src_stride, EbByte dst, uint32_t dst_stride, uint32_t pu_width, uint32_t pu_height, EbByte temp_buf, EbBool skip, uint32_t frac_pos, uint8_t choice);

    int32_t eb_cdef_find_dir_c(const uint16_t *img, int32_t stride, int32_t *var, int32_t coeff_shift);
    int32_t eb_cdef_find_dir_avx2(const uint16_t *img, int32_t stride, int32_t *var, int32_t coeff_shift);
    RTCD_EXTERN int32_t(*eb_cdef_find_dir)(const uint16_t *img, int32_t stride, int32_t *var, int32_t coeff_shift);

    void eb_cdef_filter_block_c(uint8_t *dst8, uint16_t *dst16, int32_t dstride, const uint16_t *in, int32_t pri_strength, int32_t sec_strength, int32_t dir, int32_t pri_damping, int32_t sec_damping, int32_t bsize, int32_t max, int32_t coeff_shift);
    void eb_cdef_filter_block_avx2(uint8_t *dst8, uint16_t *dst16, int32_t dstride, const uint16_t *in, int32_t pri_strength, int32_t sec_strength, int32_t dir, int32_t pri_damping, int32_t sec_damping, int32_t bsize, int32_t max, int32_t coeff_shift);
    RTCD_EXTERN void(*eb_cdef_filter_block)(uint8_t *dst8, uint16_t *dst16, int32_t dstride, const uint16_t *in, int32_t pri_strength, int32_t sec_strength, int32_t dir, int32_t pri_damping, int32_t sec_damping, int32_t bsize, int32_t max, int32_t coeff_shift);

    uint64_t compute_cdef_dist_c(const uint16_t *dst, int32_t dstride, const uint16_t *src, const cdef_list *dlist, int32_t cdef_count, BlockSize bsize, int32_t coeff_shift, int32_t pli);
    uint64_t compute_cdef_dist_avx2(const uint16_t *dst, int32_t dstride, const uint16_t *src, const cdef_list *dlist, int32_t cdef_count, BlockSize bsize, int32_t coeff_shift, int32_t pli);
    RTCD_EXTERN uint64_t(*eb_compute_cdef_dist)(const uint16_t *dst, int32_t dstride, const uint16_t *src, const cdef_list *dlist, int32_t cdef_count, BlockSize bsize, int32_t coeff_shift, int32_t pli);
    void eb_copy_rect8_8bit_to_16bit_c(uint16_t *dst, int32_t dstride, const uint8_t *src, int32_t sstride, int32_t v, int32_t h);
    void eb_copy_rect8_8bit_to_16bit_avx2(uint16_t *dst, int32_t dstride, const uint8_t *src, int32_t sstride, int32_t v, int32_t h);
    RTCD_EXTERN void(*eb_copy_rect8_8bit_to_16bit)(uint16_t *dst, int32_t dstride, const uint8_t *src, int32_t sstride, int32_t v, int32_t h);

    void eb_av1_compute_stats_c(int32_t wiener_win, const uint8_t *dgd8, const uint8_t *src8, int32_t h_start, int32_t h_end, int32_t v_start, int32_t v_end, int32_t dgd_stride, int32_t src_stride, int64_t *M, int64_t *H);
    void eb_av1_compute_stats_avx2(int32_t wiener_win, const uint8_t *dgd8, const uint8_t *src8, int32_t h_start, int32_t h_end, int32_t v_start, int32_t v_end, int32_t dgd_stride, int32_t src_stride, int64_t *M, int64_t *H);
    RTCD_EXTERN void(*eb_av1_compute_stats)(int32_t wiener_win, const uint8_t *dgd8, const uint8_t *src8, int32_t h_start, int32_t h_end, int32_t v_start, int32_t v_end, int32_t dgd_stride, int32_t src_stride, int64_t *M, int64_t *H);

    void eb_av1_compute_stats_highbd_c(int32_t wiener_win, const uint8_t *dgd8, const uint8_t *src8, int32_t h_start, int32_t h_end, int32_t v_start, int32_t v_end, int32_t dgd_stride, int32_t src_stride, int64_t *M, int64_t *H, AomBitDepth bit_depth);
    void eb_av1_compute_stats_highbd_avx2(int32_t wiener_win, const uint8_t *dgd8, const uint8_t *src8, int32_t h_start, int32_t h_end, int32_t v_start, int32_t v_end, int32_t dgd_stride, int32_t src_stride, int64_t *M, int64_t *H, AomBitDepth bit_depth);
    RTCD_EXTERN void(*eb_av1_compute_stats_highbd)(int32_t wiener_win, const uint8_t *dgd8, const uint8_t *src8, int32_t h_start, int32_t h_end, int32_t v_start, int32_t v_end, int32_t dgd_stride, int32_t src_stride, int64_t *M, int64_t *H, AomBitDepth bit_depth);

    int64_t eb_av1_lowbd_pixel_proj_error_c(const uint8_t *src8, int32_t width, int32_t height, int32_t src_stride, const uint8_t *dat8, int32_t dat_stride, int32_t *flt0, int32_t flt0_stride, int32_t *flt1, int32_t flt1_stride, int32_t xq[2], const SgrParamsType *params);
    int64_t eb_av1_lowbd_pixel_proj_error_avx2(const uint8_t *src8, int32_t width, int32_t height, int32_t src_stride, const uint8_t *dat8, int32_t dat_stride, int32_t *flt0, int32_t flt0_stride, int32_t *flt1, int32_t flt1_stride, int32_t xq[2], const SgrParamsType *params);
    RTCD_EXTERN int64_t(*eb_av1_lowbd_pixel_proj_error)(const uint8_t *src8, int32_t width, int32_t height, int32_t src_stride, const uint8_t *dat8, int32_t dat_stride, int32_t *flt0, int32_t flt0_stride, int32_t *flt1, int32_t flt1_stride, int32_t xq[2], const SgrParamsType *params);

    int64_t eb_av1_highbd_pixel_proj_error_c(const uint8_t *src8, int32_t width, int32_t height, int32_t src_stride, const uint8_t *dat8, int32_t dat_stride, int32_t *flt0, int32_t flt0_stride, int32_t *flt1, int32_t flt1_stride, int32_t xq[2], const SgrParamsType *params);
    int64_t eb_av1_highbd_pixel_proj_error_avx2(const uint8_t *src8, int32_t width, int32_t height, int32_t src_stride, const uint8_t *dat8, int32_t dat_stride, int32_t *flt0, int32_t flt0_stride, int32_t *flt1, int32_t flt1_stride, int32_t xq[2], const SgrParamsType *params);
    RTCD_EXTERN int64_t(*eb_av1_highbd_pixel_proj_error)(const uint8_t *src8, int32_t width, int32_t height, int32_t src_stride, const uint8_t *dat8, int32_t dat_stride, int32_t *flt0, int32_t flt0_stride, int32_t *flt1, int32_t flt1_stride, int32_t xq[2], const SgrParamsType *params);

    void eb_av1_highbd_convolve_2d_copy_sr_c(const uint16_t *src, int32_t src_stride, uint16_t *dst, int32_t dst_stride, int32_t w, int32_t h, const InterpFilterParams *filter_params_x, const InterpFilterParams *filter_params_y, const int32_t subpel_x_q4, const int32_t subpel_y_q4, ConvolveParams *conv_params, int32_t bd);
    void eb_av1_highbd_convolve_2d_copy_sr_avx2(const uint16_t *src, int32_t src_stride, uint16_t *dst, int32_t dst_stride, int32_t w, int32_t h, const InterpFilterParams *filter_params_x, const InterpFilterParams *filter_params_y, const int32_t subpel_x_q4, const int32_t subpel_y_q4, ConvolveParams *conv_params, int32_t bd);
    RTCD_EXTERN void(*eb_av1_highbd_convolve_2d_copy_sr)(const uint16_t *src, int32_t src_stride, uint16_t *dst, int32_t dst_stride, int32_t w, int32_t h, const InterpFilterParams *filter_params_x, const InterpFilterParams *filter_params_y, const int32_t subpel_x_q4, const int32_t subpel_y_q4, ConvolveParams *conv_params, int32_t bd);

    void eb_av1_highbd_jnt_convolve_2d_copy_c(const uint16_t *src, int32_t src_stride, uint16_t *dst, int32_t dst_stride, int32_t w, int32_t h, const InterpFilterParams *filter_params_x, const InterpFilterParams *filter_params_y, const int32_t subpel_x_q4, const int32_t subpel_y_q4, ConvolveParams *conv_params, int32_t bd);
    void eb_av1_highbd_jnt_convolve_2d_copy_avx2(const uint16_t *src, int32_t src_stride, uint16_t *dst, int32_t dst_stride, int32_t w, int32_t h, const InterpFilterParams *filter_params_x, const InterpFilterParams *filter_params_y, const int32_t subpel_x_q4, const int32_t subpel_y_q4, ConvolveParams *conv_params, int32_t bd);
    RTCD_EXTERN void(*eb_av1_highbd_jnt_convolve_2d_copy)(const uint16_t *src, int32_t src_stride, uint16_t *dst, int32_t dst_stride, int32_t w, int32_t h, const InterpFilterParams *filter_params_x, const InterpFilterParams *filter_params_y, const int32_t subpel_x_q4, const int32_t subpel_y_q4, ConvolveParams *conv_params, int32_t bd);

    void eb_av1_highbd_convolve_y_sr_c(const uint16_t *src, int32_t src_stride, uint16_t *dst, int32_t dst_stride, int32_t w, int32_t h, const InterpFilterParams *filter_params_x, const InterpFilterParams *filter_params_y, const int32_t subpel_x_q4, const int32_t subpel_y_q4, ConvolveParams *conv_params, int32_t bd);
    void eb_av1_highbd_convolve_y_sr_avx2(const uint16_t *src, int32_t src_stride, uint16_t *dst, int32_t dst_stride, int32_t w, int32_t h, const InterpFilterParams *filter_params_x, const InterpFilterParams *filter_params_y, const int32_t subpel_x_q4, const int32_t subpel_y_q4, ConvolveParams *conv_params, int32_t bd);
    RTCD_EXTERN void(*eb_av1_highbd_convolve_y_sr)(const uint16_t *src, int32_t src_stride, uint16_t *dst, int32_t dst_stride, int32_t w, int32_t h, const InterpFilterParams *filter_params_x, const InterpFilterParams *filter_params_y, const int32_t subpel_x_q4, const int32_t subpel_y_q4, ConvolveParams *conv_params, int32_t bd);

    void eb_av1_highbd_convolve_2d_sr_c(const uint16_t *src, int32_t src_stride, uint16_t *dst, int32_t dst_stride, int32_t w, int32_t h, const InterpFilterParams *filter_params_x, const InterpFilterParams *filter_params_y, const int32_t subpel_x_q4, const int32_t subpel_y_q4, ConvolveParams *conv_params, int32_t bd);
    void eb_av1_highbd_convolve_2d_sr_avx2(const uint16_t *src, int32_t src_stride, uint16_t *dst, int32_t dst_stride, int32_t w, int32_t h, const InterpFilterParams *filter_params_x, const InterpFilterParams *filter_params_y, const int32_t subpel_x_q4, const int32_t subpel_y_q4, ConvolveParams *conv_params, int32_t bd);
    RTCD_EXTERN void(*eb_av1_highbd_convolve_2d_sr)(const uint16_t *src, int32_t src_stride, uint16_t *dst, int32_t dst_stride, int32_t w, int32_t h, const InterpFilterParams *filter_params_x, const InterpFilterParams *filter_params_y, const int32_t subpel_x_q4, const int32_t subpel_y_q4, ConvolveParams *conv_params, int32_t bd);

    void eb_av1_highbd_jnt_convolve_2d_c(const uint16_t *src, int32_t src_stride, uint16_t *dst, int32_t dst_stride, int32_t w, int32_t h, const InterpFilterParams *filter_params_x, const InterpFilterParams *filter_params_y, const int32_t subpel_x_q4, const int32_t subpel_y_q4, ConvolveParams *conv_params, int32_t bd);
    void eb_av1_highbd_jnt_convolve_2d_avx2(const uint16_t *src, int32_t src_stride, uint16_t *dst, int32_t dst_stride, int32_t w, int32_t h, const InterpFilterParams *filter_params_x, const InterpFilterParams *filter_params_y, const int32_t subpel_x_q4, const int32_t subpel_y_q4, ConvolveParams *conv_params, int32_t bd);
    RTCD_EXTERN void(*eb_av1_highbd_jnt_convolve_2d)(const uint16_t *src, int32_t src_stride, uint16_t *dst, int32_t dst_stride, int32_t w, int32_t h, const InterpFilterParams *filter_params_x, const InterpFilterParams *filter_params_y, const int32_t subpel_x_q4, const int32_t subpel_y_q4, ConvolveParams *conv_params, int32_t bd);

    void eb_av1_highbd_jnt_convolve_x_c(const uint16_t *src, int32_t src_stride, uint16_t *dst, int32_t dst_stride, int32_t w, int32_t h, const InterpFilterParams *filter_params_x, const InterpFilterParams *filter_params_y, const int32_t subpel_x_q4, const int32_t subpel_y_q4, ConvolveParams *conv_params, int32_t bd);
    void eb_av1_highbd_jnt_convolve_x_avx2(const uint16_t *src, int32_t src_stride, uint16_t *dst, int32_t dst_stride, int32_t w, int32_t h, const InterpFilterParams *filter_params_x, const InterpFilterParams *filter_params_y, const int32_t subpel_x_q4, const int32_t subpel_y_q4, ConvolveParams *conv_params, int32_t bd);
    RTCD_EXTERN void(*eb_av1_highbd_jnt_convolve_x)(const uint16_t *src, int32_t src_stride, uint16_t *dst, int32_t dst_stride, int32_t w, int32_t h, const InterpFilterParams *filter_params_x, const InterpFilterParams *filter_params_y, const int32_t subpel_x_q4, const int32_t subpel_y_q4, ConvolveParams *conv_params, int32_t bd);

    void eb_av1_highbd_jnt_convolve_y_c(const uint16_t *src, int32_t src_stride, uint16_t *dst, int32_t dst_stride, int32_t w, int32_t h, const InterpFilterParams *filter_params_x, const InterpFilterParams *filter_params_y, const int32_t subpel_x_q4, const int32_t subpel_y_q4, ConvolveParams *conv_params, int32_t bd);
    void eb_av1_highbd_jnt_convolve_y_avx2(const uint16_t *src, int32_t src_stride, uint16_t *dst, int32_t dst_stride, int32_t w, int32_t h, const InterpFilterParams *filter_params_x, const InterpFilterParams *filter_params_y, const int32_t subpel_x_q4, const int32_t subpel_y_q4, ConvolveParams *conv_params, int32_t bd);
    RTCD_EXTERN void(*eb_av1_highbd_jnt_convolve_y)(const uint16_t *src, int32_t src_stride, uint16_t *dst, int32_t dst_stride, int32_t w, int32_t h, const InterpFilterParams *filter_params_x, const InterpFilterParams *filter_params_y, const int32_t subpel_x_q4, const int32_t subpel_y_q4, ConvolveParams *conv_params, int32_t bd);

    void eb_av1_highbd_convolve_x_sr_c(const uint16_t *src, int32_t src_stride, uint16_t *dst, int32_t dst_stride, int32_t w, int32_t h, const InterpFilterParams *filter_params_x, const InterpFilterParams *filter_params_y, const int32_t subpel_x_q4, const int32_t subpel_y_q4, ConvolveParams *conv_params, int32_t bd);
    void eb_av1_highbd_convolve_x_sr_avx2(const uint16_t *src, int32_t src_stride, uint16_t *dst, int32_t dst_stride, int32_t w, int32_t h, const InterpFilterParams *filter_params_x, const InterpFilterParams *filter_params_y, const int32_t subpel_x_q4, const int32_t subpel_y_q4, ConvolveParams *conv_params, int32_t bd);
    RTCD_EXTERN void(*eb_av1_highbd_convolve_x_sr)(const uint16_t *src, int32_t src_stride, uint16_t *dst, int32_t dst_stride, int32_t w, int32_t h, const InterpFilterParams *filter_params_x, const InterpFilterParams *filter_params_y, const int32_t subpel_x_q4, const int32_t subpel_y_q4, ConvolveParams *conv_params, int32_t bd);

    void eb_subtract_average_c(int16_t *pred_buf_q3, int32_t width, int32_t height, int32_t round_offset, int32_t num_pel_log2);
    void eb_subtract_average_avx2(int16_t *pred_buf_q3, int32_t width, int32_t height, int32_t round_offset, int32_t num_pel_log2);
    RTCD_EXTERN void(*eb_subtract_average)(int16_t *pred_buf_q3, int32_t width, int32_t height, int32_t round_offset, int32_t num_pel_log2);

    void eb_cfl_predict_lbd_c(const int16_t *pred_buf_q3, uint8_t *pred, int32_t pred_stride, uint8_t *dst, int32_t dst_stride, int32_t alpha_q3, int32_t bit_depth, int32_t width, int32_t height);
    void eb_cfl_predict_lbd_avx2(const int16_t *pred_buf_q3, uint8_t *pred, int32_t pred_stride, uint8_t *dst, int32_t dst_stride, int32_t alpha_q3, int32_t bit_depth, int32_t width, int32_t height);
    RTCD_EXTERN void(*eb_cfl_predict_lbd)(const int16_t *pred_buf_q3, uint8_t *pred, int32_t pred_stride, uint8_t *dst, int32_t dst_stride, int32_t alpha_q3, int32_t bit_depth, int32_t width, int32_t height);

    void eb_cfl_predict_hbd_c(const int16_t *pred_buf_q3, uint16_t *pred, int32_t pred_stride, uint16_t *dst, int32_t dst_stride, int32_t alpha_q3, int32_t bit_depth, int32_t width, int32_t height);
    void eb_cfl_predict_hbd_avx2(const int16_t *pred_buf_q3, uint16_t *pred, int32_t pred_stride, uint16_t *dst, int32_t dst_stride, int32_t alpha_q3, int32_t bit_depth, int32_t width, int32_t height);
    RTCD_EXTERN void(*eb_cfl_predict_hbd)(const int16_t *pred_buf_q3, uint16_t *pred, int32_t pred_stride, uint16_t *dst, int32_t dst_stride, int32_t alpha_q3, int32_t bit_depth, int32_t width, int32_t height);

    void eb_av1_filter_intra_edge_high_c_old(uint8_t *p, int32_t sz, int32_t strength);
    void eb_av1_filter_intra_edge_sse4_1(uint8_t *p, int32_t sz, int32_t strength);
    RTCD_EXTERN void(*eb_av1_filter_intra_edge)(uint8_t *p, int32_t sz, int32_t strength);

    void eb_av1_filter_intra_edge_high_c(uint16_t *p, int32_t sz, int32_t strength);
    void eb_av1_filter_intra_edge_high_sse4_1(uint16_t *p, int32_t sz, int32_t strength);
    RTCD_EXTERN void(*eb_av1_filter_intra_edge_high)(uint16_t *p, int32_t sz, int32_t strength);

    int64_t eb_av1_calc_frame_error_c(const uint8_t *const ref, int stride, const uint8_t *const dst, int p_width, int p_height, int p_stride);
    int64_t eb_av1_calc_frame_error_avx2(const uint8_t *const ref, int stride, const uint8_t *const dst, int p_width, int p_height, int p_stride);
    RTCD_EXTERN int64_t (*eb_av1_calc_frame_error)(const uint8_t *const ref, int stride, const uint8_t *const dst, int p_width, int p_height, int p_stride);

    void eb_av1_fwd_txfm2d_4x16_c(int16_t *input, int32_t *output, uint32_t input_stride, TxType transform_type, uint8_t  bit_depth);
    void eb_av1_fwd_txfm2d_4x16_avx2(int16_t *input, int32_t *output, uint32_t inputStride, TxType transform_type, uint8_t  bit_depth);
    RTCD_EXTERN void(*eb_av1_fwd_txfm2d_4x16)(int16_t *input, int32_t *output, uint32_t inputStride, TxType transform_type, uint8_t  bit_depth);

    void eb_av1_fwd_txfm2d_16x4_c(int16_t *input, int32_t *output, uint32_t input_stride, TxType transform_type, uint8_t  bit_depth);
    void eb_av1_fwd_txfm2d_16x4_avx2(int16_t *input, int32_t *output, uint32_t inputStride, TxType transform_type, uint8_t  bit_depth);
    RTCD_EXTERN void(*eb_av1_fwd_txfm2d_16x4)(int16_t *input, int32_t *output, uint32_t inputStride, TxType transform_type, uint8_t  bit_depth);

    void eb_av1_fwd_txfm2d_4x8_c(int16_t *input, int32_t *output, uint32_t input_stride, TxType transform_type, uint8_t  bit_depth);
    void eb_av1_fwd_txfm2d_4x8_avx2(int16_t *input, int32_t *output, uint32_t inputStride, TxType transform_type, uint8_t  bit_depth);
    RTCD_EXTERN void(*eb_av1_fwd_txfm2d_4x8)(int16_t *input, int32_t *output, uint32_t inputStride, TxType transform_type, uint8_t  bit_depth);

    void eb_av1_fwd_txfm2d_8x4_c(int16_t *input, int32_t *output, uint32_t input_stride, TxType transform_type, uint8_t  bit_depth);
    void eb_av1_fwd_txfm2d_8x4_avx2(int16_t *input, int32_t *output, uint32_t inputStride, TxType transform_type, uint8_t  bit_depth);
    RTCD_EXTERN void(*eb_av1_fwd_txfm2d_8x4)(int16_t *input, int32_t *output, uint32_t inputStride, TxType transform_type, uint8_t  bit_depth);

    void eb_av1_fwd_txfm2d_8x16_c(int16_t *input, int32_t *output, uint32_t input_stride, TxType transform_type, uint8_t  bit_depth);
    void eb_av1_fwd_txfm2d_8x16_avx2(int16_t *input, int32_t *output, uint32_t input_stride, TxType transform_type, uint8_t  bit_depth);
    RTCD_EXTERN void(*eb_av1_fwd_txfm2d_8x16)(int16_t *input, int32_t *output, uint32_t input_stride, TxType transform_type, uint8_t  bit_depth);

    void eb_av1_fwd_txfm2d_16x8_c(int16_t *input, int32_t *output, uint32_t input_stride, TxType transform_type, uint8_t  bit_depth);
    void eb_av1_fwd_txfm2d_16x8_avx2(int16_t *input, int32_t *output, uint32_t input_stride, TxType transform_type, uint8_t  bit_depth);
    RTCD_EXTERN void(*eb_av1_fwd_txfm2d_16x8)(int16_t *input, int32_t *output, uint32_t input_stride, TxType transform_type, uint8_t  bit_depth);

    void eb_av1_fwd_txfm2d_4x16_c(int16_t *input, int32_t *output, uint32_t input_stride, TxType transform_type, uint8_t  bit_depth);
    void eb_av1_fwd_txfm2d_4x16_sse4_1(int16_t *input, int32_t *output, uint32_t input_stride, TxType transform_type, uint8_t  bit_depth);
    RTCD_EXTERN void(*eb_av1_fwd_txfm2d_4x16)(int16_t *input, int32_t *output, uint32_t input_stride, TxType transform_type, uint8_t  bit_depth);

    void eb_av1_fwd_txfm2d_16x4_c(int16_t *input, int32_t *output, uint32_t input_stride, TxType transform_type, uint8_t  bit_depth);
    void eb_av1_fwd_txfm2d_16x4_sse4_1(int16_t *input, int32_t *output, uint32_t input_stride, TxType transform_type, uint8_t  bit_depth);
    RTCD_EXTERN void(*eb_av1_fwd_txfm2d_16x4)(int16_t *input, int32_t *output, uint32_t input_stride, TxType transform_type, uint8_t  bit_depth);

    void eb_av1_fwd_txfm2d_4x8_c(int16_t *input, int32_t *output, uint32_t input_stride, TxType transform_type, uint8_t  bit_depth);
    void eb_av1_fwd_txfm2d_4x8_sse4_1(int16_t *input, int32_t *output, uint32_t input_stride, TxType transform_type, uint8_t  bit_depth);
    RTCD_EXTERN void(*eb_av1_fwd_txfm2d_4x8)(int16_t *input, int32_t *output, uint32_t input_stride, TxType transform_type, uint8_t  bit_depth);

    void eb_av1_fwd_txfm2d_8x4_c(int16_t *input, int32_t *output, uint32_t input_stride, TxType transform_type, uint8_t  bit_depth);
    void eb_av1_fwd_txfm2d_8x4_sse4_1(int16_t *input, int32_t *output, uint32_t input_stride, TxType transform_type, uint8_t  bit_depth);
    RTCD_EXTERN void(*eb_av1_fwd_txfm2d_8x4)(int16_t *input, int32_t *output, uint32_t input_stride, TxType transform_type, uint8_t  bit_depth);

    void eb_av1_fwd_txfm2d_32x16_c(int16_t *input, int32_t *output, uint32_t input_stride, TxType transform_type, uint8_t  bit_depth);
    void eb_av1_fwd_txfm2d_32x16_avx2(int16_t *input, int32_t *output, uint32_t input_stride, TxType transform_type, uint8_t  bit_depth);
    RTCD_EXTERN void(*eb_av1_fwd_txfm2d_32x16)(int16_t *input, int32_t *output, uint32_t input_stride, TxType transform_type, uint8_t  bit_depth);

    void eb_av1_fwd_txfm2d_32x8_c(int16_t *input, int32_t *output, uint32_t input_stride, TxType transform_type, uint8_t  bit_depth);
    void eb_av1_fwd_txfm2d_32x8_avx2(int16_t *input, int32_t *output, uint32_t input_stride, TxType transform_type, uint8_t  bit_depth);
    RTCD_EXTERN void(*eb_av1_fwd_txfm2d_32x8)(int16_t *input, int32_t *output, uint32_t input_stride, TxType transform_type, uint8_t  bit_depth);

    void eb_av1_fwd_txfm2d_8x32_c(int16_t *input, int32_t *output, uint32_t input_stride, TxType transform_type, uint8_t  bit_depth);
    void eb_av1_fwd_txfm2d_8x32_avx2(int16_t *input, int32_t *output, uint32_t input_stride, TxType transform_type, uint8_t  bit_depth);
    RTCD_EXTERN void(*eb_av1_fwd_txfm2d_8x32)(int16_t *input, int32_t *output, uint32_t input_stride, TxType transform_type, uint8_t  bit_depth);

    void eb_av1_fwd_txfm2d_16x32_c(int16_t *input, int32_t *output, uint32_t input_stride, TxType transform_type, uint8_t  bit_depth);
    void eb_av1_fwd_txfm2d_16x32_avx2(int16_t *input, int32_t *output, uint32_t input_stride, TxType transform_type, uint8_t  bit_depth);
    RTCD_EXTERN void(*eb_av1_fwd_txfm2d_16x32)(int16_t *input, int32_t *output, uint32_t input_stride, TxType transform_type, uint8_t  bit_depth);

    void eb_av1_fwd_txfm2d_32x64_c(int16_t *input, int32_t *output, uint32_t input_stride, TxType transform_type, uint8_t  bit_depth);
    void eb_av1_fwd_txfm2d_32x64_avx2(int16_t *input, int32_t *output, uint32_t input_stride, TxType transform_type, uint8_t  bit_depth);
    RTCD_EXTERN void(*eb_av1_fwd_txfm2d_32x64)(int16_t *input, int32_t *output, uint32_t input_stride, TxType transform_type, uint8_t  bit_depth);

    void eb_av1_fwd_txfm2d_64x32_c(int16_t *input, int32_t *output, uint32_t input_stride, TxType transform_type, uint8_t  bit_depth);
    void eb_av1_fwd_txfm2d_64x32_avx2(int16_t *input, int32_t *output, uint32_t input_stride, TxType transform_type, uint8_t  bit_depth);
    RTCD_EXTERN void(*eb_av1_fwd_txfm2d_64x32)(int16_t *input, int32_t *output, uint32_t input_stride, TxType transform_type, uint8_t  bit_depth);

    void eb_av1_fwd_txfm2d_16x64_c(int16_t *input, int32_t *output, uint32_t input_stride, TxType transform_type, uint8_t  bit_depth);
    void eb_av1_fwd_txfm2d_16x64_avx2(int16_t *input, int32_t *output, uint32_t input_stride, TxType transform_type, uint8_t  bit_depth);
    RTCD_EXTERN void(*eb_av1_fwd_txfm2d_16x64)(int16_t *input, int32_t *output, uint32_t input_stride, TxType transform_type, uint8_t  bit_depth);

    void eb_av1_fwd_txfm2d_64x16_c(int16_t *input, int32_t *output, uint32_t input_stride, TxType transform_type, uint8_t  bit_depth);
    void eb_av1_fwd_txfm2d_64x16_avx2(int16_t *input, int32_t *output, uint32_t input_stride, TxType transform_type, uint8_t  bit_depth);
    RTCD_EXTERN void(*eb_av1_fwd_txfm2d_64x16)(int16_t *input, int32_t *output, uint32_t input_stride, TxType transform_type, uint8_t  bit_depth);

    void Av1TransformTwoD_64x64_c(int16_t *input, int32_t *output, uint32_t input_stride, TxType transform_type, uint8_t  bit_depth);
    void eb_av1_fwd_txfm2d_64x64_avx2(int16_t *input, int32_t *output, uint32_t input_stride, TxType transform_type, uint8_t  bit_depth);
    RTCD_EXTERN void(*eb_av1_fwd_txfm2d_64x64)(int16_t *input, int32_t *output, uint32_t input_stride, TxType transform_type, uint8_t  bit_depth);

    void Av1TransformTwoD_32x32_c(int16_t *input, int32_t *output, uint32_t input_stride, TxType transform_type, uint8_t  bit_depth);
    void eb_av1_fwd_txfm2d_32x32_avx2(int16_t *input, int32_t *output, uint32_t input_stride, TxType transform_type, uint8_t  bit_depth);
    RTCD_EXTERN void(*eb_av1_fwd_txfm2d_32x32)(int16_t *input, int32_t *output, uint32_t input_stride, TxType transform_type, uint8_t  bit_depth);

    void eb_av1_fwd_txfm2d_pf_32x32_c(int16_t *input, int32_t *output, uint32_t inputStride, TxType transform_type, uint8_t  bit_depth);
    void eb_av1_fwd_txfm2d_pf_32x32_avx2(int16_t *input, int32_t *output, uint32_t inputStride, TxType transform_type, uint8_t  bit_depth);
    #define eb_av1_fwd_txfm2d_pf_32x32 eb_av1_fwd_txfm2d_pf_32x32_avx2

    void Av1TransformTwoD_16x16_c(int16_t *input, int32_t *output, uint32_t input_stride, TxType transform_type, uint8_t  bit_depth);
    void eb_av1_fwd_txfm2d_16x16_avx2(int16_t *input, int32_t *output, uint32_t input_stride, TxType transform_type, uint8_t  bit_depth);
    RTCD_EXTERN void(*eb_av1_fwd_txfm2d_16x16)(int16_t *input, int32_t *output, uint32_t input_stride, TxType transform_type, uint8_t  bit_depth);

    void Av1TransformTwoD_8x8_c(int16_t *input, int32_t *output, uint32_t input_stride, TxType transform_type, uint8_t  bit_depth);
    void eb_av1_fwd_txfm2d_8x8_avx2(int16_t *input, int32_t *output, uint32_t input_stride, TxType transform_type, uint8_t  bit_depth);
    RTCD_EXTERN void(*eb_av1_fwd_txfm2d_8x8)(int16_t *input, int32_t *output, uint32_t input_stride, TxType transform_type, uint8_t  bit_depth);

    void Av1TransformTwoD_4x4_c(int16_t *input, int32_t *output, uint32_t input_stride, TxType transform_type, uint8_t  bit_depth);
    void eb_av1_fwd_txfm2d_4x4_sse4_1(int16_t *input, int32_t *output, uint32_t input_stride, TxType transform_type, uint8_t  bit_depth);
    RTCD_EXTERN void(*eb_av1_fwd_txfm2d_4x4)(int16_t *input, int32_t *output, uint32_t input_stride, TxType transform_type, uint8_t  bit_depth);

    void smooth_v_predictor_c(uint8_t *dst, ptrdiff_t stride, int32_t bw, int32_t bh, const uint8_t *above, const uint8_t *left);
    void eb_smooth_v_predictor_all_ssse3(uint8_t *dst, ptrdiff_t stride, int32_t bw, int32_t bh, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_smooth_v_predictor)(uint8_t *dst, ptrdiff_t stride, int32_t bw, int32_t bh, const uint8_t *above, const uint8_t *left);

    void smooth_h_predictor_c(uint8_t *dst, ptrdiff_t stride, int32_t bw, int32_t bh, const uint8_t *above, const uint8_t *left);
    void eb_smooth_h_predictor_all_ssse3(uint8_t *dst, ptrdiff_t stride, int32_t bw, int32_t bh, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_smooth_h_predictor)(uint8_t *dst, ptrdiff_t stride, int32_t bw, int32_t bh, const uint8_t *above, const uint8_t *left);

    void get_proj_subspace_c(const uint8_t *src8, int width, int height, int src_stride, const uint8_t *dat8, int dat_stride, int use_highbitdepth, int32_t *flt0, int flt0_stride, int32_t *flt1, int flt1_stride, int *xq, const SgrParamsType *params);
    void get_proj_subspace_avx2(const uint8_t *src8, int width, int height, int src_stride, const uint8_t *dat8, int dat_stride, int use_highbitdepth, int32_t *flt0, int flt0_stride, int32_t *flt1, int flt1_stride, int *xq, const SgrParamsType *params);
    RTCD_EXTERN void(*get_proj_subspace)(const uint8_t *src8, int width, int height, int src_stride, const uint8_t *dat8, int dat_stride, int use_highbitdepth, int32_t *flt0, int flt0_stride, int32_t *flt1, int flt1_stride, int *xq, const SgrParamsType *params);

    uint64_t HandleTransform16x64_c(int32_t *output);
    uint64_t HandleTransform16x64_avx2(int32_t *output);
    RTCD_EXTERN uint64_t(*HandleTransform16x64)(int32_t *output);

    uint64_t HandleTransform32x64_c(int32_t *output);
    uint64_t HandleTransform32x64_avx2(int32_t *output);
    RTCD_EXTERN uint64_t(*HandleTransform32x64)(int32_t *output);

    uint64_t HandleTransform64x16_c(int32_t *output);
    uint64_t HandleTransform64x16_avx2(int32_t *output);
    RTCD_EXTERN uint64_t(*HandleTransform64x16)(int32_t *output);

    uint64_t HandleTransform64x32_c(int32_t *output);
    uint64_t HandleTransform64x32_avx2(int32_t *output);
    RTCD_EXTERN uint64_t(*HandleTransform64x32)(int32_t *output);

    uint64_t HandleTransform64x64_c(int32_t *output);
    uint64_t HandleTransform64x64_avx2(int32_t *output);
    RTCD_EXTERN uint64_t(*HandleTransform64x64)(int32_t *output);

    uint64_t search_one_dual_c(int *lev0, int *lev1, int nb_strengths, uint64_t(**mse)[64], int sb_count, int fast, int start_gi, int end_gi);
    uint64_t search_one_dual_avx2(int *lev0, int *lev1, int nb_strengths, uint64_t(**mse)[64], int sb_count, int fast, int start_gi, int end_gi);
    RTCD_EXTERN uint64_t(*search_one_dual)(int *lev0, int *lev1, int nb_strengths, uint64_t(**mse)[64], int sb_count, int fast, int start_gi, int end_gi);

    uint32_t eb_aom_mse16x16_c(const uint8_t *src_ptr, int32_t  source_stride, const uint8_t *ref_ptr, int32_t  recon_stride, uint32_t *sse);
    uint32_t eb_aom_mse16x16_avx2(const uint8_t *src_ptr, int32_t  source_stride, const uint8_t *ref_ptr, int32_t  recon_stride, uint32_t *sse);
    RTCD_EXTERN uint32_t (*eb_aom_mse16x16)(const uint8_t *src_ptr, int32_t  source_stride, const uint8_t *ref_ptr, int32_t  recon_stride, uint32_t *sse);

    void eb_av1_convolve_2d_copy_sr_c(const uint8_t *src, int32_t src_stride, uint8_t *dst, int32_t dst_stride, int32_t w, int32_t h, InterpFilterParams *filter_params_x, InterpFilterParams *filter_params_y, const int32_t subpel_x_q4, const int32_t subpel_y_q4, ConvolveParams *conv_params);
    void eb_av1_convolve_2d_copy_sr_avx2(const uint8_t *src, int32_t src_stride, uint8_t *dst, int32_t dst_stride, int32_t w, int32_t h, InterpFilterParams *filter_params_x, InterpFilterParams *filter_params_y, const int32_t subpel_x_q4, const int32_t subpel_y_q4, ConvolveParams *conv_params);
    RTCD_EXTERN void(*eb_av1_convolve_2d_copy_sr)(const uint8_t *src, int32_t src_stride, uint8_t *dst, int32_t dst_stride, int32_t w, int32_t h, InterpFilterParams *filter_params_x, InterpFilterParams *filter_params_y, const int32_t subpel_x_q4, const int32_t subpel_y_q4, ConvolveParams *conv_params);

    void eb_av1_convolve_2d_sr_c(const uint8_t *src, int32_t src_stride, uint8_t *dst, int32_t dst_stride, int32_t w, int32_t h, InterpFilterParams *filter_params_x, InterpFilterParams *filter_params_y, const int32_t subpel_x_q4, const int32_t subpel_y_q4, ConvolveParams *conv_params);
    void eb_av1_convolve_2d_sr_avx2(const uint8_t *src, int32_t src_stride, uint8_t *dst, int32_t dst_stride, int32_t w, int32_t h, InterpFilterParams *filter_params_x, InterpFilterParams *filter_params_y, const int32_t subpel_x_q4, const int32_t subpel_y_q4, ConvolveParams *conv_params);
    RTCD_EXTERN void(*eb_av1_convolve_2d_sr)(const uint8_t *src, int32_t src_stride, uint8_t *dst, int32_t dst_stride, int32_t w, int32_t h, InterpFilterParams *filter_params_x, InterpFilterParams *filter_params_y, const int32_t subpel_x_q4, const int32_t subpel_y_q4, ConvolveParams *conv_params);

    void eb_av1_jnt_convolve_2d_copy_c(const uint8_t *src, int32_t src_stride, uint8_t *dst, int32_t dst_stride, int32_t w, int32_t h, InterpFilterParams *filter_params_x, InterpFilterParams *filter_params_y, const int32_t subpel_x_q4, const int32_t subpel_y_q4, ConvolveParams *conv_params);
    void eb_av1_jnt_convolve_2d_copy_avx2(const uint8_t *src, int32_t src_stride, uint8_t *dst, int32_t dst_stride, int32_t w, int32_t h, InterpFilterParams *filter_params_x, InterpFilterParams *filter_params_y, const int32_t subpel_x_q4, const int32_t subpel_y_q4, ConvolveParams *conv_params);
    RTCD_EXTERN void(*eb_av1_jnt_convolve_2d_copy)(const uint8_t *src, int32_t src_stride, uint8_t *dst, int32_t dst_stride, int32_t w, int32_t h, InterpFilterParams *filter_params_x, InterpFilterParams *filter_params_y, const int32_t subpel_x_q4, const int32_t subpel_y_q4, ConvolveParams *conv_params);

    void eb_av1_convolve_x_sr_c(const uint8_t *src, int32_t src_stride, uint8_t *dst, int32_t dst_stride, int32_t w, int32_t h, InterpFilterParams *filter_params_x, InterpFilterParams *filter_params_y, const int32_t subpel_x_q4, const int32_t subpel_y_q4, ConvolveParams *conv_params);
    void eb_av1_convolve_x_sr_avx2(const uint8_t *src, int32_t src_stride, uint8_t *dst, int32_t dst_stride, int32_t w, int32_t h, InterpFilterParams *filter_params_x, InterpFilterParams *filter_params_y, const int32_t subpel_x_q4, const int32_t subpel_y_q4, ConvolveParams *conv_params);
    RTCD_EXTERN void(*eb_av1_convolve_x_sr)(const uint8_t *src, int32_t src_stride, uint8_t *dst, int32_t dst_stride, int32_t w, int32_t h, InterpFilterParams *filter_params_x, InterpFilterParams *filter_params_y, const int32_t subpel_x_q4, const int32_t subpel_y_q4, ConvolveParams *conv_params);

    void eb_av1_convolve_y_sr_c(const uint8_t *src, int32_t src_stride, uint8_t *dst, int32_t dst_stride, int32_t w, int32_t h, InterpFilterParams *filter_params_x, InterpFilterParams *filter_params_y, const int32_t subpel_x_q4, const int32_t subpel_y_q4, ConvolveParams *conv_params);
    void eb_av1_convolve_y_sr_avx2(const uint8_t *src, int32_t src_stride, uint8_t *dst, int32_t dst_stride, int32_t w, int32_t h, InterpFilterParams *filter_params_x, InterpFilterParams *filter_params_y, const int32_t subpel_x_q4, const int32_t subpel_y_q4, ConvolveParams *conv_params);
    RTCD_EXTERN void(*eb_av1_convolve_y_sr)(const uint8_t *src, int32_t src_stride, uint8_t *dst, int32_t dst_stride, int32_t w, int32_t h, InterpFilterParams *filter_params_x, InterpFilterParams *filter_params_y, const int32_t subpel_x_q4, const int32_t subpel_y_q4, ConvolveParams *conv_params);

    void eb_av1_jnt_convolve_x_c(const uint8_t *src, int32_t src_stride, uint8_t *dst, int32_t dst_stride, int32_t w, int32_t h, InterpFilterParams *filter_params_x, InterpFilterParams *filter_params_y, const int32_t subpel_x_q4, const int32_t subpel_y_q4, ConvolveParams *conv_params);
    void eb_av1_jnt_convolve_x_avx2(const uint8_t *src, int32_t src_stride, uint8_t *dst, int32_t dst_stride, int32_t w, int32_t h, InterpFilterParams *filter_params_x, InterpFilterParams *filter_params_y, const int32_t subpel_x_q4, const int32_t subpel_y_q4, ConvolveParams *conv_params);
    RTCD_EXTERN void(*eb_av1_jnt_convolve_x)(const uint8_t *src, int32_t src_stride, uint8_t *dst, int32_t dst_stride, int32_t w, int32_t h, InterpFilterParams *filter_params_x, InterpFilterParams *filter_params_y, const int32_t subpel_x_q4, const int32_t subpel_y_q4, ConvolveParams *conv_params);

    void eb_av1_jnt_convolve_y_c(const uint8_t *src, int32_t src_stride, uint8_t *dst, int32_t dst_stride, int32_t w, int32_t h, InterpFilterParams *filter_params_x, InterpFilterParams *filter_params_y, const int32_t subpel_x_q4, const int32_t subpel_y_q4, ConvolveParams *conv_params);
    void eb_av1_jnt_convolve_y_avx2(const uint8_t *src, int32_t src_stride, uint8_t *dst, int32_t dst_stride, int32_t w, int32_t h, InterpFilterParams *filter_params_x, InterpFilterParams *filter_params_y, const int32_t subpel_x_q4, const int32_t subpel_y_q4, ConvolveParams *conv_params);
    RTCD_EXTERN void(*eb_av1_jnt_convolve_y)(const uint8_t *src, int32_t src_stride, uint8_t *dst, int32_t dst_stride, int32_t w, int32_t h, InterpFilterParams *filter_params_x, InterpFilterParams *filter_params_y, const int32_t subpel_x_q4, const int32_t subpel_y_q4, ConvolveParams *conv_params);

    void eb_av1_jnt_convolve_2d_c(const uint8_t *src, int32_t src_stride, uint8_t *dst, int32_t dst_stride, int32_t w, int32_t h, InterpFilterParams *filter_params_x, InterpFilterParams *filter_params_y, const int32_t subpel_x_q4, const int32_t subpel_y_q4, ConvolveParams *conv_params);
    void eb_av1_jnt_convolve_2d_avx2(const uint8_t *src, int32_t src_stride, uint8_t *dst, int32_t dst_stride, int32_t w, int32_t h, InterpFilterParams *filter_params_x, InterpFilterParams *filter_params_y, const int32_t subpel_x_q4, const int32_t subpel_y_q4, ConvolveParams *conv_params);
    RTCD_EXTERN void(*eb_av1_jnt_convolve_2d)(const uint8_t *src, int32_t src_stride, uint8_t *dst, int32_t dst_stride, int32_t w, int32_t h, InterpFilterParams *filter_params_x, InterpFilterParams *filter_params_y, const int32_t subpel_x_q4, const int32_t subpel_y_q4, ConvolveParams *conv_params);

    void eb_aom_quantize_b_c_II(const TranLow *coeff_ptr, intptr_t n_coeffs, int32_t skip_block, const int16_t *zbin_ptr, const int16_t *round_ptr, const int16_t *quant_ptr, const int16_t *quant_shift_ptr, TranLow *qcoeff_ptr, TranLow *dqcoeff_ptr, const int16_t *dequant_ptr, uint16_t *eob_ptr, const int16_t *scan, const int16_t *iscan);
    void eb_aom_highbd_quantize_b_avx2(const TranLow *coeff_ptr, intptr_t n_coeffs, int32_t skip_block, const int16_t *zbin_ptr, const int16_t *round_ptr, const int16_t *quant_ptr, const int16_t *quant_shift_ptr, TranLow *qcoeff_ptr, TranLow *dqcoeff_ptr, const int16_t *dequant_ptr, uint16_t *eob_ptr, const int16_t *scan, const int16_t *iscan);
    RTCD_EXTERN void(*eb_aom_quantize_b)(const TranLow *coeff_ptr, intptr_t n_coeffs, int32_t skip_block, const int16_t *zbin_ptr, const int16_t *round_ptr, const int16_t *quant_ptr, const int16_t *quant_shift_ptr, TranLow *qcoeff_ptr, TranLow *dqcoeff_ptr, const int16_t *dequant_ptr, uint16_t *eob_ptr, const int16_t *scan, const int16_t *iscan);

    void eb_aom_quantize_b_32x32_c_II(const TranLow *coeff_ptr, intptr_t n_coeffs, int32_t skip_block, const int16_t *zbin_ptr, const int16_t *round_ptr, const int16_t *quant_ptr, const int16_t *quant_shift_ptr, TranLow *qcoeff_ptr, TranLow *dqcoeff_ptr, const int16_t *dequant_ptr, uint16_t *eob_ptr, const int16_t *scan, const int16_t *iscan);
    void eb_aom_highbd_quantize_b_32x32_avx2(const TranLow *coeff_ptr, intptr_t n_coeffs, int skip_block, const int16_t *zbin_ptr, const int16_t *round_ptr, const int16_t *quant_ptr, const int16_t *quant_shift_ptr, TranLow *qcoeff_ptr, TranLow *dqcoeff_ptr, const int16_t *dequant_ptr, uint16_t *eob_ptr, const int16_t *scan, const int16_t *iscan);
    RTCD_EXTERN void(*eb_aom_quantize_b_32x32)(const TranLow *coeff_ptr, intptr_t n_coeffs, int32_t skip_block, const int16_t *zbin_ptr, const int16_t *round_ptr, const int16_t *quant_ptr, const int16_t *quant_shift_ptr, TranLow *qcoeff_ptr, TranLow *dqcoeff_ptr, const int16_t *dequant_ptr, uint16_t *eob_ptr, const int16_t *scan, const int16_t *iscan);

    void eb_aom_quantize_b_64x64_c_II(const TranLow *coeff_ptr, intptr_t n_coeffs, int32_t skip_block, const int16_t *zbin_ptr, const int16_t *round_ptr, const int16_t *quant_ptr, const int16_t *quant_shift_ptr, TranLow *qcoeff_ptr, TranLow *dqcoeff_ptr, const int16_t *dequant_ptr, uint16_t *eob_ptr, const int16_t *scan, const int16_t *iscan);
    void eb_aom_highbd_quantize_b_64x64_avx2(const TranLow *coeff_ptr, intptr_t n_coeffs, int32_t skip_block, const int16_t *zbin_ptr, const int16_t *round_ptr, const int16_t *quant_ptr, const int16_t *quant_shift_ptr, TranLow *qcoeff_ptr, TranLow *dqcoeff_ptr, const int16_t *dequant_ptr, uint16_t *eob_ptr, const int16_t *scan, const int16_t *iscan);
    RTCD_EXTERN void(*eb_aom_quantize_b_64x64)(const TranLow *coeff_ptr, intptr_t n_coeffs, int32_t skip_block, const int16_t *zbin_ptr, const int16_t *round_ptr, const int16_t *quant_ptr, const int16_t *quant_shift_ptr, TranLow *qcoeff_ptr, TranLow *dqcoeff_ptr, const int16_t *dequant_ptr, uint16_t *eob_ptr, const int16_t *scan, const int16_t *iscan);

    void eb_aom_highbd_quantize_b_c(const TranLow *coeff_ptr, intptr_t n_coeffs, int32_t skip_block, const int16_t *zbin_ptr, const int16_t *round_ptr, const int16_t *quant_ptr, const int16_t *quant_shift_ptr, TranLow *qcoeff_ptr, TranLow *dqcoeff_ptr, const int16_t *dequant_ptr, uint16_t *eob_ptr, const int16_t *scan, const int16_t *iscan);
    RTCD_EXTERN void(*eb_aom_highbd_quantize_b)(const TranLow *coeff_ptr, intptr_t n_coeffs, int32_t skip_block, const int16_t *zbin_ptr, const int16_t *round_ptr, const int16_t *quant_ptr, const int16_t *quant_shift_ptr, TranLow *qcoeff_ptr, TranLow *dqcoeff_ptr, const int16_t *dequant_ptr, uint16_t *eob_ptr, const int16_t *scan, const int16_t *iscan);

    void eb_aom_highbd_quantize_b_32x32_c(const TranLow *coeff_ptr, intptr_t n_coeffs, int32_t skip_block, const int16_t *zbin_ptr, const int16_t *round_ptr, const int16_t *quant_ptr, const int16_t *quant_shift_ptr, TranLow *qcoeff_ptr, TranLow *dqcoeff_ptr, const int16_t *dequant_ptr, uint16_t *eob_ptr, const int16_t *scan, const int16_t *iscan);
    RTCD_EXTERN void(*eb_aom_highbd_quantize_b_32x32)(const TranLow *coeff_ptr, intptr_t n_coeffs, int32_t skip_block, const int16_t *zbin_ptr, const int16_t *round_ptr, const int16_t *quant_ptr, const int16_t *quant_shift_ptr, TranLow *qcoeff_ptr, TranLow *dqcoeff_ptr, const int16_t *dequant_ptr, uint16_t *eob_ptr, const int16_t *scan, const int16_t *iscan);

    void eb_aom_highbd_quantize_b_64x64_c(const TranLow *coeff_ptr, intptr_t n_coeffs, int32_t skip_block, const int16_t *zbin_ptr, const int16_t *round_ptr, const int16_t *quant_ptr, const int16_t *quant_shift_ptr, TranLow *qcoeff_ptr, TranLow *dqcoeff_ptr, const int16_t *dequant_ptr, uint16_t *eob_ptr, const int16_t *scan, const int16_t *iscan);
    RTCD_EXTERN void(*eb_aom_highbd_quantize_b_64x64)(const TranLow *coeff_ptr, intptr_t n_coeffs, int32_t skip_block, const int16_t *zbin_ptr, const int16_t *round_ptr, const int16_t *quant_ptr, const int16_t *quant_shift_ptr, TranLow *qcoeff_ptr, TranLow *dqcoeff_ptr, const int16_t *dequant_ptr, uint16_t *eob_ptr, const int16_t *scan, const int16_t *iscan);

    void eb_av1_highbd_warp_affine_c(const int32_t *mat, const uint16_t *ref, int width, int height, int stride, uint16_t *pred, int p_col, int p_row, int p_width, int p_height, int p_stride, int subsampling_x, int subsampling_y, int bd, ConvolveParams *conv_params, int16_t alpha, int16_t beta, int16_t gamma, int16_t delta);

    void eb_av1_inv_txfm2d_add_4x4_c(const int32_t *input, uint16_t *output, int32_t stride, TxType tx_type, int32_t bd);
    void eb_av1_inv_txfm2d_add_4x4_sse4_1(const int32_t *input, uint16_t *output, int32_t stride, TxType tx_type, int32_t bd);
    void eb_av1_inv_txfm2d_add_4x4_avx2(const int32_t *input, uint16_t *output, int32_t stride, TxType tx_type, int32_t bd);
    RTCD_EXTERN void(*eb_av1_inv_txfm2d_add_4x4)(const int32_t *input, uint16_t *output, int32_t stride, TxType tx_type, int32_t bd);

    void eb_av1_inv_txfm2d_add_8x8_c(const int32_t *input, uint16_t *output, int32_t stride, TxType tx_type, int32_t bd);
    void eb_av1_inv_txfm2d_add_8x8_sse4_1(const int32_t *input, uint16_t *output, int32_t stride, TxType tx_type, int32_t bd);
    void eb_av1_inv_txfm2d_add_8x8_avx2(const int32_t *input, uint16_t *output, int32_t stride, TxType tx_type, int32_t bd);
    RTCD_EXTERN void(*eb_av1_inv_txfm2d_add_8x8)(const int32_t *input, uint16_t *output, int32_t stride, TxType tx_type, int32_t bd);

    void eb_av1_inv_txfm2d_add_16x16_c(const int32_t *input, uint16_t *output, int32_t stride, TxType tx_type, int32_t bd);
    void eb_av1_inv_txfm2d_add_16x16_sse4_1(const int32_t *input, uint16_t *output, int32_t stride, TxType tx_type, int32_t bd);
    void eb_av1_inv_txfm2d_add_16x16_avx2(const int32_t *input, uint16_t *output, int32_t stride, TxType tx_type, int32_t bd);
    RTCD_EXTERN void(*eb_av1_inv_txfm2d_add_16x16)(const int32_t *input, uint16_t *output, int32_t stride, TxType tx_type, int32_t bd);

    void eb_av1_inv_txfm2d_add_32x32_c(const int32_t *input, uint16_t *output, int32_t stride, TxType tx_type, int32_t bd);
    void eb_av1_inv_txfm2d_add_32x32_avx2(const int32_t *input, uint16_t *output, int32_t stride, TxType tx_type, int32_t bd);
    RTCD_EXTERN void(*eb_av1_inv_txfm2d_add_32x32)(const int32_t *input, uint16_t *output, int32_t stride, TxType tx_type, int32_t bd);

    void eb_av1_inv_txfm2d_add_64x64_c(const int32_t *input, uint16_t *output, int32_t stride, TxType tx_type, int32_t bd);
    void eb_av1_inv_txfm2d_add_64x64_sse4_1(const int32_t *input, uint16_t *output, int32_t stride, TxType tx_type, int32_t bd);
    RTCD_EXTERN void(*eb_av1_inv_txfm2d_add_64x64)(const int32_t *input, uint16_t *output, int32_t stride, TxType tx_type, int32_t bd);

    void eb_av1_highbd_inv_txfm_add_avx2(const int32_t *input, uint16_t *output, int32_t stride, TxType tx_type, TxSize tx_size, int32_t eob, int32_t bd);

    void eb_av1_inv_txfm2d_add_8x16_c(const int32_t *input, uint16_t *output, int32_t stride, TxType tx_type, TxSize tx_size, int32_t eob, int32_t bd);
    RTCD_EXTERN void(*eb_av1_inv_txfm2d_add_8x16)(const int32_t *input, uint16_t *output, int32_t stride, const TxType tx_type, TxSize tx_size, int32_t eob, int32_t bd);

    void eb_av1_inv_txfm2d_add_16x8_c(const int32_t *input, uint16_t *output, int32_t stride, TxType tx_type, TxSize tx_size, int32_t eob, int32_t bd);
    RTCD_EXTERN void(*eb_av1_inv_txfm2d_add_16x8)(const int32_t *input, uint16_t *output, int32_t stride, TxType tx_type, TxSize tx_size, int32_t eob, int32_t bd);

    void eb_av1_inv_txfm2d_add_16x32_c(const int32_t *input, uint16_t *output, int32_t stride, TxType tx_type, TxSize tx_size, int32_t eob, int32_t bd);
    RTCD_EXTERN void(*eb_av1_inv_txfm2d_add_16x32)(const int32_t *input, uint16_t *output, int32_t stride, TxType tx_type, TxSize tx_size, int32_t eob, int32_t bd);

    void eb_av1_inv_txfm2d_add_32x16_c(const int32_t *input, uint16_t *output, int32_t stride, TxType tx_type, TxSize tx_size, int32_t eob, int32_t bd);
    RTCD_EXTERN void(*eb_av1_inv_txfm2d_add_32x16)(const int32_t *input, uint16_t *output, int32_t stride, TxType tx_type, TxSize tx_size, int32_t eob, int32_t bd);

    void eb_av1_inv_txfm2d_add_32x8_c(const int32_t *input, uint16_t *output, int32_t stride, TxType tx_type, TxSize tx_size, int32_t eob, int32_t bd);
    RTCD_EXTERN void(*eb_av1_inv_txfm2d_add_32x8)(const int32_t *input, uint16_t *output, int32_t stride, TxType tx_type, TxSize tx_size, int32_t eob, int32_t bd);

    void eb_av1_inv_txfm2d_add_8x32_c(const int32_t *input, uint16_t *output, int32_t stride, TxType tx_type, TxSize tx_size, int32_t eob, int32_t bd);
    RTCD_EXTERN void(*eb_av1_inv_txfm2d_add_8x32)(const int32_t *input, uint16_t *output, int32_t stride, TxType tx_type, TxSize tx_size, int32_t eob, int32_t bd);

    void eb_av1_inv_txfm2d_add_32x64_c(const int32_t *input, uint16_t *output, int32_t stride, TxType tx_type, TxSize tx_size, int32_t eob, int32_t bd);
    RTCD_EXTERN void(*eb_av1_inv_txfm2d_add_32x64)(const int32_t *input, uint16_t *output, int32_t stride, TxType tx_type, TxSize tx_size, int32_t eob, int32_t bd);

    void eb_av1_inv_txfm2d_add_64x32_c(const int32_t *input, uint16_t *output, int32_t stride, TxType tx_type, TxSize tx_size, int32_t eob, int32_t bd);
    RTCD_EXTERN void(*eb_av1_inv_txfm2d_add_64x32)(const int32_t *input, uint16_t *output, int32_t stride, TxType tx_type, TxSize tx_size, int32_t eob, int32_t bd);

    void eb_av1_inv_txfm2d_add_16x64_c(const int32_t *input, uint16_t *output, int32_t stride, TxType tx_type, TxSize tx_size, int32_t eob, int32_t bd);
    RTCD_EXTERN void(*eb_av1_inv_txfm2d_add_16x64)(const int32_t *input, uint16_t *output, int32_t stride, TxType tx_type, TxSize tx_size, int32_t eob, int32_t bd);

    void eb_av1_inv_txfm2d_add_64x16_c(const int32_t *input, uint16_t *output, int32_t stride, TxType tx_type, TxSize tx_size, int32_t eob, int32_t bd);
    RTCD_EXTERN void(*eb_av1_inv_txfm2d_add_64x16)(const int32_t *input, uint16_t *output, int32_t stride, TxType tx_type, TxSize tx_size, int32_t eob, int32_t bd);

    void eb_av1_inv_txfm2d_add_4x8_c(const int32_t *input, uint16_t *output, int32_t stride, TxType tx_type, TxSize tx_size, int32_t bd);
    void eb_av1_inv_txfm2d_add_4x8_sse4_1(const int32_t *input, uint16_t *output, int32_t stride, TxType tx_type, TxSize tx_size, int32_t bd);
    RTCD_EXTERN void(*eb_av1_inv_txfm2d_add_4x8)(const int32_t *input, uint16_t *output, int32_t stride, TxType tx_type, TxSize tx_size, int32_t bd);

    void eb_av1_inv_txfm2d_add_8x4_c(const int32_t *input, uint16_t *output, int32_t stride, TxType tx_type, TxSize tx_size, int32_t bd);
    void eb_av1_inv_txfm2d_add_8x4_sse4_1(const int32_t *input, uint16_t *output, int32_t stride, TxType tx_type, TxSize tx_size, int32_t bd);
    RTCD_EXTERN void(*eb_av1_inv_txfm2d_add_8x4)(const int32_t *input, uint16_t *output, int32_t stride, TxType tx_type, TxSize tx_size, int32_t bd);

    void eb_av1_inv_txfm2d_add_4x16_c(const int32_t *input, uint16_t *output, int32_t stride, TxType tx_type, TxSize tx_size, int32_t bd);
    void eb_av1_inv_txfm2d_add_4x16_sse4_1(const int32_t *input, uint16_t *output, int32_t stride, TxType tx_type, TxSize tx_size, int32_t bd);
    RTCD_EXTERN void(*eb_av1_inv_txfm2d_add_4x16)(const int32_t *input, uint16_t *output, int32_t stride, TxType tx_type, TxSize tx_size, int32_t bd);

    void eb_av1_inv_txfm2d_add_16x4_c(const int32_t *input, uint16_t *output, int32_t stride, TxType tx_type, TxSize tx_size, int32_t bd);
    void eb_av1_inv_txfm2d_add_16x4_sse4_1(const int32_t *input, uint16_t *output, int32_t stride, TxType tx_type, TxSize tx_size, int32_t bd);
    RTCD_EXTERN void(*eb_av1_inv_txfm2d_add_16x4)(const int32_t *input, uint16_t *output, int32_t stride, TxType tx_type, TxSize tx_size, int32_t bd);

    void eb_av1_inv_txfm_add_c(const TranLow *dqcoeff, uint8_t *dst, int32_t stride, const TxfmParam *txfm_param);
    void eb_av1_inv_txfm_add_ssse3(const TranLow *dqcoeff, uint8_t *dst, int32_t stride, const TxfmParam *txfm_param);
    void eb_av1_inv_txfm_add_avx2(const TranLow *dqcoeff, uint8_t *dst, int32_t stride, const TxfmParam *txfm_param);
    RTCD_EXTERN void(*eb_av1_inv_txfm_add)(const TranLow *dqcoeff, uint8_t *dst, int32_t stride, const TxfmParam *txfm_param);

    void eb_av1_quantize_fp_c(const TranLow *coeff_ptr, intptr_t n_coeffs, const int16_t *zbin_ptr, const int16_t *round_ptr, const int16_t *quant_ptr, const int16_t *quant_shift_ptr, TranLow *qcoeff_ptr, TranLow *dqcoeff_ptr, const int16_t *dequant_ptr, uint16_t *eob_ptr, const int16_t *scan, const int16_t *iscan);
    void eb_av1_quantize_fp_avx2(const TranLow *coeff_ptr, intptr_t n_coeffs, const int16_t *zbin_ptr, const int16_t *round_ptr, const int16_t *quant_ptr, const int16_t *quant_shift_ptr, TranLow *qcoeff_ptr, TranLow *dqcoeff_ptr, const int16_t *dequant_ptr, uint16_t *eob_ptr, const int16_t *scan, const int16_t *iscan);
    RTCD_EXTERN void (*eb_av1_quantize_fp)(const TranLow *coeff_ptr, intptr_t n_coeffs, const int16_t *zbin_ptr, const int16_t *round_ptr, const int16_t *quant_ptr, const int16_t *quant_shift_ptr, TranLow *qcoeff_ptr, TranLow *dqcoeff_ptr, const int16_t *dequant_ptr, uint16_t *eob_ptr, const int16_t *scan, const int16_t *iscan);

    void eb_av1_quantize_fp_32x32_c(const TranLow *coeff_ptr, intptr_t n_coeffs, const int16_t *zbin_ptr, const int16_t *round_ptr, const int16_t *quant_ptr, const int16_t *quant_shift_ptr, TranLow *qcoeff_ptr, TranLow *dqcoeff_ptr, const int16_t *dequant_ptr, uint16_t *eob_ptr, const int16_t *scan, const int16_t *iscan);
    void eb_av1_quantize_fp_32x32_avx2(const TranLow *coeff_ptr, intptr_t n_coeffs, const int16_t *zbin_ptr, const int16_t *round_ptr, const int16_t *quant_ptr, const int16_t *quant_shift_ptr, TranLow *qcoeff_ptr, TranLow *dqcoeff_ptr, const int16_t *dequant_ptr, uint16_t *eob_ptr, const int16_t *scan, const int16_t *iscan);
    RTCD_EXTERN void(*eb_av1_quantize_fp_32x32)(const TranLow *coeff_ptr, intptr_t n_coeffs, const int16_t *zbin_ptr, const int16_t *round_ptr, const int16_t *quant_ptr, const int16_t *quant_shift_ptr, TranLow *qcoeff_ptr, TranLow *dqcoeff_ptr, const int16_t *dequant_ptr, uint16_t *eob_ptr, const int16_t *scan, const int16_t *iscan);

    void eb_av1_quantize_fp_64x64_c(const TranLow *coeff_ptr, intptr_t n_coeffs, const int16_t *zbin_ptr, const int16_t *round_ptr, const int16_t *quant_ptr, const int16_t *quant_shift_ptr, TranLow *qcoeff_ptr, TranLow *dqcoeff_ptr, const int16_t *dequant_ptr, uint16_t *eob_ptr, const int16_t *scan, const int16_t *iscan);
    void eb_av1_quantize_fp_64x64_avx2(const TranLow *coeff_ptr, intptr_t n_coeffs, const int16_t *zbin_ptr, const int16_t *round_ptr, const int16_t *quant_ptr, const int16_t *quant_shift_ptr, TranLow *qcoeff_ptr, TranLow *dqcoeff_ptr, const int16_t *dequant_ptr, uint16_t *eob_ptr, const int16_t *scan, const int16_t *iscan);
    RTCD_EXTERN void(*eb_av1_quantize_fp_64x64)(const TranLow *coeff_ptr, intptr_t n_coeffs, const int16_t *zbin_ptr, const int16_t *round_ptr, const int16_t *quant_ptr, const int16_t *quant_shift_ptr, TranLow *qcoeff_ptr, TranLow *dqcoeff_ptr, const int16_t *dequant_ptr, uint16_t *eob_ptr, const int16_t *scan, const int16_t *iscan);

    //uint32_t eb_aom_highbd_8_mse16x16_c(const uint8_t *src_ptr, int32_t  source_stride, const uint8_t *ref_ptr, int32_t  recon_stride, uint32_t *sse);
    void      eb_aom_highbd_8_mse16x16_sse2(const uint8_t *src_ptr, int32_t  source_stride, const uint8_t *ref_ptr, int32_t  recon_stride, uint32_t *sse);
    RTCD_EXTERN void(*eb_aom_highbd_8_mse16x16)(const uint8_t *src_ptr, int32_t  source_stride, const uint8_t *ref_ptr, int32_t  recon_stride, uint32_t *sse);

    void eb_av1_upsample_intra_edge_c(uint8_t *p, int32_t sz);
    void eb_av1_upsample_intra_edge_sse4_1(uint8_t *p, int32_t sz);
    RTCD_EXTERN void(*eb_av1_upsample_intra_edge)(uint8_t *p, int32_t sz);

    // AMIR
    void eb_av1_upsample_intra_edge_high_c(uint16_t *p, int32_t sz, int32_t bd);
    //void eb_av1_upsample_intra_edge_high_sse4_1(uint16_t *p, int32_t sz, int32_t bd);
    //RTCD_EXTERN void(*eb_av1_upsample_intra_edge_high)(uint16_t *p, int32_t sz, int32_t bd);
/* DC_PRED top */

    void eb_av1_warp_affine_c(const int32_t *mat, const uint8_t *ref, int width, int height, int stride, uint8_t *pred, int p_col, int p_row, int p_width, int p_height, int p_stride, int subsampling_x, int subsampling_y, ConvolveParams *conv_params, int16_t alpha, int16_t beta, int16_t gamma, int16_t delta);
    void eb_av1_warp_affine_avx2(const int32_t *mat, const uint8_t *ref, int width, int height, int stride, uint8_t *pred, int p_col, int p_row, int p_width, int p_height, int p_stride, int subsampling_x, int subsampling_y, ConvolveParams *conv_params, int16_t alpha, int16_t beta, int16_t gamma, int16_t delta);
    RTCD_EXTERN void (*eb_av1_warp_affine)(const int32_t *mat, const uint8_t *ref, int width, int height, int stride, uint8_t *pred, int p_col, int p_row, int p_width, int p_height, int p_stride, int subsampling_x, int subsampling_y, ConvolveParams *conv_params, int16_t alpha, int16_t beta, int16_t gamma, int16_t delta);

    void eb_aom_dc_predictor_4x4_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_dc_predictor_4x4_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_dc_predictor_4x4)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_dc_predictor_8x8_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_dc_predictor_8x8_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_dc_predictor_8x8)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_dc_predictor_16x16_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_dc_predictor_16x16_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_dc_predictor_16x16)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_dc_predictor_32x32_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_dc_predictor_32x32_avx2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_dc_predictor_32x32)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_dc_predictor_64x64_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_dc_predictor_64x64_avx2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_dc_predictor_64x64)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_dc_predictor_16x32_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_dc_predictor_16x32_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_dc_predictor_16x32)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_dc_predictor_16x4_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_dc_predictor_16x4_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_dc_predictor_16x4)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_dc_predictor_16x64_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_dc_predictor_16x64_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_dc_predictor_16x64)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_dc_predictor_16x8_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_dc_predictor_16x8_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_dc_predictor_16x8)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_dc_predictor_32x16_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_dc_predictor_32x16_avx2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_dc_predictor_32x16)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_dc_predictor_32x64_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_dc_predictor_32x64_avx2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_dc_predictor_32x64)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_dc_predictor_32x8_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_dc_predictor_32x8_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_dc_predictor_32x8)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_dc_predictor_4x16_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_dc_predictor_4x16_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_dc_predictor_4x16)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_dc_predictor_4x8_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_dc_predictor_4x8_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_dc_predictor_4x8)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_dc_predictor_64x16_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_dc_predictor_64x16_avx2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_dc_predictor_64x16)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_dc_predictor_64x32_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_dc_predictor_64x32_avx2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_dc_predictor_64x32)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_dc_predictor_8x16_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_dc_predictor_8x16_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_dc_predictor_8x16)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_dc_predictor_8x32_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_dc_predictor_8x32_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_dc_predictor_8x32)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_dc_predictor_8x4_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_dc_predictor_8x4_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_dc_predictor_8x4)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    /* DC_PRED top */

    void eb_aom_dc_left_predictor_4x4_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_dc_left_predictor_4x4_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_dc_left_predictor_4x4)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_dc_left_predictor_8x8_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_dc_left_predictor_8x8_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_dc_left_predictor_8x8)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_dc_left_predictor_16x16_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_dc_left_predictor_16x16_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_dc_left_predictor_16x16)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_dc_left_predictor_32x32_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_dc_left_predictor_32x32_avx2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_dc_left_predictor_32x32)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_dc_left_predictor_64x64_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_dc_left_predictor_64x64_avx2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_dc_left_predictor_64x64)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_dc_left_predictor_16x32_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_dc_left_predictor_16x32_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_dc_left_predictor_16x32)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_dc_left_predictor_16x4_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_dc_left_predictor_16x4_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_dc_left_predictor_16x4)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_dc_left_predictor_16x64_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_dc_left_predictor_16x64_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_dc_left_predictor_16x64)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_dc_left_predictor_16x8_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_dc_left_predictor_16x8_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_dc_left_predictor_16x8)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_dc_left_predictor_32x16_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_dc_left_predictor_32x16_avx2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_dc_left_predictor_32x16)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_dc_left_predictor_32x64_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_dc_left_predictor_32x64_avx2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_dc_left_predictor_32x64)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_dc_left_predictor_32x8_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_dc_left_predictor_32x8_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_dc_left_predictor_32x8)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_dc_left_predictor_4x16_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_dc_left_predictor_4x16_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_dc_left_predictor_4x16)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_dc_left_predictor_4x8_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_dc_left_predictor_4x8_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_dc_left_predictor_4x8)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_dc_left_predictor_64x16_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_dc_left_predictor_64x16_avx2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_dc_left_predictor_64x16)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_dc_left_predictor_64x32_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_dc_left_predictor_64x32_avx2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_dc_left_predictor_64x32)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_dc_left_predictor_8x16_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_dc_left_predictor_8x16_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_dc_left_predictor_8x16)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_dc_left_predictor_8x32_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_dc_left_predictor_8x32_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_dc_left_predictor_8x32)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_dc_left_predictor_8x4_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_dc_left_predictor_8x4_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_dc_left_predictor_8x4)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    /* DC_PRED top */

    void eb_aom_dc_top_predictor_4x4_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_dc_top_predictor_4x4_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_dc_top_predictor_4x4)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_dc_top_predictor_8x8_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_dc_top_predictor_8x8_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_dc_top_predictor_8x8)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_dc_top_predictor_16x16_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_dc_top_predictor_16x16_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_dc_top_predictor_16x16)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_dc_top_predictor_32x32_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_dc_top_predictor_32x32_avx2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_dc_top_predictor_32x32)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_dc_top_predictor_64x64_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_dc_top_predictor_64x64_avx2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_dc_top_predictor_64x64)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_dc_top_predictor_16x32_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_dc_top_predictor_16x32_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_dc_top_predictor_16x32)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_dc_top_predictor_16x4_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_dc_top_predictor_16x4_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_dc_top_predictor_16x4)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_dc_top_predictor_16x64_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_dc_top_predictor_16x64_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_dc_top_predictor_16x64)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_dc_top_predictor_16x8_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_dc_top_predictor_16x8_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_dc_top_predictor_16x8)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_dc_top_predictor_32x16_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_dc_top_predictor_32x16_avx2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_dc_top_predictor_32x16)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_dc_top_predictor_32x64_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_dc_top_predictor_32x64_avx2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_dc_top_predictor_32x64)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_dc_top_predictor_32x8_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_dc_top_predictor_32x8_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_dc_top_predictor_32x8)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_dc_top_predictor_4x16_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_dc_top_predictor_4x16_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_dc_top_predictor_4x16)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_dc_top_predictor_4x8_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_dc_top_predictor_4x8_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_dc_top_predictor_4x8)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_dc_top_predictor_64x16_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_dc_top_predictor_64x16_avx2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_dc_top_predictor_64x16)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_dc_top_predictor_64x32_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_dc_top_predictor_64x32_avx2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_dc_top_predictor_64x32)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_dc_top_predictor_8x16_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_dc_top_predictor_8x16_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_dc_top_predictor_8x16)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_dc_top_predictor_8x32_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_dc_top_predictor_8x32_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_dc_top_predictor_8x32)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_dc_top_predictor_8x4_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_dc_top_predictor_8x4_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_dc_top_predictor_8x4)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    /* DC_PRED 128 */
    void eb_aom_dc_128_predictor_4x4_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_dc_128_predictor_4x4_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_dc_128_predictor_4x4)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_dc_128_predictor_8x8_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_dc_128_predictor_8x8_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_dc_128_predictor_8x8)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_dc_128_predictor_16x16_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_dc_128_predictor_16x16_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_dc_128_predictor_16x16)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_dc_128_predictor_32x32_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_dc_128_predictor_32x32_avx2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_dc_128_predictor_32x32)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_dc_128_predictor_64x64_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_dc_128_predictor_64x64_avx2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_dc_128_predictor_64x64)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_dc_128_predictor_16x32_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_dc_128_predictor_16x32_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_dc_128_predictor_16x32)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_dc_128_predictor_16x4_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_dc_128_predictor_16x4_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_dc_128_predictor_16x4)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_dc_128_predictor_16x64_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_dc_128_predictor_16x64_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_dc_128_predictor_16x64)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_dc_128_predictor_16x8_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_dc_128_predictor_16x8_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_dc_128_predictor_16x8)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_dc_128_predictor_32x16_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_dc_128_predictor_32x16_avx2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_dc_128_predictor_32x16)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_dc_128_predictor_32x64_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_dc_128_predictor_32x64_avx2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_dc_128_predictor_32x64)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_dc_128_predictor_32x8_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_dc_128_predictor_32x8_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_dc_128_predictor_32x8)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_dc_128_predictor_4x16_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_dc_128_predictor_4x16_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_dc_128_predictor_4x16)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_dc_128_predictor_4x8_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_dc_128_predictor_4x8_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_dc_128_predictor_4x8)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_dc_128_predictor_64x16_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_dc_128_predictor_64x16_avx2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_dc_128_predictor_64x16)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_dc_128_predictor_64x32_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_dc_128_predictor_64x32_avx2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_dc_128_predictor_64x32)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_dc_128_predictor_8x16_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_dc_128_predictor_8x16_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_dc_128_predictor_8x16)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_dc_128_predictor_8x32_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_dc_128_predictor_8x32_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_dc_128_predictor_8x32)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_dc_128_predictor_8x4_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_dc_128_predictor_8x4_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_dc_128_predictor_8x4)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    /* SMOOTH_H_PRED */
    void eb_aom_smooth_h_predictor_64x64_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_smooth_h_predictor_64x64_ssse3(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_smooth_h_predictor_64x64)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_smooth_h_predictor_32x32_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_smooth_h_predictor_32x32_ssse3(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_smooth_h_predictor_32x32)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_smooth_h_predictor_16x16_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_smooth_h_predictor_16x16_ssse3(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_smooth_h_predictor_16x16)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_smooth_h_predictor_8x8_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_smooth_h_predictor_8x8_ssse3(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_smooth_h_predictor_8x8)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_smooth_h_predictor_4x4_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_smooth_h_predictor_4x4_ssse3(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_smooth_h_predictor_4x4)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_smooth_h_predictor_16x32_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_smooth_h_predictor_16x32_ssse3(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_smooth_h_predictor_16x32)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_smooth_h_predictor_16x4_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_smooth_h_predictor_16x4_ssse3(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_smooth_h_predictor_16x4)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_smooth_h_predictor_16x64_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_smooth_h_predictor_16x64_ssse3(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_smooth_h_predictor_16x64)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_smooth_h_predictor_16x8_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_smooth_h_predictor_16x8_ssse3(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_smooth_h_predictor_16x8)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_smooth_h_predictor_32x16_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_smooth_h_predictor_32x16_ssse3(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_smooth_h_predictor_32x16)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_smooth_h_predictor_32x64_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_smooth_h_predictor_32x64_ssse3(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_smooth_h_predictor_32x64)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_smooth_h_predictor_32x8_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_smooth_h_predictor_32x8_ssse3(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_smooth_h_predictor_32x8)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_smooth_h_predictor_4x16_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_smooth_h_predictor_4x16_ssse3(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_smooth_h_predictor_4x16)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_smooth_h_predictor_4x8_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_smooth_h_predictor_4x8_ssse3(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_smooth_h_predictor_4x8)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_smooth_h_predictor_64x16_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_smooth_h_predictor_64x16_ssse3(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_smooth_h_predictor_64x16)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_smooth_h_predictor_64x32_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_smooth_h_predictor_64x32_ssse3(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_smooth_h_predictor_64x32)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_smooth_h_predictor_8x16_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_smooth_h_predictor_8x16_ssse3(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_smooth_h_predictor_8x16)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_smooth_h_predictor_8x32_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_smooth_h_predictor_8x32_ssse3(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_smooth_h_predictor_8x32)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_smooth_h_predictor_8x4_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_smooth_h_predictor_8x4_ssse3(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_smooth_h_predictor_8x4)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    /* SMOOTH_V_PRED */

    void eb_aom_smooth_v_predictor_64x64_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_smooth_v_predictor_64x64_ssse3(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_smooth_v_predictor_64x64)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_smooth_v_predictor_32x32_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_smooth_v_predictor_32x32_ssse3(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_smooth_v_predictor_32x32)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_smooth_v_predictor_16x16_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_smooth_v_predictor_16x16_ssse3(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_smooth_v_predictor_16x16)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_smooth_v_predictor_8x8_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_smooth_v_predictor_8x8_ssse3(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_smooth_v_predictor_8x8)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_smooth_v_predictor_4x4_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_smooth_v_predictor_4x4_ssse3(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_smooth_v_predictor_4x4)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_smooth_v_predictor_16x32_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_smooth_v_predictor_16x32_ssse3(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_smooth_v_predictor_16x32)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_smooth_v_predictor_16x4_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_smooth_v_predictor_16x4_ssse3(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_smooth_v_predictor_16x4)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_smooth_v_predictor_16x64_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_smooth_v_predictor_16x64_ssse3(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_smooth_v_predictor_16x64)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_smooth_v_predictor_16x8_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_smooth_v_predictor_16x8_ssse3(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_smooth_v_predictor_16x8)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_smooth_v_predictor_32x16_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_smooth_v_predictor_32x16_ssse3(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_smooth_v_predictor_32x16)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_smooth_v_predictor_32x64_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_smooth_v_predictor_32x64_ssse3(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_smooth_v_predictor_32x64)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_smooth_v_predictor_32x8_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_smooth_v_predictor_32x8_ssse3(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_smooth_v_predictor_32x8)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_smooth_v_predictor_4x16_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_smooth_v_predictor_4x16_ssse3(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_smooth_v_predictor_4x16)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_smooth_v_predictor_4x8_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_smooth_v_predictor_4x8_ssse3(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_smooth_v_predictor_4x8)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_smooth_v_predictor_64x16_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_smooth_v_predictor_64x16_ssse3(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_smooth_v_predictor_64x16)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_smooth_v_predictor_64x32_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_smooth_v_predictor_64x32_ssse3(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_smooth_v_predictor_64x32)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_smooth_v_predictor_8x16_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_smooth_v_predictor_8x16_ssse3(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_smooth_v_predictor_8x16)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_smooth_v_predictor_8x32_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_smooth_v_predictor_8x32_ssse3(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_smooth_v_predictor_8x32)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_smooth_v_predictor_8x4_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_smooth_v_predictor_8x4_ssse3(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_smooth_v_predictor_8x4)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    /* SMOOTH_PRED */

    void eb_aom_smooth_predictor_64x64_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_smooth_predictor_64x64_ssse3(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_smooth_predictor_64x64)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_smooth_predictor_32x32_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_smooth_predictor_32x32_ssse3(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_smooth_predictor_32x32)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_smooth_predictor_16x16_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_smooth_predictor_16x16_ssse3(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_smooth_predictor_16x16)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_smooth_predictor_8x8_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_smooth_predictor_8x8_ssse3(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_smooth_predictor_8x8)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_smooth_predictor_4x4_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_smooth_predictor_4x4_ssse3(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_smooth_predictor_4x4)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_smooth_predictor_16x32_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_smooth_predictor_16x32_ssse3(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_smooth_predictor_16x32)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_smooth_predictor_16x4_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_smooth_predictor_16x4_ssse3(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_smooth_predictor_16x4)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_smooth_predictor_16x64_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_smooth_predictor_16x64_ssse3(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_smooth_predictor_16x64)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_smooth_predictor_16x8_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_smooth_predictor_16x8_ssse3(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_smooth_predictor_16x8)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_smooth_predictor_32x16_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_smooth_predictor_32x16_ssse3(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_smooth_predictor_32x16)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_smooth_predictor_32x64_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_smooth_predictor_32x64_ssse3(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_smooth_predictor_32x64)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_smooth_predictor_32x8_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_smooth_predictor_32x8_ssse3(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_smooth_predictor_32x8)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_smooth_predictor_4x16_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_smooth_predictor_4x16_ssse3(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_smooth_predictor_4x16)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_smooth_predictor_4x8_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_smooth_predictor_4x8_ssse3(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_smooth_predictor_4x8)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_smooth_predictor_64x16_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_smooth_predictor_64x16_ssse3(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_smooth_predictor_64x16)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_smooth_predictor_64x32_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_smooth_predictor_64x32_ssse3(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_smooth_predictor_64x32)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_smooth_predictor_8x16_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_smooth_predictor_8x16_ssse3(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_smooth_predictor_8x16)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_smooth_predictor_8x32_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_smooth_predictor_8x32_ssse3(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_smooth_predictor_8x32)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_smooth_predictor_8x4_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_smooth_predictor_8x4_ssse3(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_smooth_predictor_8x4)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    /* V_PRED */
    void eb_aom_v_predictor_4x4_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_v_predictor_4x4_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_v_predictor_4x4)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_v_predictor_8x8_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_v_predictor_8x8_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_v_predictor_8x8)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_v_predictor_16x16_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_v_predictor_16x16_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_v_predictor_16x16)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_v_predictor_32x32_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_v_predictor_32x32_avx2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_v_predictor_32x32)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_v_predictor_64x64_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_v_predictor_64x64_avx2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_v_predictor_64x64)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_v_predictor_16x32_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_v_predictor_16x32_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_v_predictor_16x32)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_v_predictor_16x4_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_v_predictor_16x4_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_v_predictor_16x4)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_v_predictor_16x64_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_v_predictor_16x64_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_v_predictor_16x64)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_v_predictor_16x8_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_v_predictor_16x8_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_v_predictor_16x8)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_v_predictor_32x16_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_v_predictor_32x16_avx2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_v_predictor_32x16)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_v_predictor_32x64_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_v_predictor_32x64_avx2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_v_predictor_32x64)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_v_predictor_32x8_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_v_predictor_32x8_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_v_predictor_32x8)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_v_predictor_4x16_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_v_predictor_4x16_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_v_predictor_4x16)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_v_predictor_4x8_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_v_predictor_4x8_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_v_predictor_4x8)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_v_predictor_64x16_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_v_predictor_64x16_avx2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_v_predictor_64x16)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_v_predictor_64x32_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_v_predictor_64x32_avx2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_v_predictor_64x32)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_v_predictor_8x16_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_v_predictor_8x16_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_v_predictor_8x16)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_v_predictor_8x32_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_v_predictor_8x32_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_v_predictor_8x32)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_v_predictor_8x4_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_v_predictor_8x4_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_v_predictor_8x4)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    /* H_PRED */

    void eb_aom_h_predictor_4x4_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_h_predictor_4x4_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_h_predictor_4x4)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_h_predictor_8x8_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_h_predictor_8x8_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_h_predictor_8x8)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_h_predictor_16x16_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_h_predictor_16x16_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_h_predictor_16x16)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_h_predictor_32x32_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_h_predictor_32x32_avx2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_h_predictor_32x32)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_h_predictor_64x64_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_h_predictor_64x64_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_h_predictor_64x64)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_h_predictor_16x32_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_h_predictor_16x32_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_h_predictor_16x32)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_h_predictor_16x4_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_h_predictor_16x4_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_h_predictor_16x4)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_h_predictor_16x64_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_h_predictor_16x64_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_h_predictor_16x64)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_h_predictor_16x8_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_h_predictor_16x8_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_h_predictor_16x8)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_h_predictor_32x16_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_h_predictor_32x16_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_h_predictor_32x16)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_h_predictor_32x64_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_h_predictor_32x64_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_h_predictor_32x64)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_h_predictor_32x8_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_h_predictor_32x8_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_h_predictor_32x8)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_h_predictor_4x16_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_h_predictor_4x16_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_h_predictor_4x16)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_h_predictor_4x8_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_h_predictor_4x8_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_h_predictor_4x8)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_h_predictor_64x16_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_h_predictor_64x16_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_h_predictor_64x16)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_h_predictor_64x32_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_h_predictor_64x32_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_h_predictor_64x32)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_h_predictor_8x16_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_h_predictor_8x16_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_h_predictor_8x16)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_h_predictor_8x32_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_h_predictor_8x32_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_h_predictor_8x32)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_h_predictor_8x4_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_h_predictor_8x4_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_aom_h_predictor_8x4)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_paeth_predictor_16x16_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_paeth_predictor_16x16_ssse3(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_paeth_predictor_16x16_avx2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void (*eb_aom_paeth_predictor_16x16)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_paeth_predictor_16x32_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_paeth_predictor_16x32_ssse3(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_paeth_predictor_16x32_avx2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void (*eb_aom_paeth_predictor_16x32)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_paeth_predictor_16x4_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_paeth_predictor_16x4_ssse3(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void (*eb_aom_paeth_predictor_16x4)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_paeth_predictor_16x64_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_paeth_predictor_16x64_ssse3(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_paeth_predictor_16x64_avx2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void (*eb_aom_paeth_predictor_16x64)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_paeth_predictor_16x8_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_paeth_predictor_16x8_ssse3(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_paeth_predictor_16x8_avx2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void (*eb_aom_paeth_predictor_16x8)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_paeth_predictor_2x2_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    #define eb_aom_paeth_predictor_2x2 eb_aom_paeth_predictor_2x2_c

    void eb_aom_paeth_predictor_32x16_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_paeth_predictor_32x16_ssse3(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_paeth_predictor_32x16_avx2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void (*eb_aom_paeth_predictor_32x16)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_paeth_predictor_32x32_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_paeth_predictor_32x32_ssse3(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_paeth_predictor_32x32_avx2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void (*eb_aom_paeth_predictor_32x32)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_paeth_predictor_32x64_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_paeth_predictor_32x64_ssse3(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_paeth_predictor_32x64_avx2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void (*eb_aom_paeth_predictor_32x64)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_paeth_predictor_32x8_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_paeth_predictor_32x8_ssse3(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void (*eb_aom_paeth_predictor_32x8)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_paeth_predictor_4x16_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_paeth_predictor_4x16_ssse3(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void (*eb_aom_paeth_predictor_4x16)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_paeth_predictor_4x4_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_paeth_predictor_4x4_ssse3(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void (*eb_aom_paeth_predictor_4x4)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_paeth_predictor_4x8_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_paeth_predictor_4x8_ssse3(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void (*eb_aom_paeth_predictor_4x8)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_paeth_predictor_64x16_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_paeth_predictor_64x16_ssse3(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_paeth_predictor_64x16_avx2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void (*eb_aom_paeth_predictor_64x16)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_paeth_predictor_64x32_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_paeth_predictor_64x32_ssse3(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_paeth_predictor_64x32_avx2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void (*eb_aom_paeth_predictor_64x32)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_paeth_predictor_64x64_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_paeth_predictor_64x64_ssse3(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_paeth_predictor_64x64_avx2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void (*eb_aom_paeth_predictor_64x64)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_paeth_predictor_8x16_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_paeth_predictor_8x16_ssse3(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void (*eb_aom_paeth_predictor_8x16)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_paeth_predictor_8x32_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_paeth_predictor_8x32_ssse3(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void (*eb_aom_paeth_predictor_8x32)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_paeth_predictor_8x4_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_paeth_predictor_8x4_ssse3(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void (*eb_aom_paeth_predictor_8x4)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_paeth_predictor_8x8_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void eb_aom_paeth_predictor_8x8_ssse3(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void (*eb_aom_paeth_predictor_8x8)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void eb_aom_highbd_paeth_predictor_16x16_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int bd);
    void eb_aom_highbd_paeth_predictor_16x16_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int bd);
    RTCD_EXTERN void(*eb_aom_highbd_paeth_predictor_16x16)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int bd);

    void eb_aom_highbd_paeth_predictor_16x32_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int bd);
    void eb_aom_highbd_paeth_predictor_16x32_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int bd);
    RTCD_EXTERN void(*eb_aom_highbd_paeth_predictor_16x32)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int bd);

    void eb_aom_highbd_paeth_predictor_16x4_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int bd);
    void eb_aom_highbd_paeth_predictor_16x4_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int bd);
    RTCD_EXTERN void(*eb_aom_highbd_paeth_predictor_16x4)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int bd);

    void eb_aom_highbd_paeth_predictor_16x64_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int bd);
    void eb_aom_highbd_paeth_predictor_16x64_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int bd);
    RTCD_EXTERN void(*eb_aom_highbd_paeth_predictor_16x64)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int bd);

    void eb_aom_highbd_paeth_predictor_16x8_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int bd);
    void eb_aom_highbd_paeth_predictor_16x8_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int bd);
    RTCD_EXTERN void(*eb_aom_highbd_paeth_predictor_16x8)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int bd);

    void eb_aom_highbd_paeth_predictor_2x2_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int bd);
    void eb_aom_highbd_paeth_predictor_2x2_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int bd);
    RTCD_EXTERN void(*eb_aom_highbd_paeth_predictor_2x2)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int bd);

    void eb_aom_highbd_paeth_predictor_32x16_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int bd);
    void eb_aom_highbd_paeth_predictor_32x16_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int bd);
    RTCD_EXTERN void(*eb_aom_highbd_paeth_predictor_32x16)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int bd);

    void eb_aom_highbd_paeth_predictor_32x32_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int bd);
    void eb_aom_highbd_paeth_predictor_32x32_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int bd);
    RTCD_EXTERN void(*eb_aom_highbd_paeth_predictor_32x32)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int bd);

    void eb_aom_highbd_paeth_predictor_32x64_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int bd);
    void eb_aom_highbd_paeth_predictor_32x64_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int bd);
    RTCD_EXTERN void(*eb_aom_highbd_paeth_predictor_32x64)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int bd);

    void eb_aom_highbd_paeth_predictor_32x8_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int bd);
    void eb_aom_highbd_paeth_predictor_32x8_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int bd);
    RTCD_EXTERN void(*eb_aom_highbd_paeth_predictor_32x8)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int bd);

    void eb_aom_highbd_paeth_predictor_4x16_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int bd);
    void eb_aom_highbd_paeth_predictor_4x16_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int bd);
    RTCD_EXTERN void(*eb_aom_highbd_paeth_predictor_4x16)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int bd);

    void eb_aom_highbd_paeth_predictor_4x4_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int bd);
    void eb_aom_highbd_paeth_predictor_4x4_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int bd);
    RTCD_EXTERN void(*eb_aom_highbd_paeth_predictor_4x4)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int bd);

    void eb_aom_highbd_paeth_predictor_4x8_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int bd);
    void eb_aom_highbd_paeth_predictor_4x8_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int bd);
    RTCD_EXTERN void(*eb_aom_highbd_paeth_predictor_4x8)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int bd);

    void eb_aom_highbd_paeth_predictor_64x16_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int bd);
    void eb_aom_highbd_paeth_predictor_64x16_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int bd);
    RTCD_EXTERN void(*eb_aom_highbd_paeth_predictor_64x16)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int bd);

    void eb_aom_highbd_paeth_predictor_64x32_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int bd);
    void eb_aom_highbd_paeth_predictor_64x32_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int bd);
    RTCD_EXTERN void(*eb_aom_highbd_paeth_predictor_64x32)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int bd);

    void eb_aom_highbd_paeth_predictor_64x64_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int bd);
    void eb_aom_highbd_paeth_predictor_64x64_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int bd);
    RTCD_EXTERN void(*eb_aom_highbd_paeth_predictor_64x64)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int bd);

    void eb_aom_highbd_paeth_predictor_8x16_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int bd);
    void eb_aom_highbd_paeth_predictor_8x16_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int bd);
    RTCD_EXTERN void(*eb_aom_highbd_paeth_predictor_8x16)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int bd);

    void eb_aom_highbd_paeth_predictor_8x32_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int bd);
    void eb_aom_highbd_paeth_predictor_8x32_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int bd);
    RTCD_EXTERN void(*eb_aom_highbd_paeth_predictor_8x32)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int bd);

    void eb_aom_highbd_paeth_predictor_8x4_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int bd);
    void eb_aom_highbd_paeth_predictor_8x4_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int bd);
    RTCD_EXTERN void(*eb_aom_highbd_paeth_predictor_8x4)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int bd);

    void eb_aom_highbd_paeth_predictor_8x8_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int bd);
    void eb_aom_highbd_paeth_predictor_8x8_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int bd);
    RTCD_EXTERN void(*eb_aom_highbd_paeth_predictor_8x8)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int bd);

    uint32_t eb_aom_sad128x128_c(const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride);
    uint32_t eb_aom_sad128x128_avx2(const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride);
    RTCD_EXTERN uint32_t(*eb_aom_sad128x128)(const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride);

    void eb_aom_sad128x128x4d_c(const uint8_t *src_ptr, int src_stride, const uint8_t * const ref_ptr[], int ref_stride, uint32_t *sad_array);
    void eb_aom_sad128x128x4d_avx2(const uint8_t *src_ptr, int src_stride, const uint8_t * const ref_ptr[], int ref_stride, uint32_t *sad_array);
    RTCD_EXTERN void(*eb_aom_sad128x128x4d)(const uint8_t *src_ptr, int src_stride, const uint8_t * const ref_ptr[], int ref_stride, uint32_t *sad_array);

    uint32_t eb_aom_sad128x64_c(const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride);
    uint32_t eb_aom_sad128x64_avx2(const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride);
    RTCD_EXTERN uint32_t(*eb_aom_sad128x64)(const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride);

    void eb_aom_sad128x64x4d_c(const uint8_t *src_ptr, int src_stride, const uint8_t * const ref_ptr[], int ref_stride, uint32_t *sad_array);
    void eb_aom_sad128x64x4d_avx2(const uint8_t *src_ptr, int src_stride, const uint8_t * const ref_ptr[], int ref_stride, uint32_t *sad_array);
    RTCD_EXTERN void(*eb_aom_sad128x64x4d)(const uint8_t *src_ptr, int src_stride, const uint8_t * const ref_ptr[], int ref_stride, uint32_t *sad_array);

    uint32_t eb_aom_sad16x16_c(const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride);
    uint32_t eb_aom_sad16x16_avx2(const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride);
    RTCD_EXTERN uint32_t(*eb_aom_sad16x16)(const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride);

    void eb_aom_sad16x16x4d_c(const uint8_t *src_ptr, int src_stride, const uint8_t * const ref_ptr[], int ref_stride, uint32_t *sad_array);
    void eb_aom_sad16x16x4d_avx2(const uint8_t *src_ptr, int src_stride, const uint8_t * const ref_ptr[], int ref_stride, uint32_t *sad_array);
    RTCD_EXTERN void(*eb_aom_sad16x16x4d)(const uint8_t *src_ptr, int src_stride, const uint8_t * const ref_ptr[], int ref_stride, uint32_t *sad_array);

    uint32_t eb_aom_sad16x32_c(const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride);
    uint32_t eb_aom_sad16x32_avx2(const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride);
    RTCD_EXTERN uint32_t(*eb_aom_sad16x32)(const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride);

    void eb_aom_sad16x32x4d_c(const uint8_t *src_ptr, int src_stride, const uint8_t * const ref_ptr[], int ref_stride, uint32_t *sad_array);
    void eb_aom_sad16x32x4d_avx2(const uint8_t *src_ptr, int src_stride, const uint8_t * const ref_ptr[], int ref_stride, uint32_t *sad_array);
    RTCD_EXTERN void(*eb_aom_sad16x32x4d)(const uint8_t *src_ptr, int src_stride, const uint8_t * const ref_ptr[], int ref_stride, uint32_t *sad_array);

    uint32_t eb_aom_sad16x4_c(const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride);
    uint32_t eb_aom_sad16x4_avx2(const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride);
    RTCD_EXTERN uint32_t(*eb_aom_sad16x4)(const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride);

    void eb_aom_sad16x4x4d_c(const uint8_t *src_ptr, int src_stride, const uint8_t * const ref_ptr[], int ref_stride, uint32_t *sad_array);
    void eb_aom_sad16x4x4d_avx2(const uint8_t *src_ptr, int src_stride, const uint8_t * const ref_ptr[], int ref_stride, uint32_t *sad_array);
    RTCD_EXTERN void(*eb_aom_sad16x4x4d)(const uint8_t *src_ptr, int src_stride, const uint8_t * const ref_ptr[], int ref_stride, uint32_t *sad_array);

    uint32_t eb_aom_sad16x64_c(const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride);
    uint32_t eb_aom_sad16x64_avx2(const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride);
    RTCD_EXTERN uint32_t(*eb_aom_sad16x64)(const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride);

    void eb_aom_sad16x64x4d_c(const uint8_t *src_ptr, int src_stride, const uint8_t * const ref_ptr[], int ref_stride, uint32_t *sad_array);
    void eb_aom_sad16x64x4d_avx2(const uint8_t *src_ptr, int src_stride, const uint8_t * const ref_ptr[], int ref_stride, uint32_t *sad_array);
    RTCD_EXTERN void(*eb_aom_sad16x64x4d)(const uint8_t *src_ptr, int src_stride, const uint8_t * const ref_ptr[], int ref_stride, uint32_t *sad_array);

    uint32_t eb_aom_sad16x8_c(const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride);
    uint32_t eb_aom_sad16x8_avx2(const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride);
    RTCD_EXTERN uint32_t(*eb_aom_sad16x8)(const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride);

    void eb_aom_sad16x8x4d_c(const uint8_t *src_ptr, int src_stride, const uint8_t * const ref_ptr[], int ref_stride, uint32_t *sad_array);
    void eb_aom_sad16x8x4d_avx2(const uint8_t *src_ptr, int src_stride, const uint8_t * const ref_ptr[], int ref_stride, uint32_t *sad_array);
    RTCD_EXTERN void(*eb_aom_sad16x8x4d)(const uint8_t *src_ptr, int src_stride, const uint8_t * const ref_ptr[], int ref_stride, uint32_t *sad_array);

    uint32_t eb_aom_sad32x16_c(const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride);
    uint32_t eb_aom_sad32x16_avx2(const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride);
    RTCD_EXTERN uint32_t(*eb_aom_sad32x16)(const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride);

    void eb_aom_sad32x16x4d_c(const uint8_t *src_ptr, int src_stride, const uint8_t * const ref_ptr[], int ref_stride, uint32_t *sad_array);
    void eb_aom_sad32x16x4d_avx2(const uint8_t *src_ptr, int src_stride, const uint8_t * const ref_ptr[], int ref_stride, uint32_t *sad_array);
    RTCD_EXTERN void(*eb_aom_sad32x16x4d)(const uint8_t *src_ptr, int src_stride, const uint8_t * const ref_ptr[], int ref_stride, uint32_t *sad_array);

    uint32_t eb_aom_sad32x32_c(const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride);
    uint32_t eb_aom_sad32x32_avx2(const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride);
    RTCD_EXTERN uint32_t(*eb_aom_sad32x32)(const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride);

    void eb_aom_sad32x32x4d_c(const uint8_t *src_ptr, int src_stride, const uint8_t * const ref_ptr[], int ref_stride, uint32_t *sad_array);
    void eb_aom_sad32x32x4d_avx2(const uint8_t *src_ptr, int src_stride, const uint8_t * const ref_ptr[], int ref_stride, uint32_t *sad_array);
    RTCD_EXTERN void(*eb_aom_sad32x32x4d)(const uint8_t *src_ptr, int src_stride, const uint8_t * const ref_ptr[], int ref_stride, uint32_t *sad_array);

    uint32_t eb_aom_sad32x64_c(const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride);
    uint32_t eb_aom_sad32x64_avx2(const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride);
    RTCD_EXTERN uint32_t(*eb_aom_sad32x64)(const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride);

    void eb_aom_sad32x64x4d_c(const uint8_t *src_ptr, int src_stride, const uint8_t * const ref_ptr[], int ref_stride, uint32_t *sad_array);
    void eb_aom_sad32x64x4d_avx2(const uint8_t *src_ptr, int src_stride, const uint8_t * const ref_ptr[], int ref_stride, uint32_t *sad_array);
    RTCD_EXTERN void(*eb_aom_sad32x64x4d)(const uint8_t *src_ptr, int src_stride, const uint8_t * const ref_ptr[], int ref_stride, uint32_t *sad_array);

    uint32_t eb_aom_sad32x8_c(const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride);
    uint32_t eb_aom_sad32x8_avx2(const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride);
    RTCD_EXTERN uint32_t(*eb_aom_sad32x8)(const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride);

    void eb_aom_sad32x8x4d_c(const uint8_t *src_ptr, int src_stride, const uint8_t * const ref_ptr[], int ref_stride, uint32_t *sad_array);
    void eb_aom_sad32x8x4d_avx2(const uint8_t *src_ptr, int src_stride, const uint8_t * const ref_ptr[], int ref_stride, uint32_t *sad_array);
    RTCD_EXTERN void(*eb_aom_sad32x8x4d)(const uint8_t *src_ptr, int src_stride, const uint8_t * const ref_ptr[], int ref_stride, uint32_t *sad_array);

    uint32_t eb_aom_sad4x16_c(const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride);
    uint32_t eb_aom_sad4x16_avx2(const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride);
    RTCD_EXTERN uint32_t(*eb_aom_sad4x16)(const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride);

    void eb_aom_sad4x16x4d_c(const uint8_t *src_ptr, int src_stride, const uint8_t * const ref_ptr[], int ref_stride, uint32_t *sad_array);
    void eb_aom_sad4x16x4d_avx2(const uint8_t *src_ptr, int src_stride, const uint8_t * const ref_ptr[], int ref_stride, uint32_t *sad_array);
    RTCD_EXTERN void(*eb_aom_sad4x16x4d)(const uint8_t *src_ptr, int src_stride, const uint8_t * const ref_ptr[], int ref_stride, uint32_t *sad_array);

    uint32_t eb_aom_sad4x4_c(const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride);
    uint32_t eb_aom_sad4x4_avx2(const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride);
    RTCD_EXTERN uint32_t(*eb_aom_sad4x4)(const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride);

    void eb_aom_sad4x4x4d_c(const uint8_t *src_ptr, int src_stride, const uint8_t * const ref_ptr[], int ref_stride, uint32_t *sad_array);
    void eb_aom_sad4x4x4d_avx2(const uint8_t *src_ptr, int src_stride, const uint8_t * const ref_ptr[], int ref_stride, uint32_t *sad_array);
    RTCD_EXTERN void(*eb_aom_sad4x4x4d)(const uint8_t *src_ptr, int src_stride, const uint8_t * const ref_ptr[], int ref_stride, uint32_t *sad_array);

    uint32_t eb_aom_sad4x8_c(const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride);
    uint32_t eb_aom_sad4x8_avx2(const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride);
    RTCD_EXTERN uint32_t(*eb_aom_sad4x8)(const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride);

    void eb_aom_sad4x8x4d_c(const uint8_t *src_ptr, int src_stride, const uint8_t * const ref_ptr[], int ref_stride, uint32_t *sad_array);
    void eb_aom_sad4x8x4d_avx2(const uint8_t *src_ptr, int src_stride, const uint8_t * const ref_ptr[], int ref_stride, uint32_t *sad_array);
    RTCD_EXTERN void(*eb_aom_sad4x8x4d)(const uint8_t *src_ptr, int src_stride, const uint8_t * const ref_ptr[], int ref_stride, uint32_t *sad_array);

    uint32_t eb_aom_sad64x128_c(const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride);
    uint32_t eb_aom_sad64x128_avx2(const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride);
    RTCD_EXTERN uint32_t(*eb_aom_sad64x128)(const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride);

    void eb_aom_sad64x128x4d_c(const uint8_t *src_ptr, int src_stride, const uint8_t * const ref_ptr[], int ref_stride, uint32_t *sad_array);
    void eb_aom_sad64x128x4d_avx2(const uint8_t *src_ptr, int src_stride, const uint8_t * const ref_ptr[], int ref_stride, uint32_t *sad_array);
    RTCD_EXTERN void(*eb_aom_sad64x128x4d)(const uint8_t *src_ptr, int src_stride, const uint8_t * const ref_ptr[], int ref_stride, uint32_t *sad_array);

    uint32_t eb_aom_sad64x16_c(const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride);
    uint32_t eb_aom_sad64x16_avx2(const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride);
    RTCD_EXTERN uint32_t(*eb_aom_sad64x16)(const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride);

    void eb_aom_sad64x16x4d_c(const uint8_t *src_ptr, int src_stride, const uint8_t * const ref_ptr[], int ref_stride, uint32_t *sad_array);
    void eb_aom_sad64x16x4d_avx2(const uint8_t *src_ptr, int src_stride, const uint8_t * const ref_ptr[], int ref_stride, uint32_t *sad_array);
    RTCD_EXTERN void(*eb_aom_sad64x16x4d)(const uint8_t *src_ptr, int src_stride, const uint8_t * const ref_ptr[], int ref_stride, uint32_t *sad_array);

    uint32_t eb_aom_sad64x32_c(const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride);
    uint32_t eb_aom_sad64x32_avx2(const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride);
    RTCD_EXTERN uint32_t(*eb_aom_sad64x32)(const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride);

    void eb_aom_sad64x32x4d_c(const uint8_t *src_ptr, int src_stride, const uint8_t * const ref_ptr[], int ref_stride, uint32_t *sad_array);
    void eb_aom_sad64x32x4d_avx2(const uint8_t *src_ptr, int src_stride, const uint8_t * const ref_ptr[], int ref_stride, uint32_t *sad_array);
    RTCD_EXTERN void(*eb_aom_sad64x32x4d)(const uint8_t *src_ptr, int src_stride, const uint8_t * const ref_ptr[], int ref_stride, uint32_t *sad_array);

    uint32_t eb_aom_sad64x64_c(const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride);
    uint32_t eb_aom_sad64x64_avx2(const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride);
    RTCD_EXTERN uint32_t(*eb_aom_sad64x64)(const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride);

    void eb_aom_sad64x64x4d_c(const uint8_t *src_ptr, int src_stride, const uint8_t * const ref_ptr[], int ref_stride, uint32_t *sad_array);
    void eb_aom_sad64x64x4d_avx2(const uint8_t *src_ptr, int src_stride, const uint8_t * const ref_ptr[], int ref_stride, uint32_t *sad_array);
    RTCD_EXTERN void(*eb_aom_sad64x64x4d)(const uint8_t *src_ptr, int src_stride, const uint8_t * const ref_ptr[], int ref_stride, uint32_t *sad_array);

    uint32_t eb_aom_sad8x16_c(const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride);
    uint32_t eb_aom_sad8x16_avx2(const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride);
    RTCD_EXTERN uint32_t(*eb_aom_sad8x16)(const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride);

    void eb_aom_sad8x16x4d_c(const uint8_t *src_ptr, int src_stride, const uint8_t * const ref_ptr[], int ref_stride, uint32_t *sad_array);
    void eb_aom_sad8x16x4d_avx2(const uint8_t *src_ptr, int src_stride, const uint8_t * const ref_ptr[], int ref_stride, uint32_t *sad_array);
    RTCD_EXTERN void(*eb_aom_sad8x16x4d)(const uint8_t *src_ptr, int src_stride, const uint8_t * const ref_ptr[], int ref_stride, uint32_t *sad_array);

    uint32_t eb_aom_sad8x32_c(const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride);
    uint32_t eb_aom_sad8x32_avx2(const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride);
    RTCD_EXTERN uint32_t(*eb_aom_sad8x32)(const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride);

    void eb_aom_sad8x32x4d_c(const uint8_t *src_ptr, int src_stride, const uint8_t * const ref_ptr[], int ref_stride, uint32_t *sad_array);
    void eb_aom_sad8x32x4d_avx2(const uint8_t *src_ptr, int src_stride, const uint8_t * const ref_ptr[], int ref_stride, uint32_t *sad_array);
    RTCD_EXTERN void(*eb_aom_sad8x32x4d)(const uint8_t *src_ptr, int src_stride, const uint8_t * const ref_ptr[], int ref_stride, uint32_t *sad_array);

    uint32_t eb_aom_sad8x4_c(const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride);
    uint32_t eb_aom_sad8x4_avx2(const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride);
    RTCD_EXTERN uint32_t(*eb_aom_sad8x4)(const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride);

    void eb_aom_sad8x4x4d_c(const uint8_t *src_ptr, int src_stride, const uint8_t * const ref_ptr[], int ref_stride, uint32_t *sad_array);
    void eb_aom_sad8x4x4d_avx2(const uint8_t *src_ptr, int src_stride, const uint8_t * const ref_ptr[], int ref_stride, uint32_t *sad_array);
    RTCD_EXTERN void(*eb_aom_sad8x4x4d)(const uint8_t *src_ptr, int src_stride, const uint8_t * const ref_ptr[], int ref_stride, uint32_t *sad_array);

    uint32_t eb_aom_sad8x8_c(const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride);
    uint32_t eb_aom_sad8x8_avx2(const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride);
    RTCD_EXTERN uint32_t(*eb_aom_sad8x8)(const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride);

    void eb_aom_sad8x8x4d_c(const uint8_t *src_ptr, int src_stride, const uint8_t * const ref_ptr[], int ref_stride, uint32_t *sad_array);
    void eb_aom_sad8x8x4d_avx2(const uint8_t *src_ptr, int src_stride, const uint8_t * const ref_ptr[], int ref_stride, uint32_t *sad_array);
    RTCD_EXTERN void(*eb_aom_sad8x8x4d)(const uint8_t *src_ptr, int src_stride, const uint8_t * const ref_ptr[], int ref_stride, uint32_t *sad_array);

    unsigned int eb_aom_variance4x4_c(const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse);
    unsigned int eb_aom_variance4x4_sse2(const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse);
    RTCD_EXTERN unsigned int(*eb_aom_variance4x4)(const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse);

    unsigned int eb_aom_variance4x8_c(const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse);
    unsigned int eb_aom_variance4x8_sse2(const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse);
    RTCD_EXTERN unsigned int(*eb_aom_variance4x8)(const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse);

    unsigned int eb_aom_variance4x16_c(const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse);
    unsigned int eb_aom_variance4x16_sse2(const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse);
    RTCD_EXTERN unsigned int(*eb_aom_variance4x16)(const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse);

    unsigned int eb_aom_variance8x4_c(const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse);
    unsigned int eb_aom_variance8x4_sse2(const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse);
    RTCD_EXTERN unsigned int(*eb_aom_variance8x4)(const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse);

    unsigned int eb_aom_variance8x8_c(const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse);
    unsigned int eb_aom_variance8x8_sse2(const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse);
    RTCD_EXTERN unsigned int(*eb_aom_variance8x8)(const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse);

    unsigned int eb_aom_variance8x16_c(const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse);
    unsigned int eb_aom_variance8x16_sse2(const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse);
    RTCD_EXTERN unsigned int(*eb_aom_variance8x16)(const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse);

    unsigned int eb_aom_variance8x32_c(const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse);
    unsigned int eb_aom_variance8x32_sse2(const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse);
    RTCD_EXTERN unsigned int(*eb_aom_variance8x32)(const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse);

    unsigned int eb_aom_variance16x4_c(const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse);
    unsigned int eb_aom_variance16x4_avx2(const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse);
    RTCD_EXTERN unsigned int(*eb_aom_variance16x4)(const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse);

    unsigned int eb_aom_variance16x8_c(const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse);
    unsigned int eb_aom_variance16x8_avx2(const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse);
    RTCD_EXTERN unsigned int(*eb_aom_variance16x8)(const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse);

    unsigned int eb_aom_variance16x16_c(const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse);
    unsigned int eb_aom_variance16x16_avx2(const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse);
    RTCD_EXTERN unsigned int(*eb_aom_variance16x16)(const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse);

    unsigned int eb_aom_variance16x32_c(const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse);
    unsigned int eb_aom_variance16x32_avx2(const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse);
    RTCD_EXTERN unsigned int(*eb_aom_variance16x32)(const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse);

    unsigned int eb_aom_variance16x64_c(const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse);
    unsigned int eb_aom_variance16x64_avx2(const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse);
    RTCD_EXTERN unsigned int(*eb_aom_variance16x64)(const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse);

    unsigned int eb_aom_variance32x8_c(const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse);
    unsigned int eb_aom_variance32x8_avx2(const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse);
    RTCD_EXTERN unsigned int(*eb_aom_variance32x8)(const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse);

    unsigned int eb_aom_variance32x16_c(const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse);
    unsigned int eb_aom_variance32x16_avx2(const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse);
    RTCD_EXTERN unsigned int(*eb_aom_variance32x16)(const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse);

    unsigned int eb_aom_variance32x32_c(const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse);
    unsigned int eb_aom_variance32x32_avx2(const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse);
    RTCD_EXTERN unsigned int(*eb_aom_variance32x32)(const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse);

    unsigned int eb_aom_variance32x64_c(const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse);
    unsigned int eb_aom_variance32x64_avx2(const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse);
    RTCD_EXTERN unsigned int(*eb_aom_variance32x64)(const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse);

    unsigned int eb_aom_variance64x16_c(const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse);
    unsigned int eb_aom_variance64x16_avx2(const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse);
    RTCD_EXTERN unsigned int(*eb_aom_variance64x16)(const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse);

    unsigned int eb_aom_variance64x32_c(const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse);
    unsigned int eb_aom_variance64x32_avx2(const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse);
    RTCD_EXTERN unsigned int(*eb_aom_variance64x32)(const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse);

    unsigned int eb_aom_variance64x64_c(const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse);
    unsigned int eb_aom_variance64x64_avx2(const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse);
    RTCD_EXTERN unsigned int(*eb_aom_variance64x64)(const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse);

    unsigned int eb_aom_variance64x128_c(const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse);
    unsigned int eb_aom_variance64x128_avx2(const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse);
    RTCD_EXTERN unsigned int(*eb_aom_variance64x128)(const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse);

    unsigned int eb_aom_variance128x64_c(const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse);
    unsigned int eb_aom_variance128x64_avx2(const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse);
    RTCD_EXTERN unsigned int(*eb_aom_variance128x64)(const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse);

    unsigned int eb_aom_variance128x128_c(const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse);
    unsigned int eb_aom_variance128x128_avx2(const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse);
    RTCD_EXTERN unsigned int(*eb_aom_variance128x128)(const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse);

    void eb_aom_ifft16x16_float_avx2(const float *input, float *temp, float *output);
    RTCD_EXTERN void(*eb_aom_ifft16x16_float)(const float *input, float *temp, float *output);

    void eb_aom_ifft2x2_float_c(const float *input, float *temp, float *output);
    RTCD_EXTERN void(*eb_aom_ifft2x2_float)(const float *input, float *temp, float *output);

    void eb_aom_ifft32x32_float_avx2(const float *input, float *temp, float *output);
    RTCD_EXTERN void(*eb_aom_ifft32x32_float)(const float *input, float *temp, float *output);

    void eb_aom_ifft4x4_float_sse2(const float *input, float *temp, float *output);
    RTCD_EXTERN void(*eb_aom_ifft4x4_float)(const float *input, float *temp, float *output);

    void eb_aom_ifft8x8_float_avx2(const float *input, float *temp, float *output);
    RTCD_EXTERN void(*eb_aom_ifft8x8_float)(const float *input, float *temp, float *output);

    void eb_aom_fft16x16_float_c(const float *input, float *temp, float *output);
    void eb_aom_fft16x16_float_avx2(const float *input, float *temp, float *output);
    RTCD_EXTERN void(*eb_aom_fft16x16_float)(const float *input, float *temp, float *output);

    void eb_aom_fft2x2_float_c(const float *input, float *temp, float *output);
    RTCD_EXTERN void(*eb_aom_fft2x2_float)(const float *input, float *temp, float *output);

    void eb_aom_fft32x32_float_c(const float *input, float *temp, float *output);
    void eb_aom_fft32x32_float_avx2(const float *input, float *temp, float *output);
    RTCD_EXTERN void(*eb_aom_fft32x32_float)(const float *input, float *temp, float *output);

    void eb_aom_fft4x4_float_c(const float *input, float *temp, float *output);
    void eb_aom_fft4x4_float_sse2(const float *input, float *temp, float *output);
    RTCD_EXTERN void(*eb_aom_fft4x4_float)(const float *input, float *temp, float *output);

    void eb_aom_fft8x8_float_c(const float *input, float *temp, float *output);
    void eb_aom_fft8x8_float_avx2(const float *input, float *temp, float *output);
    RTCD_EXTERN void(*eb_aom_fft8x8_float)(const float *input, float *temp, float *output);

    void eb_aom_highbd_dc_128_predictor_16x16_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_dc_128_predictor_16x16_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_dc_128_predictor_16x16)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_dc_128_predictor_16x32_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_dc_128_predictor_16x32_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_dc_128_predictor_16x32)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_dc_128_predictor_16x4_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_dc_128_predictor_16x4_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_dc_128_predictor_16x4)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_dc_128_predictor_16x64_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_dc_128_predictor_16x64_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_dc_128_predictor_16x64)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_dc_128_predictor_16x8_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_dc_128_predictor_16x8_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_dc_128_predictor_16x8)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_dc_128_predictor_2x2_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_dc_128_predictor_2x2)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_dc_128_predictor_32x16_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_dc_128_predictor_32x16_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_dc_128_predictor_32x16)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_dc_128_predictor_32x32_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_dc_128_predictor_32x32_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_dc_128_predictor_32x32)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_dc_128_predictor_32x64_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_dc_128_predictor_32x64_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_dc_128_predictor_32x64)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_dc_128_predictor_32x8_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_dc_128_predictor_32x8_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_dc_128_predictor_32x8)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_dc_128_predictor_4x16_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_dc_128_predictor_4x16_sse2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_dc_128_predictor_4x16)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_dc_128_predictor_4x4_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_dc_128_predictor_4x4_sse2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_dc_128_predictor_4x4)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_dc_128_predictor_4x8_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_dc_128_predictor_4x8_sse2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_dc_128_predictor_4x8)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_dc_128_predictor_64x16_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_dc_128_predictor_64x16_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_dc_128_predictor_64x16)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_dc_128_predictor_64x32_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_dc_128_predictor_64x32_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_dc_128_predictor_64x32)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_dc_128_predictor_64x64_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_dc_128_predictor_64x64_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_dc_128_predictor_64x64)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_dc_128_predictor_8x16_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_dc_128_predictor_8x16_sse2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_dc_128_predictor_8x16)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_dc_128_predictor_8x32_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_dc_128_predictor_8x32_sse2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_dc_128_predictor_8x32)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_dc_128_predictor_8x4_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_dc_128_predictor_8x4_sse2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_dc_128_predictor_8x4)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_dc_128_predictor_8x8_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_dc_128_predictor_8x8_sse2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_dc_128_predictor_8x8)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_dc_left_predictor_16x16_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_dc_left_predictor_16x16_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_dc_left_predictor_16x16)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_dc_left_predictor_16x32_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_dc_left_predictor_16x32_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_dc_left_predictor_16x32)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_dc_left_predictor_16x4_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_dc_left_predictor_16x4_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_dc_left_predictor_16x4)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_dc_left_predictor_16x64_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_dc_left_predictor_16x64_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_dc_left_predictor_16x64)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_dc_left_predictor_16x8_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_dc_left_predictor_16x8_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_dc_left_predictor_16x8)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_dc_left_predictor_2x2_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_dc_left_predictor_2x2)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_dc_left_predictor_32x16_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_dc_left_predictor_32x16_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_dc_left_predictor_32x16_avx512(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_dc_left_predictor_32x16)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_dc_left_predictor_32x32_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_dc_left_predictor_32x32_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_dc_left_predictor_32x32_avx512(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_dc_left_predictor_32x32)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_dc_left_predictor_32x64_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_dc_left_predictor_32x64_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_dc_left_predictor_32x64_avx512(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_dc_left_predictor_32x64)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_dc_left_predictor_32x8_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_dc_left_predictor_32x8_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_dc_left_predictor_32x8_avx512(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_dc_left_predictor_32x8)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_dc_left_predictor_4x16_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_dc_left_predictor_4x16_sse2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_dc_left_predictor_4x16)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_dc_left_predictor_4x4_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_dc_left_predictor_4x4_sse2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_dc_left_predictor_4x4)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_dc_left_predictor_4x8_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_dc_left_predictor_4x8_sse2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_dc_left_predictor_4x8)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_dc_left_predictor_64x16_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_dc_left_predictor_64x16_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_dc_left_predictor_64x16_avx512(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_dc_left_predictor_64x16)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_dc_left_predictor_64x32_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_dc_left_predictor_64x32_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_dc_left_predictor_64x32_avx512(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_dc_left_predictor_64x32)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_dc_left_predictor_64x64_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_dc_left_predictor_64x64_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_dc_left_predictor_64x64_avx512(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_dc_left_predictor_64x64)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_dc_left_predictor_8x16_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_dc_left_predictor_8x16_sse2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_dc_left_predictor_8x16)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_dc_left_predictor_8x32_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_dc_left_predictor_8x32_sse2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_dc_left_predictor_8x32)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_dc_left_predictor_8x4_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_dc_left_predictor_8x4_sse2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_dc_left_predictor_8x4)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_dc_left_predictor_8x8_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_dc_left_predictor_8x8_sse2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_dc_left_predictor_8x8)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_dc_predictor_16x16_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_dc_predictor_16x16_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_dc_predictor_16x16)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_dc_predictor_16x32_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_dc_predictor_16x32_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_dc_predictor_16x32)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_dc_predictor_16x4_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_dc_predictor_16x4_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_dc_predictor_16x4)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_dc_predictor_16x64_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_dc_predictor_16x64_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_dc_predictor_16x64)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_dc_predictor_16x8_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_dc_predictor_16x8_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_dc_predictor_16x8)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_dc_predictor_2x2_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_dc_predictor_2x2)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_dc_predictor_32x16_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_dc_predictor_32x16_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_dc_predictor_32x16_avx512(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_dc_predictor_32x16)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_dc_predictor_32x32_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_dc_predictor_32x32_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_dc_predictor_32x32_avx512(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_dc_predictor_32x32)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_dc_predictor_32x64_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_dc_predictor_32x64_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_dc_predictor_32x64_avx512(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_dc_predictor_32x64)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_dc_predictor_32x8_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_dc_predictor_32x8_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_dc_predictor_32x8_avx512(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_dc_predictor_32x8)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_dc_predictor_4x16_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_dc_predictor_4x16_sse2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_dc_predictor_4x16)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_dc_predictor_4x4_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_dc_predictor_4x4_sse2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_dc_predictor_4x4)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_dc_predictor_4x8_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_dc_predictor_4x8_sse2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_dc_predictor_4x8)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_dc_predictor_64x16_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_dc_predictor_64x16_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_dc_predictor_64x16_avx512(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_dc_predictor_64x16)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_dc_predictor_64x32_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_dc_predictor_64x32_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_dc_predictor_64x32_avx512(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_dc_predictor_64x32)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_dc_predictor_64x64_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_dc_predictor_64x64_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_dc_predictor_64x64_avx512(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_dc_predictor_64x64)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_dc_predictor_8x16_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_dc_predictor_8x16_sse2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_dc_predictor_8x16)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_dc_predictor_8x32_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_dc_predictor_8x32_sse2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_dc_predictor_8x32)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_dc_predictor_8x4_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_dc_predictor_8x4_sse2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_dc_predictor_8x4)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_dc_predictor_8x8_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_dc_predictor_8x8_sse2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_dc_predictor_8x8)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_dc_top_predictor_16x16_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_dc_top_predictor_16x16_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_dc_top_predictor_16x16)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_dc_top_predictor_16x32_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_dc_top_predictor_16x32_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_dc_top_predictor_16x32)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_dc_top_predictor_16x4_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_dc_top_predictor_16x4_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_dc_top_predictor_16x4)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_dc_top_predictor_16x64_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_dc_top_predictor_16x64_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_dc_top_predictor_16x64)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_dc_top_predictor_16x8_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_dc_top_predictor_16x8_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_dc_top_predictor_16x8)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_dc_top_predictor_2x2_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_dc_top_predictor_2x2)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_dc_top_predictor_32x16_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_dc_top_predictor_32x16_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_dc_top_predictor_32x16_avx512(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_dc_top_predictor_32x16)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_dc_top_predictor_32x32_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_dc_top_predictor_32x32_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_dc_top_predictor_32x32_avx512(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_dc_top_predictor_32x32)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_dc_top_predictor_32x64_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_dc_top_predictor_32x64_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_dc_top_predictor_32x64_avx512(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_dc_top_predictor_32x64)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_dc_top_predictor_32x8_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_dc_top_predictor_32x8_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_dc_top_predictor_32x8_avx512(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_dc_top_predictor_32x8)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_dc_top_predictor_4x16_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_dc_top_predictor_4x16_sse2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_dc_top_predictor_4x16)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_dc_top_predictor_4x4_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_dc_top_predictor_4x4_sse2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_dc_top_predictor_4x4)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_dc_top_predictor_4x8_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_dc_top_predictor_4x8_sse2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_dc_top_predictor_4x8)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_dc_top_predictor_64x16_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_dc_top_predictor_64x16_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_dc_top_predictor_64x16_avx512(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_dc_top_predictor_64x16)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_dc_top_predictor_64x32_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_dc_top_predictor_64x32_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_dc_top_predictor_64x32_avx512(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_dc_top_predictor_64x32)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_dc_top_predictor_64x64_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_dc_top_predictor_64x64_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_dc_top_predictor_64x64_avx512(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_dc_top_predictor_64x64)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_dc_top_predictor_8x16_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_dc_top_predictor_8x16_sse2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_dc_top_predictor_8x16)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_dc_top_predictor_8x32_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_dc_top_predictor_8x32)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_dc_top_predictor_8x4_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_dc_top_predictor_8x4_sse2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_dc_top_predictor_8x4)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_dc_top_predictor_8x8_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_dc_top_predictor_8x8_sse2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_dc_top_predictor_8x8)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_h_predictor_16x16_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_h_predictor_16x16_sse2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_h_predictor_16x16)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_h_predictor_16x32_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_h_predictor_16x32_sse2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_h_predictor_16x32)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_h_predictor_16x4_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_h_predictor_16x4_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_h_predictor_16x4)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_h_predictor_16x64_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_h_predictor_16x64_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_h_predictor_16x64)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_h_predictor_16x8_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_h_predictor_16x8_sse2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_h_predictor_16x8)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_h_predictor_2x2_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_h_predictor_2x2)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_h_predictor_32x16_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_h_predictor_32x16_sse2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_h_predictor_32x16)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_h_predictor_32x32_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_h_predictor_32x32_sse2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_h_predictor_32x32)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_h_predictor_32x64_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_h_predictor_32x64_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_h_predictor_32x64)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_h_predictor_32x8_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_h_predictor_32x8_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_h_predictor_32x8)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_h_predictor_4x16_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_h_predictor_4x16_sse2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_h_predictor_4x16)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_h_predictor_4x4_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_h_predictor_4x4_sse2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_h_predictor_4x4)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_h_predictor_4x8_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_h_predictor_4x8_sse2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_h_predictor_4x8)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_h_predictor_64x16_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_h_predictor_64x16_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_h_predictor_64x16)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_h_predictor_64x32_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_h_predictor_64x32_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_h_predictor_64x32)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_h_predictor_64x64_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_h_predictor_64x64_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_h_predictor_64x64)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_h_predictor_8x16_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_h_predictor_8x16_sse2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_h_predictor_8x16)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_h_predictor_8x32_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_h_predictor_8x32_sse2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_h_predictor_8x32)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_h_predictor_8x4_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_h_predictor_8x4_sse2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_h_predictor_8x4)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_h_predictor_8x8_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_h_predictor_8x8_sse2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_h_predictor_8x8)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_smooth_h_predictor_16x16_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_smooth_h_predictor_16x16_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_smooth_h_predictor_16x16)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_smooth_h_predictor_16x32_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_smooth_h_predictor_16x32_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_smooth_h_predictor_16x32)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_smooth_h_predictor_16x4_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_smooth_h_predictor_16x4_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_smooth_h_predictor_16x4)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_smooth_h_predictor_16x64_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_smooth_h_predictor_16x64_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_smooth_h_predictor_16x64)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_smooth_h_predictor_16x8_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_smooth_h_predictor_16x8_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_smooth_h_predictor_16x8)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_smooth_h_predictor_2x2_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_smooth_h_predictor_2x2)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_smooth_h_predictor_32x16_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_smooth_h_predictor_32x16_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_smooth_h_predictor_32x16)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_smooth_h_predictor_32x32_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_smooth_h_predictor_32x32_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_smooth_h_predictor_32x32)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_smooth_h_predictor_32x64_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_smooth_h_predictor_32x64_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_smooth_h_predictor_32x64)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_smooth_h_predictor_32x8_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_smooth_h_predictor_32x8_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_smooth_h_predictor_32x8)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_smooth_h_predictor_4x16_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_smooth_h_predictor_4x16_ssse3(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_smooth_h_predictor_4x16)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_smooth_h_predictor_4x4_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_smooth_h_predictor_4x4_ssse3(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_smooth_h_predictor_4x4)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_smooth_h_predictor_4x8_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_smooth_h_predictor_4x8_ssse3(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_smooth_h_predictor_4x8)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_smooth_h_predictor_64x16_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_smooth_h_predictor_64x16_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_smooth_h_predictor_64x16)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_smooth_h_predictor_64x32_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_smooth_h_predictor_64x32_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_smooth_h_predictor_64x32)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_smooth_h_predictor_64x64_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_smooth_h_predictor_64x64_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_smooth_h_predictor_64x64)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_smooth_h_predictor_8x16_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_smooth_h_predictor_8x16_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_smooth_h_predictor_8x16)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_smooth_h_predictor_8x32_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_smooth_h_predictor_8x32_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_smooth_h_predictor_8x32)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_smooth_h_predictor_8x4_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_smooth_h_predictor_8x4_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_smooth_h_predictor_8x4)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_smooth_h_predictor_8x8_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_smooth_h_predictor_8x8_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_smooth_h_predictor_8x8)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_smooth_predictor_16x16_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_smooth_predictor_16x16_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_smooth_predictor_16x16)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_smooth_predictor_16x32_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_smooth_predictor_16x32_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_smooth_predictor_16x32)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_smooth_predictor_16x4_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_smooth_predictor_16x4_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_smooth_predictor_16x4)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_smooth_predictor_16x64_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_smooth_predictor_16x64_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_smooth_predictor_16x64)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_smooth_predictor_16x8_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_smooth_predictor_16x8_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_smooth_predictor_16x8)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_smooth_predictor_2x2_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_smooth_predictor_2x2)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_smooth_predictor_32x16_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_smooth_predictor_32x16_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_smooth_predictor_32x16)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_smooth_predictor_32x32_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_smooth_predictor_32x32_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_smooth_predictor_32x32)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_smooth_predictor_32x64_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_smooth_predictor_32x64_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_smooth_predictor_32x64)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_smooth_predictor_32x8_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_smooth_predictor_32x8_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_smooth_predictor_32x8)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_smooth_predictor_4x16_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_smooth_predictor_4x16_ssse3(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_smooth_predictor_4x16)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_smooth_predictor_4x4_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_smooth_predictor_4x4_ssse3(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_smooth_predictor_4x4)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_smooth_predictor_4x8_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_smooth_predictor_4x8_ssse3(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_smooth_predictor_4x8)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_smooth_predictor_64x16_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_smooth_predictor_64x16_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_smooth_predictor_64x16)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_smooth_predictor_64x32_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_smooth_predictor_64x32_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_smooth_predictor_64x32)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_smooth_predictor_64x64_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_smooth_predictor_64x64_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_smooth_predictor_64x64)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_smooth_predictor_8x16_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_smooth_predictor_8x16_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_smooth_predictor_8x16)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_smooth_predictor_8x32_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_smooth_predictor_8x32_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_smooth_predictor_8x32)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_smooth_predictor_8x4_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_smooth_predictor_8x4_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_smooth_predictor_8x4)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_smooth_predictor_8x8_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_smooth_predictor_8x8_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_smooth_predictor_8x8)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_smooth_v_predictor_16x16_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_smooth_v_predictor_16x16_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_smooth_v_predictor_16x16)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_smooth_v_predictor_16x32_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_smooth_v_predictor_16x32_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_smooth_v_predictor_16x32)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_smooth_v_predictor_16x4_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_smooth_v_predictor_16x4_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_smooth_v_predictor_16x4)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_smooth_v_predictor_16x64_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_smooth_v_predictor_16x64_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_smooth_v_predictor_16x64)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_smooth_v_predictor_16x8_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_smooth_v_predictor_16x8_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_smooth_v_predictor_16x8)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_smooth_v_predictor_2x2_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_smooth_v_predictor_2x2)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_smooth_v_predictor_32x16_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_smooth_v_predictor_32x16_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_smooth_v_predictor_32x16)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_smooth_v_predictor_32x32_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_smooth_v_predictor_32x32_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_smooth_v_predictor_32x32)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_smooth_v_predictor_32x64_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_smooth_v_predictor_32x64_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_smooth_v_predictor_32x64)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_smooth_v_predictor_32x8_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_smooth_v_predictor_32x8_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_smooth_v_predictor_32x8)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_smooth_v_predictor_4x16_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_smooth_v_predictor_4x16_ssse3(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_smooth_v_predictor_4x16)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_smooth_v_predictor_4x4_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_smooth_v_predictor_4x4_ssse3(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_smooth_v_predictor_4x4)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_smooth_v_predictor_4x8_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_smooth_v_predictor_4x8_ssse3(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_smooth_v_predictor_4x8)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_smooth_v_predictor_64x16_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_smooth_v_predictor_64x16_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_smooth_v_predictor_64x16)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_smooth_v_predictor_64x32_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_smooth_v_predictor_64x32_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_smooth_v_predictor_64x32)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_smooth_v_predictor_64x64_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_smooth_v_predictor_64x64_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_smooth_v_predictor_64x64)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_smooth_v_predictor_8x16_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_smooth_v_predictor_8x16_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_smooth_v_predictor_8x16)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_smooth_v_predictor_8x32_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_smooth_v_predictor_8x32_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_smooth_v_predictor_8x32)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_smooth_v_predictor_8x4_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_smooth_v_predictor_8x4_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_smooth_v_predictor_8x4)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_smooth_v_predictor_8x8_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_smooth_v_predictor_8x8_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_smooth_v_predictor_8x8)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_v_predictor_16x16_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_v_predictor_16x16_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_v_predictor_16x16)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_v_predictor_16x32_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_v_predictor_16x32_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_v_predictor_16x32)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_v_predictor_16x4_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_v_predictor_16x4_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_v_predictor_16x4)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_v_predictor_16x64_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_v_predictor_16x64_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_v_predictor_16x64)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_v_predictor_16x8_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_v_predictor_16x8_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_v_predictor_16x8)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_v_predictor_2x2_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_v_predictor_2x2)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_v_predictor_32x16_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_v_predictor_32x16_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_v_predictor_32x16)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_v_predictor_32x32_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_v_predictor_32x32_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_v_predictor_32x32)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_v_predictor_32x64_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_v_predictor_32x64_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_v_predictor_32x64)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_v_predictor_32x8_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_v_predictor_32x8_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_v_predictor_32x8)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_v_predictor_4x16_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_v_predictor_4x16_sse2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_v_predictor_4x16)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_v_predictor_4x4_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_v_predictor_4x4_sse2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_v_predictor_4x4)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_v_predictor_4x8_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_v_predictor_4x8_sse2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_v_predictor_4x8)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_v_predictor_64x16_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_v_predictor_64x16_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_v_predictor_64x16)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_v_predictor_64x32_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_v_predictor_64x32_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_v_predictor_64x32)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_v_predictor_64x64_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_v_predictor_64x64_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_v_predictor_64x64)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_v_predictor_8x16_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_v_predictor_8x16_sse2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_v_predictor_8x16)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_v_predictor_8x32_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_v_predictor_8x32_sse2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_v_predictor_8x32)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_v_predictor_8x4_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_v_predictor_8x4_sse2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_v_predictor_8x4)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_aom_highbd_v_predictor_8x8_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void eb_aom_highbd_v_predictor_8x8_sse2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*eb_aom_highbd_v_predictor_8x8)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void eb_av1_dr_prediction_z1_c(uint8_t *dst, ptrdiff_t stride, int32_t bw, int32_t bh, const uint8_t *above, const uint8_t *left, int32_t upsample_above, int32_t dx, int32_t dy);
    void eb_av1_dr_prediction_z1_avx2(uint8_t *dst, ptrdiff_t stride, int32_t bw, int32_t bh, const uint8_t *above, const uint8_t *left, int32_t upsample_above, int32_t dx, int32_t dy);
    RTCD_EXTERN void(*eb_av1_dr_prediction_z1)(uint8_t *dst, ptrdiff_t stride, int32_t bw, int32_t bh, const uint8_t *above, const uint8_t *left, int32_t upsample_above, int32_t dx, int32_t dy);

    void eb_av1_dr_prediction_z2_c(uint8_t *dst, ptrdiff_t stride, int32_t bw, int32_t bh, const uint8_t *above, const uint8_t *left, int32_t upsample_above, int32_t upsample_left, int32_t dx, int32_t dy);
    void eb_av1_dr_prediction_z2_avx2(uint8_t *dst, ptrdiff_t stride, int32_t bw, int32_t bh, const uint8_t *above, const uint8_t *left, int32_t upsample_above, int32_t upsample_left, int32_t dx, int32_t dy);
    RTCD_EXTERN void(*eb_av1_dr_prediction_z2)(uint8_t *dst, ptrdiff_t stride, int32_t bw, int32_t bh, const uint8_t *above, const uint8_t *left, int32_t upsample_above, int32_t upsample_left, int32_t dx, int32_t dy);

    void eb_av1_dr_prediction_z3_c(uint8_t *dst, ptrdiff_t stride, int32_t bw, int32_t bh, const uint8_t *above, const uint8_t *left, int32_t upsample_left, int32_t dx, int32_t dy);
    void eb_av1_dr_prediction_z3_avx2(uint8_t *dst, ptrdiff_t stride, int32_t bw, int32_t bh, const uint8_t *above, const uint8_t *left, int32_t upsample_left, int32_t dx, int32_t dy);
    RTCD_EXTERN void(*eb_av1_dr_prediction_z3)(uint8_t *dst, ptrdiff_t stride, int32_t bw, int32_t bh, const uint8_t *above, const uint8_t *left, int32_t upsample_left, int32_t dx, int32_t dy);

    void eb_av1_highbd_dr_prediction_z1_c(uint16_t *dst, ptrdiff_t stride, int32_t bw, int32_t bh, const uint16_t *above, const uint16_t *left, int32_t upsample_above, int32_t dx, int32_t dy, int32_t bd);
    void eb_av1_highbd_dr_prediction_z1_avx2(uint16_t *dst, ptrdiff_t stride, int32_t bw, int32_t bh, const uint16_t *above, const uint16_t *left, int32_t upsample_above, int32_t dx, int32_t dy, int32_t bd);
    RTCD_EXTERN void(*eb_av1_highbd_dr_prediction_z1)(uint16_t *dst, ptrdiff_t stride, int32_t bw, int32_t bh, const uint16_t *above, const uint16_t *left, int32_t upsample_above, int32_t dx, int32_t dy, int32_t bd);

    void eb_av1_highbd_dr_prediction_z2_c(uint16_t *dst, ptrdiff_t stride, int32_t bw, int32_t bh, const uint16_t *above, const uint16_t *left, int32_t upsample_above, int32_t upsample_left, int32_t dx, int32_t dy, int32_t bd);
    void eb_av1_highbd_dr_prediction_z2_avx2(uint16_t *dst, ptrdiff_t stride, int32_t bw, int32_t bh, const uint16_t *above, const uint16_t *left, int32_t upsample_above, int32_t upsample_left, int32_t dx, int32_t dy, int32_t bd);
    RTCD_EXTERN void(*eb_av1_highbd_dr_prediction_z2)(uint16_t *dst, ptrdiff_t stride, int32_t bw, int32_t bh, const uint16_t *above, const uint16_t *left, int32_t upsample_above, int32_t upsample_left, int32_t dx, int32_t dy, int32_t bd);

    void eb_av1_highbd_dr_prediction_z3_c(uint16_t *dst, ptrdiff_t stride, int32_t bw, int32_t bh, const uint16_t *above, const uint16_t *left, int32_t upsample_left, int32_t dx, int32_t dy, int32_t bd);
    void eb_av1_highbd_dr_prediction_z3_avx2(uint16_t *dst, ptrdiff_t stride, int32_t bw, int32_t bh, const uint16_t *above, const uint16_t *left, int32_t upsample_left, int32_t dx, int32_t dy, int32_t bd);
    RTCD_EXTERN void(*eb_av1_highbd_dr_prediction_z3)(uint16_t *dst, ptrdiff_t stride, int32_t bw, int32_t bh, const uint16_t *above, const uint16_t *left, int32_t upsample_left, int32_t dx, int32_t dy, int32_t bd);

    void eb_av1_filter_intra_predictor_c(uint8_t *dst, ptrdiff_t stride, TxSize tx_size, const uint8_t *above, const uint8_t *left, int32_t mode);
    RTCD_EXTERN void (*eb_av1_filter_intra_predictor) (uint8_t *dst, ptrdiff_t stride, TxSize tx_size, const uint8_t *above, const uint8_t *left, int32_t mode);

    //void eb_av1_get_nz_map_contexts_c(const uint8_t *const levels, const int16_t *const scan, const uint16_t eob, const TxSize tx_size, const TxClass tx_class, int8_t *const coeff_contexts);
    void eb_av1_get_nz_map_contexts_sse2(const uint8_t *const levels, const int16_t *const scan, const uint16_t eob, const TxSize tx_size, const TxClass tx_class, int8_t *const coeff_contexts);
    RTCD_EXTERN void(*eb_av1_get_nz_map_contexts)(const uint8_t *const levels, const int16_t *const scan, const uint16_t eob, const TxSize tx_size, const TxClass tx_class, int8_t *const coeff_contexts);

    void highbd_variance64_c(const uint8_t *a8, int32_t a_stride, const uint8_t *b8, int32_t b_stride, int32_t w, int32_t h, uint64_t *sse);
    void highbd_variance64_avx2(const uint8_t *a8, int32_t a_stride, const uint8_t *b8, int32_t b_stride, int32_t w, int32_t h, uint64_t *sse);
    RTCD_EXTERN void (*highbd_variance64)(const uint8_t *a8, int32_t a_stride, const uint8_t *b8, int32_t b_stride, int32_t w, int32_t h, uint64_t *sse);

    void residual_kernel_c(uint8_t *input, uint32_t input_stride, uint8_t *pred, uint32_t pred_stride, int16_t *residual, uint32_t residual_stride, uint32_t area_width, uint32_t area_height);
    void ResidualKernel_avx2(uint8_t *input, uint32_t input_stride, uint8_t *pred, uint32_t pred_stride, int16_t *residual, uint32_t residual_stride, uint32_t area_width, uint32_t area_height);
    RTCD_EXTERN void(*ResidualKernel)(uint8_t *input, uint32_t input_stride, uint8_t *pred, uint32_t pred_stride, int16_t *residual, uint32_t residual_stride, uint32_t area_width, uint32_t area_height);

    void eb_av1_txb_init_levels_c(const TranLow *const coeff, const int32_t width, const int32_t height, uint8_t *const levels);
    void eb_av1_txb_init_levels_avx2(const TranLow *const coeff, const int32_t width, const int32_t height, uint8_t *const levels);
    RTCD_EXTERN void(*eb_av1_txb_init_levels)(const TranLow *const coeff, const int32_t width, const int32_t height, uint8_t *const levels);

    void eb_aom_dsp_rtcd(void);

#ifdef RTCD_C

    static void setup_rtcd_internal(EbAsm asm_type)
    {
        int32_t flags = HAS_MMX | HAS_SSE | HAS_SSE2 | HAS_SSE3 | HAS_SSSE3 | HAS_SSE4_1 | HAS_SSE4_2 | HAS_AVX;

        if (asm_type == ASM_AVX2)
            flags |= HAS_AVX2;
        //if (asm_type == ASM_NON_AVX2)
        //    flags = ~HAS_AVX2;

        //to use C: flags=0

        picture_average = picture_average_kernel_sse2_intrin;
        if (flags & HAS_SSE2) picture_average = picture_average_kernel_sse2_intrin;

        picture_average1_line = picture_average_kernel1_line_sse2_intrin;
        if (flags & HAS_SSE2) picture_average1_line = picture_average_kernel1_line_sse2_intrin;

        nxm_sad_loop_kernel = sad_loop_kernel_sse4_1_intrin;
        if (flags & HAS_AVX2) nxm_sad_loop_kernel = sad_loop_kernel_avx2_intrin;

        nxm_sad_loop_kernel_sparse = sad_loop_kernel_sparse_sse4_1_intrin;
        if (flags & HAS_AVX2) nxm_sad_loop_kernel_sparse = sad_loop_kernel_sparse_avx2_intrin;

        get_eight_horizontal_search_point_results_8x8_16x16 = get_eight_horizontal_search_point_results_8x8_16x16_pu_sse41_intrin;
        if (flags & HAS_AVX2)get_eight_horizontal_search_point_results_8x8_16x16 = get_eight_horizontal_search_point_results_8x8_16x16_pu_avx2_intrin;

        get_eight_horizontal_search_point_results_32x32_64x64 = get_eight_horizontal_search_point_results_32x32_64x64_pu_sse41_intrin;
        if (flags & HAS_AVX2) get_eight_horizontal_search_point_results_32x32_64x64 = get_eight_horizontal_search_point_results_32x32_64x64_pu_avx2_intrin;

        combined_averaging_ssd = combined_averaging_ssd_c;
        if (flags & HAS_AVX2) combined_averaging_ssd = combined_averaging_ssd_avx2;

        ebav1_smooth_h_predictor = ebav1_smooth_h_predictor_c;
        ebav1_smooth_v_predictor = ebav1_smooth_v_predictor_c;

        noise_extract_luma_weak = noise_extract_luma_weak_c;
        if (flags & HAS_AVX2) noise_extract_luma_weak = noise_extract_luma_weak_avx2_intrin;

        noise_extract_luma_weak_lcu = noise_extract_luma_weak_lcu_c;
        if (flags & HAS_AVX2) noise_extract_luma_weak_lcu = noise_extract_luma_weak_lcu_avx2_intrin;

        strong_luma_filter = noise_extract_luma_strong_c;
        if (flags & HAS_AVX2) strong_luma_filter = noise_extract_luma_strong_avx2_intrin;

        noise_extract_chroma_strong = noise_extract_chroma_strong_c;
        if (flags & HAS_AVX2) noise_extract_chroma_strong = noise_extract_chroma_strong_avx2_intrin;

        noise_extract_chroma_weak = noise_extract_chroma_weak_c;
        if (flags & HAS_AVX2) noise_extract_chroma_weak = noise_extract_chroma_weak_avx2_intrin;

        sum_residual = sum_residual_c;
        if (flags & HAS_AVX2) sum_residual = sum_residual8bit_avx2_intrin;

        memset16bit_block = memset16bit_block_c;
        if (flags & HAS_AVX2) memset16bit_block = memset16bit_block_avx2_intrin;

        apply_temp_filtering_32x32 = apply_temp_filtering_32x32_c;
        if (flags & HAS_SSE4_1) apply_temp_filtering_32x32 = apply_temp_filtering_32x32_sse4_1;

        compute4x4_SAD = compute4x4SAD_Kernel_c;
        if (flags & HAS_AVX2) compute4x4_SAD = compute4x_m_sad_avx2_intrin;

        full_distortion_kernel32_bits = full_distortion_kernel32_bits_c;
        if (flags & HAS_AVX2) full_distortion_kernel32_bits = full_distortion_kernel32_bits_avx2;

        full_distortion_kernel_cbf_zero32_bits = full_distortion_kernel_cbf_zero32_bits_c;
        if (flags & HAS_AVX2) full_distortion_kernel_cbf_zero32_bits = full_distortion_kernel_cbf_zero32_bits_avx2;

        spatial_full_distortion_kernel = spatial_full_distortion_kernel_c;
        if (flags & HAS_AVX2) spatial_full_distortion_kernel = spatial_full_distortion_kernel_avx2;
#ifndef NON_AVX512_SUPPORT
        spatial_full_distortion_kernel = spatial_full_distortion_kernel_avx512;
#endif

        compute8x4_sad_kernel = compute8x4_sad_kernel_c;
        
        compute8x8_satd_u8 = compute8x8_satd_u8_sse4;

        picture_addition_kernel16bit = picture_addition_kernel16bit_sse2_intrin;

        avc_style_luma_interpolation_filter = avc_style_luma_interpolation_filter_ssse3_helper;

        eb_apply_selfguided_restoration = eb_apply_selfguided_restoration_c;
        if (flags & HAS_AVX2) eb_apply_selfguided_restoration = eb_apply_selfguided_restoration_avx2;

        eb_av1_wiener_convolve_add_src = eb_av1_wiener_convolve_add_src_c;
        if (flags & HAS_AVX2) eb_av1_wiener_convolve_add_src = eb_av1_wiener_convolve_add_src_avx2;

        eb_av1_highbd_wiener_convolve_add_src = eb_av1_highbd_wiener_convolve_add_src_c;
        if (flags & HAS_AVX2) eb_av1_highbd_wiener_convolve_add_src = eb_av1_highbd_wiener_convolve_add_src_avx2;

        eb_av1_selfguided_restoration = eb_av1_selfguided_restoration_c;
        if (flags & HAS_AVX2) eb_av1_selfguided_restoration = eb_av1_selfguided_restoration_avx2;

        eb_cdef_find_dir = eb_cdef_find_dir_c;
        if (flags & HAS_AVX2) eb_cdef_find_dir = eb_cdef_find_dir_avx2;

        eb_cdef_filter_block = eb_cdef_filter_block_c;
        if (flags & HAS_AVX2) eb_cdef_filter_block = eb_cdef_filter_block_avx2;
        eb_compute_cdef_dist = compute_cdef_dist_c;
        if (flags & HAS_AVX2) eb_compute_cdef_dist = compute_cdef_dist_avx2;

        eb_copy_rect8_8bit_to_16bit = eb_copy_rect8_8bit_to_16bit_c;
        if (flags & HAS_AVX2) eb_copy_rect8_8bit_to_16bit = eb_copy_rect8_8bit_to_16bit_avx2;

        eb_av1_compute_stats = eb_av1_compute_stats_c;
        if (flags & HAS_AVX2) eb_av1_compute_stats = eb_av1_compute_stats_avx2;
        eb_av1_compute_stats_highbd = eb_av1_compute_stats_highbd_c;
        if (flags & HAS_AVX2) eb_av1_compute_stats_highbd = eb_av1_compute_stats_highbd_avx2;
        eb_av1_lowbd_pixel_proj_error = eb_av1_lowbd_pixel_proj_error_c;
        if (flags & HAS_AVX2) eb_av1_lowbd_pixel_proj_error = eb_av1_lowbd_pixel_proj_error_avx2;
        eb_av1_highbd_pixel_proj_error = eb_av1_highbd_pixel_proj_error_c;
        if (flags & HAS_AVX2) eb_av1_highbd_pixel_proj_error = eb_av1_highbd_pixel_proj_error_avx2;
        eb_av1_filter_intra_edge_high = eb_av1_filter_intra_edge_high_c;
        if (flags & HAS_SSE4_1) eb_av1_filter_intra_edge_high = eb_av1_filter_intra_edge_high_sse4_1;
        eb_av1_calc_frame_error = eb_av1_calc_frame_error_c;
        if (flags & HAS_AVX2) eb_av1_calc_frame_error = eb_av1_calc_frame_error_avx2;
        eb_av1_highbd_convolve_2d_copy_sr = eb_av1_highbd_convolve_2d_copy_sr_c;
        if (flags & HAS_AVX2) eb_av1_highbd_convolve_2d_copy_sr = eb_av1_highbd_convolve_2d_copy_sr_avx2;
        eb_av1_highbd_jnt_convolve_2d_copy = eb_av1_highbd_jnt_convolve_2d_copy_c;
        if (flags & HAS_AVX2) eb_av1_highbd_jnt_convolve_2d_copy = eb_av1_highbd_jnt_convolve_2d_copy_avx2;
        eb_av1_highbd_convolve_y_sr = eb_av1_highbd_convolve_y_sr_c;
        if (flags & HAS_AVX2) eb_av1_highbd_convolve_y_sr = eb_av1_highbd_convolve_y_sr_avx2;
        eb_av1_highbd_convolve_2d_sr = eb_av1_highbd_convolve_2d_sr_c;
        if (flags & HAS_AVX2) eb_av1_highbd_convolve_2d_sr = eb_av1_highbd_convolve_2d_sr_avx2;
        eb_av1_highbd_jnt_convolve_2d = eb_av1_highbd_jnt_convolve_2d_c;
        if (flags & HAS_AVX2) eb_av1_highbd_jnt_convolve_2d = eb_av1_highbd_jnt_convolve_2d_avx2;
        eb_av1_highbd_jnt_convolve_x = eb_av1_highbd_jnt_convolve_x_c;
        if (flags & HAS_AVX2) eb_av1_highbd_jnt_convolve_x = eb_av1_highbd_jnt_convolve_x_avx2;
        eb_av1_highbd_jnt_convolve_y = eb_av1_highbd_jnt_convolve_y_c;
        if (flags & HAS_AVX2) eb_av1_highbd_jnt_convolve_y = eb_av1_highbd_jnt_convolve_y_avx2;
        eb_av1_highbd_convolve_x_sr = eb_av1_highbd_convolve_x_sr_c;
        if (flags & HAS_AVX2) eb_av1_highbd_convolve_x_sr = eb_av1_highbd_convolve_x_sr_avx2;
        eb_subtract_average = eb_subtract_average_c;
        if (flags & HAS_AVX2) eb_subtract_average = eb_subtract_average_avx2;

        eb_av1_filter_intra_edge = eb_av1_filter_intra_edge_high_c_old;

        if (flags & HAS_SSE4_1) eb_av1_filter_intra_edge = eb_av1_filter_intra_edge_sse4_1;

        eb_smooth_v_predictor = smooth_v_predictor_c;
        if (flags & HAS_SSSE3) eb_smooth_v_predictor = eb_smooth_v_predictor_all_ssse3;

        eb_smooth_h_predictor = smooth_h_predictor_c;
        if (flags & HAS_SSSE3) eb_smooth_h_predictor = eb_smooth_h_predictor_all_ssse3;

        get_proj_subspace = get_proj_subspace_c;
        if (flags & HAS_AVX2) get_proj_subspace = get_proj_subspace_avx2;

        search_one_dual = search_one_dual_c;
        if (flags & HAS_AVX2) search_one_dual = search_one_dual_avx2;

        eb_aom_mse16x16 = eb_aom_mse16x16_c;
        if (flags & HAS_AVX2) eb_aom_mse16x16 = eb_aom_mse16x16_avx2;

        eb_av1_convolve_2d_copy_sr = eb_av1_convolve_2d_copy_sr_c;
        if (flags & HAS_AVX2) eb_av1_convolve_2d_copy_sr = eb_av1_convolve_2d_copy_sr_avx2;

        eb_av1_convolve_2d_sr = eb_av1_convolve_2d_sr_c;
        if (flags & HAS_AVX2) eb_av1_convolve_2d_sr = eb_av1_convolve_2d_sr_avx2;

        eb_av1_jnt_convolve_2d_copy = eb_av1_jnt_convolve_2d_copy_c;
        if (flags & HAS_AVX2) eb_av1_jnt_convolve_2d_copy = eb_av1_jnt_convolve_2d_copy_avx2;

        eb_av1_convolve_x_sr = eb_av1_convolve_x_sr_c;
        if (flags & HAS_AVX2) eb_av1_convolve_x_sr = eb_av1_convolve_x_sr_avx2;
        eb_av1_convolve_y_sr = eb_av1_convolve_y_sr_c;
        if (flags & HAS_AVX2) eb_av1_convolve_y_sr = eb_av1_convolve_y_sr_avx2;

        eb_av1_jnt_convolve_x = eb_av1_jnt_convolve_x_c;
        if (flags & HAS_AVX2) eb_av1_jnt_convolve_x = eb_av1_jnt_convolve_x_avx2;
        eb_av1_jnt_convolve_y = eb_av1_jnt_convolve_y_c;
        if (flags & HAS_AVX2) eb_av1_jnt_convolve_y = eb_av1_jnt_convolve_y_avx2;

        eb_av1_jnt_convolve_2d = eb_av1_jnt_convolve_2d_c;
        if (flags & HAS_AVX2) eb_av1_jnt_convolve_2d = eb_av1_jnt_convolve_2d_avx2;

        eb_aom_quantize_b = eb_aom_quantize_b_c_II;
        if (flags & HAS_AVX2) eb_aom_quantize_b = eb_aom_highbd_quantize_b_avx2;

        eb_aom_quantize_b_32x32 = eb_aom_quantize_b_32x32_c_II;
        if (flags & HAS_AVX2) eb_aom_quantize_b_32x32 = eb_aom_highbd_quantize_b_32x32_avx2;

        eb_aom_highbd_quantize_b_32x32 = eb_aom_highbd_quantize_b_32x32_c;
        if (flags & HAS_AVX2) eb_aom_highbd_quantize_b_32x32 = eb_aom_highbd_quantize_b_32x32_avx2;

        eb_aom_highbd_quantize_b = eb_aom_highbd_quantize_b_c;
        if (flags & HAS_AVX2) eb_aom_highbd_quantize_b = eb_aom_highbd_quantize_b_avx2;

        eb_av1_inv_txfm2d_add_16x16 = eb_av1_inv_txfm2d_add_16x16_c;
        if (flags & HAS_AVX2) eb_av1_inv_txfm2d_add_16x16 = eb_av1_inv_txfm2d_add_16x16_avx2;
        eb_av1_inv_txfm2d_add_32x32 = eb_av1_inv_txfm2d_add_32x32_c;
        if (flags & HAS_AVX2) eb_av1_inv_txfm2d_add_32x32 = eb_av1_inv_txfm2d_add_32x32_avx2;
        eb_av1_inv_txfm2d_add_4x4 = eb_av1_inv_txfm2d_add_4x4_c;
        if (flags & HAS_AVX2) eb_av1_inv_txfm2d_add_4x4 = eb_av1_inv_txfm2d_add_4x4_avx2;
        eb_av1_inv_txfm2d_add_64x64 = eb_av1_inv_txfm2d_add_64x64_c;
        if (flags & HAS_SSE4_1) eb_av1_inv_txfm2d_add_64x64 = eb_av1_inv_txfm2d_add_64x64_sse4_1;
        eb_av1_inv_txfm2d_add_8x8 = eb_av1_inv_txfm2d_add_8x8_c;
        if (flags & HAS_AVX2) eb_av1_inv_txfm2d_add_8x8 = eb_av1_inv_txfm2d_add_8x8_avx2;

        eb_av1_inv_txfm2d_add_8x16 = eb_av1_inv_txfm2d_add_8x16_c;
        if (flags & HAS_AVX2) eb_av1_inv_txfm2d_add_8x16 = eb_av1_highbd_inv_txfm_add_avx2;
        eb_av1_inv_txfm2d_add_16x8 = eb_av1_inv_txfm2d_add_16x8_c;
        if (flags & HAS_AVX2) eb_av1_inv_txfm2d_add_16x8 = eb_av1_highbd_inv_txfm_add_avx2;
        eb_av1_inv_txfm2d_add_16x32 = eb_av1_inv_txfm2d_add_16x32_c;
        if (flags & HAS_AVX2) eb_av1_inv_txfm2d_add_16x32 = eb_av1_highbd_inv_txfm_add_avx2;
        eb_av1_inv_txfm2d_add_32x16 = eb_av1_inv_txfm2d_add_32x16_c;
        if (flags & HAS_AVX2) eb_av1_inv_txfm2d_add_32x16 = eb_av1_highbd_inv_txfm_add_avx2;
        eb_av1_inv_txfm2d_add_32x8 = eb_av1_inv_txfm2d_add_32x8_c;
        if (flags & HAS_AVX2) eb_av1_inv_txfm2d_add_32x8 = eb_av1_highbd_inv_txfm_add_avx2;
        eb_av1_inv_txfm2d_add_8x32 = eb_av1_inv_txfm2d_add_8x32_c;
        if (flags & HAS_AVX2) eb_av1_inv_txfm2d_add_8x32 = eb_av1_highbd_inv_txfm_add_avx2;
        eb_av1_inv_txfm2d_add_32x64 = eb_av1_inv_txfm2d_add_32x64_c;
        if (flags & HAS_AVX2) eb_av1_inv_txfm2d_add_32x64 = eb_av1_highbd_inv_txfm_add_avx2;
        eb_av1_inv_txfm2d_add_64x32 = eb_av1_inv_txfm2d_add_64x32_c;
        if (flags & HAS_AVX2) eb_av1_inv_txfm2d_add_64x32 = eb_av1_highbd_inv_txfm_add_avx2;
        eb_av1_inv_txfm2d_add_16x64 = eb_av1_inv_txfm2d_add_16x64_c;
        if (flags & HAS_AVX2) eb_av1_inv_txfm2d_add_16x64 = eb_av1_highbd_inv_txfm_add_avx2;
        eb_av1_inv_txfm2d_add_64x16 = eb_av1_inv_txfm2d_add_64x16_c;
        if (flags & HAS_AVX2) eb_av1_inv_txfm2d_add_64x16 = eb_av1_highbd_inv_txfm_add_avx2;
        eb_av1_inv_txfm2d_add_4x8 = eb_av1_inv_txfm2d_add_4x8_c;
        if (flags & HAS_SSE4_1) eb_av1_inv_txfm2d_add_4x8 = eb_av1_inv_txfm2d_add_4x8_sse4_1;
        eb_av1_inv_txfm2d_add_8x4 = eb_av1_inv_txfm2d_add_8x4_c;
        if (flags & HAS_SSE4_1) eb_av1_inv_txfm2d_add_8x4 = eb_av1_inv_txfm2d_add_8x4_sse4_1;
        eb_av1_inv_txfm2d_add_4x16 = eb_av1_inv_txfm2d_add_4x16_c;
        if (flags & HAS_SSE4_1) eb_av1_inv_txfm2d_add_4x16 = eb_av1_inv_txfm2d_add_4x16_sse4_1;
        eb_av1_inv_txfm2d_add_16x4 = eb_av1_inv_txfm2d_add_16x4_c;
        if (flags & HAS_SSE4_1) eb_av1_inv_txfm2d_add_16x4 = eb_av1_inv_txfm2d_add_16x4_sse4_1;

        eb_av1_inv_txfm_add = eb_av1_inv_txfm_add_c;
        if (flags & HAS_SSSE3) eb_av1_inv_txfm_add = eb_av1_inv_txfm_add_ssse3;
        if (flags & HAS_AVX2) eb_av1_inv_txfm_add = eb_av1_inv_txfm_add_avx2;

        eb_av1_quantize_fp = eb_av1_quantize_fp_c;
        if (flags & HAS_AVX2) eb_av1_quantize_fp = eb_av1_quantize_fp_avx2;

        eb_av1_quantize_fp_32x32 = eb_av1_quantize_fp_32x32_c;
        if (flags & HAS_AVX2) eb_av1_quantize_fp_32x32 = eb_av1_quantize_fp_32x32_avx2;

        eb_av1_quantize_fp_64x64 = eb_av1_quantize_fp_64x64_c;
        if (flags & HAS_AVX2) eb_av1_quantize_fp_64x64 = eb_av1_quantize_fp_64x64_avx2;

        highbd_variance64 = highbd_variance64_c;
        if (flags & HAS_AVX2) highbd_variance64 = highbd_variance64_avx2;

        //aom_highbd_8_mse16x16 = eb_aom_highbd_8_mse16x16_c;
        if (flags & HAS_SSE2) eb_aom_highbd_8_mse16x16 = eb_aom_highbd_8_mse16x16_sse2;

        eb_av1_upsample_intra_edge = eb_av1_upsample_intra_edge_c;
        if (flags & HAS_SSE4_1) eb_av1_upsample_intra_edge = eb_av1_upsample_intra_edge_sse4_1;

        eb_av1_warp_affine = eb_av1_warp_affine_c;
        if (flags & HAS_AVX2) eb_av1_warp_affine = eb_av1_warp_affine_avx2;

        eb_av1_filter_intra_predictor = eb_av1_filter_intra_predictor_c;

        eb_aom_highbd_smooth_v_predictor_16x16 = eb_aom_highbd_smooth_v_predictor_16x16_c;
        if (flags & HAS_AVX2) eb_aom_highbd_smooth_v_predictor_16x16 = eb_aom_highbd_smooth_v_predictor_16x16_avx2;
        eb_aom_highbd_smooth_v_predictor_16x32 = eb_aom_highbd_smooth_v_predictor_16x32_c;
        if (flags & HAS_AVX2) eb_aom_highbd_smooth_v_predictor_16x32 = eb_aom_highbd_smooth_v_predictor_16x32_avx2;
        eb_aom_highbd_smooth_v_predictor_16x4 = eb_aom_highbd_smooth_v_predictor_16x4_c;
        if (flags & HAS_AVX2) eb_aom_highbd_smooth_v_predictor_16x4 = eb_aom_highbd_smooth_v_predictor_16x4_avx2;
        eb_aom_highbd_smooth_v_predictor_16x64 = eb_aom_highbd_smooth_v_predictor_16x64_c;
        if (flags & HAS_AVX2) eb_aom_highbd_smooth_v_predictor_16x64 = eb_aom_highbd_smooth_v_predictor_16x64_avx2;
        eb_aom_highbd_smooth_v_predictor_16x8 = eb_aom_highbd_smooth_v_predictor_16x8_c;
        if (flags & HAS_AVX2) eb_aom_highbd_smooth_v_predictor_16x8 = eb_aom_highbd_smooth_v_predictor_16x8_avx2;
        eb_aom_highbd_smooth_v_predictor_2x2 = eb_aom_highbd_smooth_v_predictor_2x2_c;
        eb_aom_highbd_smooth_v_predictor_32x16 = eb_aom_highbd_smooth_v_predictor_32x16_c;
        if (flags & HAS_AVX2) eb_aom_highbd_smooth_v_predictor_32x16 = eb_aom_highbd_smooth_v_predictor_32x16_avx2;
        eb_aom_highbd_smooth_v_predictor_32x32 = eb_aom_highbd_smooth_v_predictor_32x32_c;
        if (flags & HAS_AVX2) eb_aom_highbd_smooth_v_predictor_32x32 = eb_aom_highbd_smooth_v_predictor_32x32_avx2;
        eb_aom_highbd_smooth_v_predictor_32x64 = eb_aom_highbd_smooth_v_predictor_32x64_c;
        if (flags & HAS_AVX2) eb_aom_highbd_smooth_v_predictor_32x64 = eb_aom_highbd_smooth_v_predictor_32x64_avx2;
        eb_aom_highbd_smooth_v_predictor_32x8 = eb_aom_highbd_smooth_v_predictor_32x8_c;
        if (flags & HAS_AVX2) eb_aom_highbd_smooth_v_predictor_32x8 = eb_aom_highbd_smooth_v_predictor_32x8_avx2;
        eb_aom_highbd_smooth_v_predictor_4x16 = eb_aom_highbd_smooth_v_predictor_4x16_c;
        if (flags & HAS_SSSE3) eb_aom_highbd_smooth_v_predictor_4x16 = eb_aom_highbd_smooth_v_predictor_4x16_ssse3;
        eb_aom_highbd_smooth_v_predictor_4x4 = eb_aom_highbd_smooth_v_predictor_4x4_c;
        if (flags & HAS_SSSE3) eb_aom_highbd_smooth_v_predictor_4x4 = eb_aom_highbd_smooth_v_predictor_4x4_ssse3;
        eb_aom_highbd_smooth_v_predictor_4x8 = eb_aom_highbd_smooth_v_predictor_4x8_c;
        if (flags & HAS_SSSE3) eb_aom_highbd_smooth_v_predictor_4x8 = eb_aom_highbd_smooth_v_predictor_4x8_ssse3;
        eb_aom_highbd_smooth_v_predictor_64x16 = eb_aom_highbd_smooth_v_predictor_64x16_c;
        if (flags & HAS_AVX2) eb_aom_highbd_smooth_v_predictor_64x16 = eb_aom_highbd_smooth_v_predictor_64x16_avx2;
        eb_aom_highbd_smooth_v_predictor_64x32 = eb_aom_highbd_smooth_v_predictor_64x32_c;
        if (flags & HAS_AVX2) eb_aom_highbd_smooth_v_predictor_64x32 = eb_aom_highbd_smooth_v_predictor_64x32_avx2;
        eb_aom_highbd_smooth_v_predictor_64x64 = eb_aom_highbd_smooth_v_predictor_64x64_c;
        if (flags & HAS_AVX2) eb_aom_highbd_smooth_v_predictor_64x64 = eb_aom_highbd_smooth_v_predictor_64x64_avx2;
        eb_aom_highbd_smooth_v_predictor_8x16 = eb_aom_highbd_smooth_v_predictor_8x16_c;
        if (flags & HAS_AVX2) eb_aom_highbd_smooth_v_predictor_8x16 = eb_aom_highbd_smooth_v_predictor_8x16_avx2;
        eb_aom_highbd_smooth_v_predictor_8x32 = eb_aom_highbd_smooth_v_predictor_8x32_c;
        if (flags & HAS_AVX2) eb_aom_highbd_smooth_v_predictor_8x32 = eb_aom_highbd_smooth_v_predictor_8x32_avx2;
        eb_aom_highbd_smooth_v_predictor_8x4 = eb_aom_highbd_smooth_v_predictor_8x4_c;
        if (flags & HAS_AVX2) eb_aom_highbd_smooth_v_predictor_8x4 = eb_aom_highbd_smooth_v_predictor_8x4_avx2;
        eb_aom_highbd_smooth_v_predictor_8x8 = eb_aom_highbd_smooth_v_predictor_8x8_c;
        if (flags & HAS_AVX2) eb_aom_highbd_smooth_v_predictor_8x8 = eb_aom_highbd_smooth_v_predictor_8x8_avx2;

        eb_cfl_predict_lbd = eb_cfl_predict_lbd_c;
        if (flags & HAS_AVX2) eb_cfl_predict_lbd = eb_cfl_predict_lbd_avx2;
        eb_cfl_predict_hbd = eb_cfl_predict_hbd_c;
        if (flags & HAS_AVX2) eb_cfl_predict_hbd = eb_cfl_predict_hbd_avx2;

        eb_av1_dr_prediction_z1 = eb_av1_dr_prediction_z1_c;
        if (flags & HAS_AVX2) eb_av1_dr_prediction_z1 = eb_av1_dr_prediction_z1_avx2;
        eb_av1_dr_prediction_z2 = eb_av1_dr_prediction_z2_c;
        if (flags & HAS_AVX2) eb_av1_dr_prediction_z2 = eb_av1_dr_prediction_z2_avx2;
        eb_av1_dr_prediction_z3 = eb_av1_dr_prediction_z3_c;
        if (flags & HAS_AVX2) eb_av1_dr_prediction_z3 = eb_av1_dr_prediction_z3_avx2;
        eb_av1_highbd_dr_prediction_z1 = eb_av1_highbd_dr_prediction_z1_c;
        if (flags & HAS_AVX2) eb_av1_highbd_dr_prediction_z1 = eb_av1_highbd_dr_prediction_z1_avx2;
        eb_av1_highbd_dr_prediction_z2 = eb_av1_highbd_dr_prediction_z2_c;
        if (flags & HAS_AVX2) eb_av1_highbd_dr_prediction_z2 = eb_av1_highbd_dr_prediction_z2_avx2;
        eb_av1_highbd_dr_prediction_z3 = eb_av1_highbd_dr_prediction_z3_c;
        if (flags & HAS_AVX2) eb_av1_highbd_dr_prediction_z3 = eb_av1_highbd_dr_prediction_z3_avx2;
        //av1_get_nz_map_contexts = eb_av1_get_nz_map_contexts_c;
        if (flags & HAS_SSE2) eb_av1_get_nz_map_contexts = eb_av1_get_nz_map_contexts_sse2;

        ResidualKernel = residual_kernel_c;
        if (flags & HAS_AVX2) ResidualKernel = ResidualKernel_avx2;

        eb_av1_txb_init_levels = eb_av1_txb_init_levels_c;
        if (flags & HAS_AVX2) eb_av1_txb_init_levels = eb_av1_txb_init_levels_avx2;
    eb_aom_paeth_predictor_16x16 = eb_aom_paeth_predictor_16x16_c;
    if (flags & HAS_SSSE3) eb_aom_paeth_predictor_16x16 = eb_aom_paeth_predictor_16x16_ssse3;
    if (flags & HAS_AVX2) eb_aom_paeth_predictor_16x16 = eb_aom_paeth_predictor_16x16_avx2;
    eb_aom_paeth_predictor_16x32 = eb_aom_paeth_predictor_16x32_c;
    if (flags & HAS_SSSE3) eb_aom_paeth_predictor_16x32 = eb_aom_paeth_predictor_16x32_ssse3;
    if (flags & HAS_AVX2) eb_aom_paeth_predictor_16x32 = eb_aom_paeth_predictor_16x32_avx2;
    eb_aom_paeth_predictor_16x4 = eb_aom_paeth_predictor_16x4_c;
    if (flags & HAS_SSSE3) eb_aom_paeth_predictor_16x4 = eb_aom_paeth_predictor_16x4_ssse3;
    eb_aom_paeth_predictor_16x64 = eb_aom_paeth_predictor_16x64_c;
    if (flags & HAS_SSSE3) eb_aom_paeth_predictor_16x64 = eb_aom_paeth_predictor_16x64_ssse3;
    if (flags & HAS_AVX2) eb_aom_paeth_predictor_16x64 = eb_aom_paeth_predictor_16x64_avx2;
    eb_aom_paeth_predictor_16x8 = eb_aom_paeth_predictor_16x8_c;
    if (flags & HAS_SSSE3) eb_aom_paeth_predictor_16x8 = eb_aom_paeth_predictor_16x8_ssse3;
    if (flags & HAS_AVX2) eb_aom_paeth_predictor_16x8 = eb_aom_paeth_predictor_16x8_avx2;
    eb_aom_paeth_predictor_32x16 = eb_aom_paeth_predictor_32x16_c;
    if (flags & HAS_SSSE3) eb_aom_paeth_predictor_32x16 = eb_aom_paeth_predictor_32x16_ssse3;
    if (flags & HAS_AVX2) eb_aom_paeth_predictor_32x16 = eb_aom_paeth_predictor_32x16_avx2;
    eb_aom_paeth_predictor_32x32 = eb_aom_paeth_predictor_32x32_c;
    if (flags & HAS_SSSE3) eb_aom_paeth_predictor_32x32 = eb_aom_paeth_predictor_32x32_ssse3;
    if (flags & HAS_AVX2) eb_aom_paeth_predictor_32x32 = eb_aom_paeth_predictor_32x32_avx2;
    eb_aom_paeth_predictor_32x64 = eb_aom_paeth_predictor_32x64_c;
    if (flags & HAS_SSSE3) eb_aom_paeth_predictor_32x64 = eb_aom_paeth_predictor_32x64_ssse3;
    if (flags & HAS_AVX2) eb_aom_paeth_predictor_32x64 = eb_aom_paeth_predictor_32x64_avx2;
    eb_aom_paeth_predictor_32x8 = eb_aom_paeth_predictor_32x8_c;
    if (flags & HAS_SSSE3) eb_aom_paeth_predictor_32x8 = eb_aom_paeth_predictor_32x8_ssse3;
    eb_aom_paeth_predictor_4x16 = eb_aom_paeth_predictor_4x16_c;
    if (flags & HAS_SSSE3) eb_aom_paeth_predictor_4x16 = eb_aom_paeth_predictor_4x16_ssse3;
    eb_aom_paeth_predictor_4x4 = eb_aom_paeth_predictor_4x4_c;
    if (flags & HAS_SSSE3) eb_aom_paeth_predictor_4x4 = eb_aom_paeth_predictor_4x4_ssse3;
    eb_aom_paeth_predictor_4x8 = eb_aom_paeth_predictor_4x8_c;
    if (flags & HAS_SSSE3) eb_aom_paeth_predictor_4x8 = eb_aom_paeth_predictor_4x8_ssse3;
    eb_aom_paeth_predictor_64x16 = eb_aom_paeth_predictor_64x16_c;
    if (flags & HAS_SSSE3) eb_aom_paeth_predictor_64x16 = eb_aom_paeth_predictor_64x16_ssse3;
    if (flags & HAS_AVX2) eb_aom_paeth_predictor_64x16 = eb_aom_paeth_predictor_64x16_avx2;
    eb_aom_paeth_predictor_64x32 = eb_aom_paeth_predictor_64x32_c;
    if (flags & HAS_SSSE3) eb_aom_paeth_predictor_64x32 = eb_aom_paeth_predictor_64x32_ssse3;
    if (flags & HAS_AVX2) eb_aom_paeth_predictor_64x32 = eb_aom_paeth_predictor_64x32_avx2;
    eb_aom_paeth_predictor_64x64 = eb_aom_paeth_predictor_64x64_c;
    if (flags & HAS_SSSE3) eb_aom_paeth_predictor_64x64 = eb_aom_paeth_predictor_64x64_ssse3;
    if (flags & HAS_AVX2) eb_aom_paeth_predictor_64x64 = eb_aom_paeth_predictor_64x64_avx2;
    eb_aom_paeth_predictor_8x16 = eb_aom_paeth_predictor_8x16_c;
    if (flags & HAS_SSSE3) eb_aom_paeth_predictor_8x16 = eb_aom_paeth_predictor_8x16_ssse3;
    eb_aom_paeth_predictor_8x32 = eb_aom_paeth_predictor_8x32_c;
    if (flags & HAS_SSSE3) eb_aom_paeth_predictor_8x32 = eb_aom_paeth_predictor_8x32_ssse3;
    eb_aom_paeth_predictor_8x4 = eb_aom_paeth_predictor_8x4_c;
    if (flags & HAS_SSSE3) eb_aom_paeth_predictor_8x4 = eb_aom_paeth_predictor_8x4_ssse3;
    eb_aom_paeth_predictor_8x8 = eb_aom_paeth_predictor_8x8_c;
    if (flags & HAS_SSSE3) eb_aom_paeth_predictor_8x8 = eb_aom_paeth_predictor_8x8_ssse3;

        eb_aom_highbd_paeth_predictor_16x16 = eb_aom_highbd_paeth_predictor_16x16_c;
        if (flags & HAS_AVX2) eb_aom_highbd_paeth_predictor_16x16 = eb_aom_highbd_paeth_predictor_16x16_avx2;
        eb_aom_highbd_paeth_predictor_16x32 = eb_aom_highbd_paeth_predictor_16x32_c;
        if (flags & HAS_AVX2) eb_aom_highbd_paeth_predictor_16x32 = eb_aom_highbd_paeth_predictor_16x32_avx2;
        eb_aom_highbd_paeth_predictor_16x4 = eb_aom_highbd_paeth_predictor_16x4_c;
        if (flags & HAS_AVX2) eb_aom_highbd_paeth_predictor_16x4 = eb_aom_highbd_paeth_predictor_16x4_avx2;
        eb_aom_highbd_paeth_predictor_16x64 = eb_aom_highbd_paeth_predictor_16x64_c;
        if (flags & HAS_AVX2) eb_aom_highbd_paeth_predictor_16x64 = eb_aom_highbd_paeth_predictor_16x64_avx2;
        eb_aom_highbd_paeth_predictor_16x8 = eb_aom_highbd_paeth_predictor_16x8_c;
        if (flags & HAS_AVX2) eb_aom_highbd_paeth_predictor_16x8 = eb_aom_highbd_paeth_predictor_16x8_avx2;
        eb_aom_highbd_paeth_predictor_2x2 = eb_aom_highbd_paeth_predictor_2x2_c;
        if (flags & HAS_AVX2) eb_aom_highbd_paeth_predictor_2x2 = eb_aom_highbd_paeth_predictor_2x2_avx2;
        eb_aom_highbd_paeth_predictor_32x16 = eb_aom_highbd_paeth_predictor_32x16_c;
        if (flags & HAS_AVX2) eb_aom_highbd_paeth_predictor_32x16 = eb_aom_highbd_paeth_predictor_32x16_avx2;
        eb_aom_highbd_paeth_predictor_32x32 = eb_aom_highbd_paeth_predictor_32x32_c;
        if (flags & HAS_AVX2) eb_aom_highbd_paeth_predictor_32x32 = eb_aom_highbd_paeth_predictor_32x32_avx2;
        eb_aom_highbd_paeth_predictor_32x64 = eb_aom_highbd_paeth_predictor_32x64_c;
        if (flags & HAS_AVX2) eb_aom_highbd_paeth_predictor_32x64 = eb_aom_highbd_paeth_predictor_32x64_avx2;
        eb_aom_highbd_paeth_predictor_32x8 = eb_aom_highbd_paeth_predictor_32x8_c;
        if (flags & HAS_AVX2) eb_aom_highbd_paeth_predictor_32x8 = eb_aom_highbd_paeth_predictor_32x8_avx2;
        eb_aom_highbd_paeth_predictor_4x16 = eb_aom_highbd_paeth_predictor_4x16_c;
        if (flags & HAS_AVX2) eb_aom_highbd_paeth_predictor_4x16 = eb_aom_highbd_paeth_predictor_4x16_avx2;
        eb_aom_highbd_paeth_predictor_4x4 = eb_aom_highbd_paeth_predictor_4x4_c;
        if (flags & HAS_AVX2) eb_aom_highbd_paeth_predictor_4x4 = eb_aom_highbd_paeth_predictor_4x4_avx2;
        eb_aom_highbd_paeth_predictor_4x8 = eb_aom_highbd_paeth_predictor_4x8_c;
        if (flags & HAS_AVX2) eb_aom_highbd_paeth_predictor_4x8 = eb_aom_highbd_paeth_predictor_4x8_avx2;
        eb_aom_highbd_paeth_predictor_64x16 = eb_aom_highbd_paeth_predictor_64x16_c;
        if (flags & HAS_AVX2) eb_aom_highbd_paeth_predictor_64x16 = eb_aom_highbd_paeth_predictor_64x16_avx2;
        eb_aom_highbd_paeth_predictor_64x32 = eb_aom_highbd_paeth_predictor_64x32_c;
        if (flags & HAS_AVX2) eb_aom_highbd_paeth_predictor_64x32 = eb_aom_highbd_paeth_predictor_64x32_avx2;
        eb_aom_highbd_paeth_predictor_64x64 = eb_aom_highbd_paeth_predictor_64x64_c;
        if (flags & HAS_AVX2) eb_aom_highbd_paeth_predictor_64x64 = eb_aom_highbd_paeth_predictor_64x64_avx2;
        eb_aom_highbd_paeth_predictor_8x16 = eb_aom_highbd_paeth_predictor_8x16_c;
        if (flags & HAS_AVX2) eb_aom_highbd_paeth_predictor_8x16 = eb_aom_highbd_paeth_predictor_8x16_avx2;
        eb_aom_highbd_paeth_predictor_8x32 = eb_aom_highbd_paeth_predictor_8x32_c;
        if (flags & HAS_AVX2) eb_aom_highbd_paeth_predictor_8x32 = eb_aom_highbd_paeth_predictor_8x32_avx2;
        eb_aom_highbd_paeth_predictor_8x4 = eb_aom_highbd_paeth_predictor_8x4_c;
        if (flags & HAS_AVX2) eb_aom_highbd_paeth_predictor_8x4 = eb_aom_highbd_paeth_predictor_8x4_avx2;
        eb_aom_highbd_paeth_predictor_8x8 = eb_aom_highbd_paeth_predictor_8x8_c;
        if (flags & HAS_AVX2) eb_aom_highbd_paeth_predictor_8x8 = eb_aom_highbd_paeth_predictor_8x8_avx2;

        eb_aom_dc_predictor_4x4 = eb_aom_dc_predictor_4x4_c;
        if (flags & HAS_SSE2) eb_aom_dc_predictor_4x4 = eb_aom_dc_predictor_4x4_sse2;
        eb_aom_dc_predictor_8x8 = eb_aom_dc_predictor_8x8_c;
        if (flags & HAS_SSE2) eb_aom_dc_predictor_8x8 = eb_aom_dc_predictor_8x8_sse2;
        eb_aom_dc_predictor_16x16 = eb_aom_dc_predictor_16x16_c;
        if (flags & HAS_SSE2) eb_aom_dc_predictor_16x16 = eb_aom_dc_predictor_16x16_sse2;
        eb_aom_dc_predictor_32x32 = eb_aom_dc_predictor_32x32_c;
        if (flags & HAS_AVX2) eb_aom_dc_predictor_32x32 = eb_aom_dc_predictor_32x32_avx2;
        eb_aom_dc_predictor_64x64 = eb_aom_dc_predictor_64x64_c;
        if (flags & HAS_AVX2) eb_aom_dc_predictor_64x64 = eb_aom_dc_predictor_64x64_avx2;
        eb_aom_dc_predictor_32x16 = eb_aom_dc_predictor_32x16_c;
        if (flags & HAS_AVX2) eb_aom_dc_predictor_32x16 = eb_aom_dc_predictor_32x16_avx2;
        eb_aom_dc_predictor_32x64 = eb_aom_dc_predictor_32x64_c;
        if (flags & HAS_AVX2) eb_aom_dc_predictor_32x64 = eb_aom_dc_predictor_32x64_avx2;
        eb_aom_dc_predictor_64x16 = eb_aom_dc_predictor_64x16_c;
        if (flags & HAS_AVX2) eb_aom_dc_predictor_64x16 = eb_aom_dc_predictor_64x16_avx2;
        eb_aom_dc_predictor_8x16 = eb_aom_dc_predictor_8x16_c;
        if (flags & HAS_SSE2) eb_aom_dc_predictor_8x16 = eb_aom_dc_predictor_8x16_sse2;
        eb_aom_dc_predictor_8x32 = eb_aom_dc_predictor_8x32_c;
        if (flags & HAS_SSE2) eb_aom_dc_predictor_8x32 = eb_aom_dc_predictor_8x32_sse2;
        eb_aom_dc_predictor_8x4 = eb_aom_dc_predictor_8x4_c;
        if (flags & HAS_SSE2) eb_aom_dc_predictor_8x4 = eb_aom_dc_predictor_8x4_sse2;
        eb_aom_dc_predictor_64x32 = eb_aom_dc_predictor_64x32_c;
        if (flags & HAS_AVX2) eb_aom_dc_predictor_64x32 = eb_aom_dc_predictor_64x32_avx2;
        eb_aom_dc_predictor_16x32 = eb_aom_dc_predictor_16x32_c;
        if (flags & HAS_SSE2) eb_aom_dc_predictor_16x32 = eb_aom_dc_predictor_16x32_sse2;
        eb_aom_dc_predictor_16x4 = eb_aom_dc_predictor_16x4_c;
        if (flags & HAS_SSE2) eb_aom_dc_predictor_16x4 = eb_aom_dc_predictor_16x4_sse2;
        eb_aom_dc_predictor_16x64 = eb_aom_dc_predictor_16x64_c;
        if (flags & HAS_SSE2) eb_aom_dc_predictor_16x64 = eb_aom_dc_predictor_16x64_sse2;
        eb_aom_dc_predictor_16x8 = eb_aom_dc_predictor_16x8_c;
        if (flags & HAS_SSE2) eb_aom_dc_predictor_16x8 = eb_aom_dc_predictor_16x8_sse2;
        eb_aom_dc_predictor_32x8 = eb_aom_dc_predictor_32x8_c;
        if (flags & HAS_SSE2) eb_aom_dc_predictor_32x8 = eb_aom_dc_predictor_32x8_sse2;
        eb_aom_dc_predictor_4x16 = eb_aom_dc_predictor_4x16_c;
        if (flags & HAS_SSE2) eb_aom_dc_predictor_4x16 = eb_aom_dc_predictor_4x16_sse2;
        eb_aom_dc_predictor_4x8 = eb_aom_dc_predictor_4x8_c;
        if (flags & HAS_SSE2) eb_aom_dc_predictor_4x8 = eb_aom_dc_predictor_4x8_sse2;

        eb_aom_dc_top_predictor_4x4 = eb_aom_dc_top_predictor_4x4_c;
        if (flags & HAS_SSE2) eb_aom_dc_top_predictor_4x4 = eb_aom_dc_top_predictor_4x4_sse2;
        eb_aom_dc_top_predictor_8x8 = eb_aom_dc_top_predictor_8x8_c;
        if (flags & HAS_SSE2) eb_aom_dc_top_predictor_8x8 = eb_aom_dc_top_predictor_8x8_sse2;
        eb_aom_dc_top_predictor_16x16 = eb_aom_dc_top_predictor_16x16_c;
        if (flags & HAS_SSE2) eb_aom_dc_top_predictor_16x16 = eb_aom_dc_top_predictor_16x16_sse2;
        eb_aom_dc_top_predictor_32x32 = eb_aom_dc_top_predictor_32x32_c;
        if (flags & HAS_AVX2) eb_aom_dc_top_predictor_32x32 = eb_aom_dc_top_predictor_32x32_avx2;
        eb_aom_dc_top_predictor_64x64 = eb_aom_dc_top_predictor_64x64_c;
        if (flags & HAS_AVX2) eb_aom_dc_top_predictor_64x64 = eb_aom_dc_top_predictor_64x64_avx2;
        eb_aom_dc_top_predictor_16x32 = eb_aom_dc_top_predictor_16x32_c;
        if (flags & HAS_SSE2) eb_aom_dc_top_predictor_16x32 = eb_aom_dc_top_predictor_16x32_sse2;
        eb_aom_dc_top_predictor_16x4 = eb_aom_dc_top_predictor_16x4_c;
        if (flags & HAS_SSE2) eb_aom_dc_top_predictor_16x4 = eb_aom_dc_top_predictor_16x4_sse2;
        eb_aom_dc_top_predictor_16x64 = eb_aom_dc_top_predictor_16x64_c;
        if (flags & HAS_SSE2) eb_aom_dc_top_predictor_16x64 = eb_aom_dc_top_predictor_16x64_sse2;
        eb_aom_dc_top_predictor_16x8 = eb_aom_dc_top_predictor_16x8_c;
        if (flags & HAS_SSE2) eb_aom_dc_top_predictor_16x8 = eb_aom_dc_top_predictor_16x8_sse2;
        eb_aom_dc_top_predictor_32x16 = eb_aom_dc_top_predictor_32x16_c;
        if (flags & HAS_AVX2) eb_aom_dc_top_predictor_32x16 = eb_aom_dc_top_predictor_32x16_avx2;
        eb_aom_dc_top_predictor_32x64 = eb_aom_dc_top_predictor_32x64_c;
        if (flags & HAS_AVX2) eb_aom_dc_top_predictor_32x64 = eb_aom_dc_top_predictor_32x64_avx2;
        eb_aom_dc_top_predictor_32x8 = eb_aom_dc_top_predictor_32x8_c;
        if (flags & HAS_SSE2) eb_aom_dc_top_predictor_32x8 = eb_aom_dc_top_predictor_32x8_sse2;
        eb_aom_dc_top_predictor_4x16 = eb_aom_dc_top_predictor_4x16_c;
        if (flags & HAS_SSE2) eb_aom_dc_top_predictor_4x16 = eb_aom_dc_top_predictor_4x16_sse2;
        eb_aom_dc_top_predictor_4x8 = eb_aom_dc_top_predictor_4x8_c;
        if (flags & HAS_SSE2) eb_aom_dc_top_predictor_4x8 = eb_aom_dc_top_predictor_4x8_sse2;
        eb_aom_dc_top_predictor_64x16 = eb_aom_dc_top_predictor_64x16_c;
        if (flags & HAS_AVX2) eb_aom_dc_top_predictor_64x16 = eb_aom_dc_top_predictor_64x16_avx2;
        eb_aom_dc_top_predictor_64x32 = eb_aom_dc_top_predictor_64x32_c;
        if (flags & HAS_AVX2) eb_aom_dc_top_predictor_64x32 = eb_aom_dc_top_predictor_64x32_avx2;
        eb_aom_dc_top_predictor_8x16 = eb_aom_dc_top_predictor_8x16_c;
        if (flags & HAS_SSE2) eb_aom_dc_top_predictor_8x16 = eb_aom_dc_top_predictor_8x16_sse2;
        eb_aom_dc_top_predictor_8x32 = eb_aom_dc_top_predictor_8x32_c;
        if (flags & HAS_SSE2) eb_aom_dc_top_predictor_8x32 = eb_aom_dc_top_predictor_8x32_sse2;
        eb_aom_dc_top_predictor_8x4 = eb_aom_dc_top_predictor_8x4_c;
        if (flags & HAS_SSE2) eb_aom_dc_top_predictor_8x4 = eb_aom_dc_top_predictor_8x4_sse2;

        eb_aom_dc_left_predictor_4x4 = eb_aom_dc_left_predictor_4x4_c;
        if (flags & HAS_SSE2) eb_aom_dc_left_predictor_4x4 = eb_aom_dc_left_predictor_4x4_sse2;
        eb_aom_dc_left_predictor_8x8 = eb_aom_dc_left_predictor_8x8_c;
        if (flags & HAS_SSE2) eb_aom_dc_left_predictor_8x8 = eb_aom_dc_left_predictor_8x8_sse2;
        eb_aom_dc_left_predictor_16x16 = eb_aom_dc_left_predictor_16x16_c;
        if (flags & HAS_SSE2) eb_aom_dc_left_predictor_16x16 = eb_aom_dc_left_predictor_16x16_sse2;
        eb_aom_dc_left_predictor_32x32 = eb_aom_dc_left_predictor_32x32_c;
        if (flags & HAS_AVX2) eb_aom_dc_left_predictor_32x32 = eb_aom_dc_left_predictor_32x32_avx2;
        eb_aom_dc_left_predictor_64x64 = eb_aom_dc_left_predictor_64x64_c;
        if (flags & HAS_AVX2) eb_aom_dc_left_predictor_64x64 = eb_aom_dc_left_predictor_64x64_avx2;
        eb_aom_dc_left_predictor_16x32 = eb_aom_dc_left_predictor_16x32_c;
        if (flags & HAS_SSE2) eb_aom_dc_left_predictor_16x32 = eb_aom_dc_left_predictor_16x32_sse2;
        eb_aom_dc_left_predictor_16x4 = eb_aom_dc_left_predictor_16x4_c;
        if (flags & HAS_SSE2) eb_aom_dc_left_predictor_16x4 = eb_aom_dc_left_predictor_16x4_sse2;
        eb_aom_dc_left_predictor_16x64 = eb_aom_dc_left_predictor_16x64_c;
        if (flags & HAS_SSE2) eb_aom_dc_left_predictor_16x64 = eb_aom_dc_left_predictor_16x64_sse2;
        eb_aom_dc_left_predictor_16x8 = eb_aom_dc_left_predictor_16x8_c;
        if (flags & HAS_SSE2)  eb_aom_dc_left_predictor_16x8 = eb_aom_dc_left_predictor_16x8_sse2;
        eb_aom_dc_left_predictor_32x16 = eb_aom_dc_left_predictor_32x16_c;
        if (flags & HAS_AVX2) eb_aom_dc_left_predictor_32x16 = eb_aom_dc_left_predictor_32x16_avx2;
        eb_aom_dc_left_predictor_32x64 = eb_aom_dc_left_predictor_32x64_c;
        if (flags & HAS_AVX2) eb_aom_dc_left_predictor_32x64 = eb_aom_dc_left_predictor_32x64_avx2;
        eb_aom_dc_left_predictor_64x16 = eb_aom_dc_left_predictor_64x16_c;
        if (flags & HAS_AVX2) eb_aom_dc_left_predictor_64x16 = eb_aom_dc_left_predictor_64x16_avx2;
        eb_aom_dc_left_predictor_64x32 = eb_aom_dc_left_predictor_64x32_c;
        if (flags & HAS_AVX2) eb_aom_dc_left_predictor_64x32 = eb_aom_dc_left_predictor_64x32_avx2;
        eb_aom_dc_left_predictor_32x8 = eb_aom_dc_left_predictor_32x8_c;
        if (flags & HAS_SSE2) eb_aom_dc_left_predictor_32x8 = eb_aom_dc_left_predictor_32x8_sse2;
        eb_aom_dc_left_predictor_4x16 = eb_aom_dc_left_predictor_4x16_c;
        if (flags & HAS_SSE2) eb_aom_dc_left_predictor_4x16 = eb_aom_dc_left_predictor_4x16_sse2;
        eb_aom_dc_left_predictor_4x8 = eb_aom_dc_left_predictor_4x8_c;
        if (flags & HAS_SSE2) eb_aom_dc_left_predictor_4x8 = eb_aom_dc_left_predictor_4x8_sse2;
        eb_aom_dc_left_predictor_8x16 = eb_aom_dc_left_predictor_8x16_c;
        if (flags & HAS_SSE2) eb_aom_dc_left_predictor_8x16 = eb_aom_dc_left_predictor_8x16_sse2;
        eb_aom_dc_left_predictor_8x32 = eb_aom_dc_left_predictor_8x32_c;
        if (flags & HAS_SSE2) eb_aom_dc_left_predictor_8x32 = eb_aom_dc_left_predictor_8x32_sse2;
        eb_aom_dc_left_predictor_8x4 = eb_aom_dc_left_predictor_8x4_c;
        if (flags & HAS_SSE2) eb_aom_dc_left_predictor_8x4 = eb_aom_dc_left_predictor_8x4_sse2;

        eb_aom_dc_128_predictor_4x4 = eb_aom_dc_128_predictor_4x4_c;
        if (flags & HAS_SSE2) eb_aom_dc_128_predictor_4x4 = eb_aom_dc_128_predictor_4x4_sse2;
        eb_aom_dc_128_predictor_8x8 = eb_aom_dc_128_predictor_8x8_c;
        if (flags & HAS_SSE2) eb_aom_dc_128_predictor_8x8 = eb_aom_dc_128_predictor_8x8_sse2;
        eb_aom_dc_128_predictor_16x16 = eb_aom_dc_128_predictor_16x16_c;
        if (flags & HAS_SSE2) eb_aom_dc_128_predictor_16x16 = eb_aom_dc_128_predictor_16x16_sse2;
        eb_aom_dc_128_predictor_32x32 = eb_aom_dc_128_predictor_32x32_c;
        if (flags & HAS_AVX2) eb_aom_dc_128_predictor_32x32 = eb_aom_dc_128_predictor_32x32_avx2;
        eb_aom_dc_128_predictor_64x64 = eb_aom_dc_128_predictor_64x64_c;
        if (flags & HAS_AVX2) eb_aom_dc_128_predictor_64x64 = eb_aom_dc_128_predictor_64x64_avx2;
        eb_aom_dc_128_predictor_16x32 = eb_aom_dc_128_predictor_16x32_c;
        if (flags & HAS_SSE2) eb_aom_dc_128_predictor_16x32 = eb_aom_dc_128_predictor_16x32_sse2;
        eb_aom_dc_128_predictor_16x4 = eb_aom_dc_128_predictor_16x4_c;
        if (flags & HAS_SSE2) eb_aom_dc_128_predictor_16x4 = eb_aom_dc_128_predictor_16x4_sse2;
        eb_aom_dc_128_predictor_16x64 = eb_aom_dc_128_predictor_16x64_c;
        if (flags & HAS_SSE2) eb_aom_dc_128_predictor_16x64 = eb_aom_dc_128_predictor_16x64_sse2;
        eb_aom_dc_128_predictor_16x8 = eb_aom_dc_128_predictor_16x8_c;
        if (flags & HAS_SSE2) eb_aom_dc_128_predictor_16x8 = eb_aom_dc_128_predictor_16x8_sse2;
        eb_aom_dc_128_predictor_32x16 = eb_aom_dc_128_predictor_32x16_c;
        if (flags & HAS_AVX2) eb_aom_dc_128_predictor_32x16 = eb_aom_dc_128_predictor_32x16_avx2;
        eb_aom_dc_128_predictor_32x64 = eb_aom_dc_128_predictor_32x64_c;
        if (flags & HAS_AVX2) eb_aom_dc_128_predictor_32x64 = eb_aom_dc_128_predictor_32x64_avx2;
        eb_aom_dc_128_predictor_32x8 = eb_aom_dc_128_predictor_32x8_c;
        if (flags & HAS_SSE2) eb_aom_dc_128_predictor_32x8 = eb_aom_dc_128_predictor_32x8_sse2;
        eb_aom_dc_128_predictor_4x16 = eb_aom_dc_128_predictor_4x16_c;
        if (flags & HAS_SSE2) eb_aom_dc_128_predictor_4x16 = eb_aom_dc_128_predictor_4x16_sse2;
        eb_aom_dc_128_predictor_4x8 = eb_aom_dc_128_predictor_4x8_c;
        if (flags & HAS_SSE2) eb_aom_dc_128_predictor_4x8 = eb_aom_dc_128_predictor_4x8_sse2;
        eb_aom_dc_128_predictor_64x16 = eb_aom_dc_128_predictor_64x16_c;
        if (flags & HAS_AVX2) eb_aom_dc_128_predictor_64x16 = eb_aom_dc_128_predictor_64x16_avx2;
        eb_aom_dc_128_predictor_64x32 = eb_aom_dc_128_predictor_64x32_c;
        if (flags & HAS_AVX2) eb_aom_dc_128_predictor_64x32 = eb_aom_dc_128_predictor_64x32_avx2;
        eb_aom_dc_128_predictor_8x16 = eb_aom_dc_128_predictor_8x16_c;
        if (flags & HAS_SSE2) eb_aom_dc_128_predictor_8x16 = eb_aom_dc_128_predictor_8x16_sse2;
        eb_aom_dc_128_predictor_8x32 = eb_aom_dc_128_predictor_8x32_c;
        if (flags & HAS_SSE2) eb_aom_dc_128_predictor_8x32 = eb_aom_dc_128_predictor_8x32_sse2;
        eb_aom_dc_128_predictor_8x4 = eb_aom_dc_128_predictor_8x4_c;
        if (flags & HAS_SSE2) eb_aom_dc_128_predictor_8x4 = eb_aom_dc_128_predictor_8x4_sse2;

        eb_aom_smooth_h_predictor_16x32 = eb_aom_smooth_h_predictor_16x32_c;
        if (flags & HAS_SSSE3) eb_aom_smooth_h_predictor_16x32 = eb_aom_smooth_h_predictor_16x32_ssse3;
        eb_aom_smooth_h_predictor_16x4 = eb_aom_smooth_h_predictor_16x4_c;
        if (flags & HAS_SSSE3) eb_aom_smooth_h_predictor_16x4 = eb_aom_smooth_h_predictor_16x4_ssse3;
        eb_aom_smooth_h_predictor_16x64 = eb_aom_smooth_h_predictor_16x64_c;
        if (flags & HAS_SSSE3) eb_aom_smooth_h_predictor_16x64 = eb_aom_smooth_h_predictor_16x64_ssse3;
        eb_aom_smooth_h_predictor_16x8 = eb_aom_smooth_h_predictor_16x8_c;
        if (flags & HAS_SSSE3) eb_aom_smooth_h_predictor_16x8 = eb_aom_smooth_h_predictor_16x8_ssse3;
        eb_aom_smooth_h_predictor_32x16 = eb_aom_smooth_h_predictor_32x16_c;
        if (flags & HAS_SSSE3) eb_aom_smooth_h_predictor_32x16 = eb_aom_smooth_h_predictor_32x16_ssse3;
        eb_aom_smooth_h_predictor_32x64 = eb_aom_smooth_h_predictor_32x64_c;
        if (flags & HAS_SSSE3) eb_aom_smooth_h_predictor_32x64 = eb_aom_smooth_h_predictor_32x64_ssse3;
        eb_aom_smooth_h_predictor_32x8 = eb_aom_smooth_h_predictor_32x8_c;
        if (flags & HAS_SSSE3) eb_aom_smooth_h_predictor_32x8 = eb_aom_smooth_h_predictor_32x8_ssse3;
        eb_aom_smooth_h_predictor_4x16 = eb_aom_smooth_h_predictor_4x16_c;
        if (flags & HAS_SSSE3) eb_aom_smooth_h_predictor_4x16 = eb_aom_smooth_h_predictor_4x16_ssse3;
        eb_aom_smooth_h_predictor_4x8 = eb_aom_smooth_h_predictor_4x8_c;
        if (flags & HAS_SSSE3) eb_aom_smooth_h_predictor_4x8 = eb_aom_smooth_h_predictor_4x8_ssse3;
        eb_aom_smooth_h_predictor_64x16 = eb_aom_smooth_h_predictor_64x16_c;
        if (flags & HAS_SSSE3) eb_aom_smooth_h_predictor_64x16 = eb_aom_smooth_h_predictor_64x16_ssse3;
        eb_aom_smooth_h_predictor_64x32 = eb_aom_smooth_h_predictor_64x32_c;
        if (flags & HAS_SSSE3) eb_aom_smooth_h_predictor_64x32 = eb_aom_smooth_h_predictor_64x32_ssse3;
        eb_aom_smooth_h_predictor_8x16 = eb_aom_smooth_h_predictor_8x16_c;
        if (flags & HAS_SSSE3) eb_aom_smooth_h_predictor_8x16 = eb_aom_smooth_h_predictor_8x16_ssse3;
        eb_aom_smooth_h_predictor_8x32 = eb_aom_smooth_h_predictor_8x32_c;
        if (flags & HAS_SSSE3) eb_aom_smooth_h_predictor_8x32 = eb_aom_smooth_h_predictor_8x32_ssse3;
        eb_aom_smooth_h_predictor_8x4 = eb_aom_smooth_h_predictor_8x4_c;
        if (flags & HAS_SSSE3) eb_aom_smooth_h_predictor_8x4 = eb_aom_smooth_h_predictor_8x4_ssse3;
        eb_aom_smooth_h_predictor_64x64 = eb_aom_smooth_h_predictor_64x64_c;
        if (flags & HAS_SSSE3) eb_aom_smooth_h_predictor_64x64 = eb_aom_smooth_h_predictor_64x64_ssse3;
        eb_aom_smooth_h_predictor_32x32 = eb_aom_smooth_h_predictor_32x32_c;
        if (flags & HAS_SSSE3) eb_aom_smooth_h_predictor_32x32 = eb_aom_smooth_h_predictor_32x32_ssse3;
        eb_aom_smooth_h_predictor_16x16 = eb_aom_smooth_h_predictor_16x16_c;
        if (flags & HAS_SSSE3) eb_aom_smooth_h_predictor_16x16 = eb_aom_smooth_h_predictor_16x16_ssse3;
        eb_aom_smooth_h_predictor_8x8 = eb_aom_smooth_h_predictor_8x8_c;
        if (flags & HAS_SSSE3) eb_aom_smooth_h_predictor_8x8 = eb_aom_smooth_h_predictor_8x8_ssse3;
        eb_aom_smooth_h_predictor_4x4 = eb_aom_smooth_h_predictor_4x4_c;
        if (flags & HAS_SSSE3) eb_aom_smooth_h_predictor_4x4 = eb_aom_smooth_h_predictor_4x4_ssse3;

        if (flags & HAS_SSSE3) eb_aom_smooth_v_predictor_16x32 = eb_aom_smooth_v_predictor_16x32_ssse3;
        eb_aom_smooth_v_predictor_16x4 = eb_aom_smooth_v_predictor_16x4_c;
        if (flags & HAS_SSSE3) eb_aom_smooth_v_predictor_16x4 = eb_aom_smooth_v_predictor_16x4_ssse3;
        eb_aom_smooth_v_predictor_16x64 = eb_aom_smooth_v_predictor_16x64_c;
        if (flags & HAS_SSSE3) eb_aom_smooth_v_predictor_16x64 = eb_aom_smooth_v_predictor_16x64_ssse3;
        eb_aom_smooth_v_predictor_16x8 = eb_aom_smooth_v_predictor_16x8_c;
        if (flags & HAS_SSSE3) eb_aom_smooth_v_predictor_16x8 = eb_aom_smooth_v_predictor_16x8_ssse3;
        eb_aom_smooth_v_predictor_32x16 = eb_aom_smooth_v_predictor_32x16_c;
        if (flags & HAS_SSSE3) eb_aom_smooth_v_predictor_32x16 = eb_aom_smooth_v_predictor_32x16_ssse3;
        eb_aom_smooth_v_predictor_32x64 = eb_aom_smooth_v_predictor_32x64_c;
        if (flags & HAS_SSSE3) eb_aom_smooth_v_predictor_32x64 = eb_aom_smooth_v_predictor_32x64_ssse3;
        eb_aom_smooth_v_predictor_32x8 = eb_aom_smooth_v_predictor_32x8_c;
        if (flags & HAS_SSSE3) eb_aom_smooth_v_predictor_32x8 = eb_aom_smooth_v_predictor_32x8_ssse3;
        eb_aom_smooth_v_predictor_4x16 = eb_aom_smooth_v_predictor_4x16_c;
        if (flags & HAS_SSSE3) eb_aom_smooth_v_predictor_4x16 = eb_aom_smooth_v_predictor_4x16_ssse3;
        eb_aom_smooth_v_predictor_4x8 = eb_aom_smooth_v_predictor_4x8_c;
        if (flags & HAS_SSSE3) eb_aom_smooth_v_predictor_4x8 = eb_aom_smooth_v_predictor_4x8_ssse3;
        eb_aom_smooth_v_predictor_64x16 = eb_aom_smooth_v_predictor_64x16_c;
        if (flags & HAS_SSSE3) eb_aom_smooth_v_predictor_64x16 = eb_aom_smooth_v_predictor_64x16_ssse3;
        eb_aom_smooth_v_predictor_64x32 = eb_aom_smooth_v_predictor_64x32_c;
        if (flags & HAS_SSSE3) eb_aom_smooth_v_predictor_64x32 = eb_aom_smooth_v_predictor_64x32_ssse3;
        eb_aom_smooth_v_predictor_8x16 = eb_aom_smooth_v_predictor_8x16_c;
        if (flags & HAS_SSSE3) eb_aom_smooth_v_predictor_8x16 = eb_aom_smooth_v_predictor_8x16_ssse3;
        eb_aom_smooth_v_predictor_8x32 = eb_aom_smooth_v_predictor_8x32_c;
        if (flags & HAS_SSSE3) eb_aom_smooth_v_predictor_8x32 = eb_aom_smooth_v_predictor_8x32_ssse3;
        eb_aom_smooth_v_predictor_8x4 = eb_aom_smooth_v_predictor_8x4_c;
        if (flags & HAS_SSSE3) eb_aom_smooth_v_predictor_8x4 = eb_aom_smooth_v_predictor_8x4_ssse3;
        eb_aom_smooth_v_predictor_64x64 = eb_aom_smooth_v_predictor_64x64_c;
        if (flags & HAS_SSSE3) eb_aom_smooth_v_predictor_64x64 = eb_aom_smooth_v_predictor_64x64_ssse3;
        eb_aom_smooth_v_predictor_32x32 = eb_aom_smooth_v_predictor_32x32_c;
        if (flags & HAS_SSSE3) eb_aom_smooth_v_predictor_32x32 = eb_aom_smooth_v_predictor_32x32_ssse3;
        eb_aom_smooth_v_predictor_16x16 = eb_aom_smooth_v_predictor_16x16_c;
        if (flags & HAS_SSSE3) eb_aom_smooth_v_predictor_16x16 = eb_aom_smooth_v_predictor_16x16_ssse3;
        eb_aom_smooth_v_predictor_8x8 = eb_aom_smooth_v_predictor_8x8_c;
        if (flags & HAS_SSSE3) eb_aom_smooth_v_predictor_8x8 = eb_aom_smooth_v_predictor_8x8_ssse3;
        eb_aom_smooth_v_predictor_4x4 = eb_aom_smooth_v_predictor_4x4_c;
        if (flags & HAS_SSSE3) eb_aom_smooth_v_predictor_4x4 = eb_aom_smooth_v_predictor_4x4_ssse3;

        eb_aom_smooth_predictor_16x32 = eb_aom_smooth_predictor_16x32_c;
        if (flags & HAS_SSSE3) eb_aom_smooth_predictor_16x32 = eb_aom_smooth_predictor_16x32_ssse3;
        eb_aom_smooth_predictor_16x4 = eb_aom_smooth_predictor_16x4_c;
        if (flags & HAS_SSSE3) eb_aom_smooth_predictor_16x4 = eb_aom_smooth_predictor_16x4_ssse3;
        eb_aom_smooth_predictor_16x64 = eb_aom_smooth_predictor_16x64_c;
        if (flags & HAS_SSSE3) eb_aom_smooth_predictor_16x64 = eb_aom_smooth_predictor_16x64_ssse3;
        eb_aom_smooth_predictor_16x8 = eb_aom_smooth_predictor_16x8_c;
        if (flags & HAS_SSSE3) eb_aom_smooth_predictor_16x8 = eb_aom_smooth_predictor_16x8_ssse3;
        eb_aom_smooth_predictor_32x16 = eb_aom_smooth_predictor_32x16_c;
        if (flags & HAS_SSSE3) eb_aom_smooth_predictor_32x16 = eb_aom_smooth_predictor_32x16_ssse3;
        eb_aom_smooth_predictor_32x64 = eb_aom_smooth_predictor_32x64_c;
        if (flags & HAS_SSSE3) eb_aom_smooth_predictor_32x64 = eb_aom_smooth_predictor_32x64_ssse3;
        eb_aom_smooth_predictor_32x8 = eb_aom_smooth_predictor_32x8_c;
        if (flags & HAS_SSSE3) eb_aom_smooth_predictor_32x8 = eb_aom_smooth_predictor_32x8_ssse3;
        eb_aom_smooth_predictor_4x16 = eb_aom_smooth_predictor_4x16_c;
        if (flags & HAS_SSSE3) eb_aom_smooth_predictor_4x16 = eb_aom_smooth_predictor_4x16_ssse3;
        eb_aom_smooth_predictor_4x8 = eb_aom_smooth_predictor_4x8_c;
        if (flags & HAS_SSSE3) eb_aom_smooth_predictor_4x8 = eb_aom_smooth_predictor_4x8_ssse3;
        eb_aom_smooth_predictor_64x16 = eb_aom_smooth_predictor_64x16_c;
        if (flags & HAS_SSSE3) eb_aom_smooth_predictor_64x16 = eb_aom_smooth_predictor_64x16_ssse3;
        eb_aom_smooth_predictor_64x32 = eb_aom_smooth_predictor_64x32_c;
        if (flags & HAS_SSSE3) eb_aom_smooth_predictor_64x32 = eb_aom_smooth_predictor_64x32_ssse3;
        eb_aom_smooth_predictor_8x16 = eb_aom_smooth_predictor_8x16_c;
        if (flags & HAS_SSSE3) eb_aom_smooth_predictor_8x16 = eb_aom_smooth_predictor_8x16_ssse3;
        eb_aom_smooth_predictor_8x32 = eb_aom_smooth_predictor_8x32_c;
        if (flags & HAS_SSSE3) eb_aom_smooth_predictor_8x32 = eb_aom_smooth_predictor_8x32_ssse3;
        eb_aom_smooth_predictor_8x4 = eb_aom_smooth_predictor_8x4_c;
        if (flags & HAS_SSSE3) eb_aom_smooth_predictor_8x4 = eb_aom_smooth_predictor_8x4_ssse3;
        eb_aom_smooth_predictor_64x64 = eb_aom_smooth_predictor_64x64_c;
        if (flags & HAS_SSSE3) eb_aom_smooth_predictor_64x64 = eb_aom_smooth_predictor_64x64_ssse3;
        eb_aom_smooth_predictor_32x32 = eb_aom_smooth_predictor_32x32_c;
        if (flags & HAS_SSSE3) eb_aom_smooth_predictor_32x32 = eb_aom_smooth_predictor_32x32_ssse3;
        eb_aom_smooth_predictor_16x16 = eb_aom_smooth_predictor_16x16_c;
        if (flags & HAS_SSSE3) eb_aom_smooth_predictor_16x16 = eb_aom_smooth_predictor_16x16_ssse3;
        eb_aom_smooth_predictor_8x8 = eb_aom_smooth_predictor_8x8_c;
        if (flags & HAS_SSSE3) eb_aom_smooth_predictor_8x8 = eb_aom_smooth_predictor_8x8_ssse3;
        eb_aom_smooth_predictor_4x4 = eb_aom_smooth_predictor_4x4_c;
        if (flags & HAS_SSSE3) eb_aom_smooth_predictor_4x4 = eb_aom_smooth_predictor_4x4_ssse3;

        eb_aom_v_predictor_4x4 = eb_aom_v_predictor_4x4_c;
        if (flags & HAS_SSE2) eb_aom_v_predictor_4x4 = eb_aom_v_predictor_4x4_sse2;
        eb_aom_v_predictor_8x8 = eb_aom_v_predictor_8x8_c;
        if (flags & HAS_SSE2) eb_aom_v_predictor_8x8 = eb_aom_v_predictor_8x8_sse2;
        eb_aom_v_predictor_16x16 = eb_aom_v_predictor_16x16_c;
        if (flags & HAS_SSE2) eb_aom_v_predictor_16x16 = eb_aom_v_predictor_16x16_sse2;
        eb_aom_v_predictor_32x32 = eb_aom_v_predictor_32x32_c;
        if (flags & HAS_AVX2) eb_aom_v_predictor_32x32 = eb_aom_v_predictor_32x32_avx2;
        eb_aom_v_predictor_64x64 = eb_aom_v_predictor_64x64_c;
        if (flags & HAS_AVX2) eb_aom_v_predictor_64x64 = eb_aom_v_predictor_64x64_avx2;
        eb_aom_v_predictor_16x32 = eb_aom_v_predictor_16x32_c;
        if (flags & HAS_SSE2) eb_aom_v_predictor_16x32 = eb_aom_v_predictor_16x32_sse2;
        eb_aom_v_predictor_16x4 = eb_aom_v_predictor_16x4_c;
        if (flags & HAS_SSE2) eb_aom_v_predictor_16x4 = eb_aom_v_predictor_16x4_sse2;
        eb_aom_v_predictor_16x64 = eb_aom_v_predictor_16x64_c;
        if (flags & HAS_SSE2) eb_aom_v_predictor_16x64 = eb_aom_v_predictor_16x64_sse2;
        eb_aom_v_predictor_16x8 = eb_aom_v_predictor_16x8_c;
        if (flags & HAS_SSE2) eb_aom_v_predictor_16x8 = eb_aom_v_predictor_16x8_sse2;
        eb_aom_v_predictor_32x16 = eb_aom_v_predictor_32x16_c;
        if (flags & HAS_AVX2) eb_aom_v_predictor_32x16 = eb_aom_v_predictor_32x16_avx2;
        eb_aom_v_predictor_32x64 = eb_aom_v_predictor_32x64_c;
        if (flags & HAS_AVX2) eb_aom_v_predictor_32x64 = eb_aom_v_predictor_32x64_avx2;
        eb_aom_v_predictor_32x8 = eb_aom_v_predictor_32x8_c;
        if (flags & HAS_SSE2) eb_aom_v_predictor_32x8 = eb_aom_v_predictor_32x8_sse2;
        eb_aom_v_predictor_4x16 = eb_aom_v_predictor_4x16_c;
        if (flags & HAS_SSE2) eb_aom_v_predictor_4x16 = eb_aom_v_predictor_4x16_sse2;
        eb_aom_v_predictor_4x8 = eb_aom_v_predictor_4x8_c;
        if (flags & HAS_SSE2) eb_aom_v_predictor_4x8 = eb_aom_v_predictor_4x8_sse2;
        eb_aom_v_predictor_64x16 = eb_aom_v_predictor_64x16_c;
        if (flags & HAS_AVX2) eb_aom_v_predictor_64x16 = eb_aom_v_predictor_64x16_avx2;
        eb_aom_v_predictor_64x32 = eb_aom_v_predictor_64x32_c;
        if (flags & HAS_AVX2) eb_aom_v_predictor_64x32 = eb_aom_v_predictor_64x32_avx2;
        eb_aom_v_predictor_8x16 = eb_aom_v_predictor_8x16_c;
        if (flags & HAS_SSE2) eb_aom_v_predictor_8x16 = eb_aom_v_predictor_8x16_sse2;
        eb_aom_v_predictor_8x32 = eb_aom_v_predictor_8x32_c;
        if (flags & HAS_SSE2) eb_aom_v_predictor_8x32 = eb_aom_v_predictor_8x32_sse2;
        eb_aom_v_predictor_8x4 = eb_aom_v_predictor_8x4_c;
        if (flags & HAS_SSE2) eb_aom_v_predictor_8x4 = eb_aom_v_predictor_8x4_sse2;

        eb_aom_h_predictor_4x4 = eb_aom_h_predictor_4x4_c;
        if (flags & HAS_SSE2) eb_aom_h_predictor_4x4 = eb_aom_h_predictor_4x4_sse2;
        eb_aom_h_predictor_8x8 = eb_aom_h_predictor_8x8_c;
        if (flags & HAS_SSE2) eb_aom_h_predictor_8x8 = eb_aom_h_predictor_8x8_sse2;
        eb_aom_h_predictor_16x16 = eb_aom_h_predictor_16x16_c;
        if (flags & HAS_SSE2) eb_aom_h_predictor_16x16 = eb_aom_h_predictor_16x16_sse2;
        eb_aom_h_predictor_32x32 = eb_aom_h_predictor_32x32_c;
        if (flags & HAS_AVX2) eb_aom_h_predictor_32x32 = eb_aom_h_predictor_32x32_avx2;
        eb_aom_h_predictor_64x64 = eb_aom_h_predictor_64x64_c;
        if (flags & HAS_SSE2) eb_aom_h_predictor_64x64 = eb_aom_h_predictor_64x64_sse2;
        eb_aom_h_predictor_16x32 = eb_aom_h_predictor_16x32_c;
        if (flags & HAS_SSE2) eb_aom_h_predictor_16x32 = eb_aom_h_predictor_16x32_sse2;
        eb_aom_h_predictor_16x4 = eb_aom_h_predictor_16x4_c;
        if (flags & HAS_SSE2) eb_aom_h_predictor_16x4 = eb_aom_h_predictor_16x4_sse2;
        eb_aom_h_predictor_16x64 = eb_aom_h_predictor_16x64_c;
        if (flags & HAS_SSE2) eb_aom_h_predictor_16x64 = eb_aom_h_predictor_16x64_sse2;
        eb_aom_h_predictor_16x8 = eb_aom_h_predictor_16x8_c;
        if (flags & HAS_SSE2) eb_aom_h_predictor_16x8 = eb_aom_h_predictor_16x8_sse2;
        eb_aom_h_predictor_32x16 = eb_aom_h_predictor_32x16_c;
        if (flags & HAS_SSE2) eb_aom_h_predictor_32x16 = eb_aom_h_predictor_32x16_sse2;
        eb_aom_h_predictor_32x64 = eb_aom_h_predictor_32x64_c;
        if (flags & HAS_SSE2) eb_aom_h_predictor_32x64 = eb_aom_h_predictor_32x64_sse2;
        eb_aom_h_predictor_32x8 = eb_aom_h_predictor_32x8_c;
        if (flags & HAS_SSE2) eb_aom_h_predictor_32x8 = eb_aom_h_predictor_32x8_sse2;
        eb_aom_h_predictor_4x16 = eb_aom_h_predictor_4x16_c;
        if (flags & HAS_SSE2) eb_aom_h_predictor_4x16 = eb_aom_h_predictor_4x16_sse2;
        eb_aom_h_predictor_4x8 = eb_aom_h_predictor_4x8_c;
        if (flags & HAS_SSE2) eb_aom_h_predictor_4x8 = eb_aom_h_predictor_4x8_sse2;
        eb_aom_h_predictor_64x16 = eb_aom_h_predictor_64x16_c;
        if (flags & HAS_SSE2) eb_aom_h_predictor_64x16 = eb_aom_h_predictor_64x16_sse2;
        eb_aom_h_predictor_64x32 = eb_aom_h_predictor_64x32_c;
        if (flags & HAS_SSE2) eb_aom_h_predictor_64x32 = eb_aom_h_predictor_64x32_sse2;
        eb_aom_h_predictor_8x16 = eb_aom_h_predictor_8x16_c;
        if (flags & HAS_SSE2) eb_aom_h_predictor_8x16 = eb_aom_h_predictor_8x16_sse2;
        eb_aom_h_predictor_8x32 = eb_aom_h_predictor_8x32_c;
        if (flags & HAS_SSE2) eb_aom_h_predictor_8x32 = eb_aom_h_predictor_8x32_sse2;
        eb_aom_h_predictor_8x4 = eb_aom_h_predictor_8x4_c;
        if (flags & HAS_SSE2) eb_aom_h_predictor_8x4 = eb_aom_h_predictor_8x4_sse2;

        //SAD
        eb_aom_sad4x4 = eb_aom_sad4x4_c;
        if (flags & HAS_AVX2) eb_aom_sad4x4 = eb_aom_sad4x4_avx2;
        eb_aom_sad4x4x4d = eb_aom_sad4x4x4d_c;
        if (flags & HAS_AVX2) eb_aom_sad4x4x4d = eb_aom_sad4x4x4d_avx2;
        eb_aom_sad4x16 = eb_aom_sad4x16_c;
        if (flags & HAS_AVX2) eb_aom_sad4x16 = eb_aom_sad4x16_avx2;
        eb_aom_sad4x16x4d = eb_aom_sad4x16x4d_c;
        if (flags & HAS_AVX2) eb_aom_sad4x16x4d = eb_aom_sad4x16x4d_avx2;
        eb_aom_sad4x8 = eb_aom_sad4x8_c;
        if (flags & HAS_AVX2) eb_aom_sad4x8 = eb_aom_sad4x8_avx2;
        eb_aom_sad4x8x4d = eb_aom_sad4x8x4d_c;
        if (flags & HAS_AVX2) eb_aom_sad4x8x4d = eb_aom_sad4x8x4d_avx2;
        eb_aom_sad64x128 = eb_aom_sad64x128_c;
        if (flags & HAS_AVX2) eb_aom_sad64x128 = eb_aom_sad64x128_avx2;
        eb_aom_sad64x128x4d = eb_aom_sad64x128x4d_c;
        if (flags & HAS_AVX2) eb_aom_sad64x128x4d = eb_aom_sad64x128x4d_avx2;
        eb_aom_sad64x16 = eb_aom_sad64x16_c;
        if (flags & HAS_AVX2) eb_aom_sad64x16 = eb_aom_sad64x16_avx2;
        eb_aom_sad64x16x4d = eb_aom_sad64x16x4d_c;
        if (flags & HAS_AVX2) eb_aom_sad64x16x4d = eb_aom_sad64x16x4d_avx2;
        eb_aom_sad64x32 = eb_aom_sad64x32_c;
        if (flags & HAS_AVX2) eb_aom_sad64x32 = eb_aom_sad64x32_avx2;
        eb_aom_sad64x32x4d = eb_aom_sad64x32x4d_c;
        if (flags & HAS_AVX2) eb_aom_sad64x32x4d = eb_aom_sad64x32x4d_avx2;
        eb_aom_sad64x64 = eb_aom_sad64x64_c;
        if (flags & HAS_AVX2) eb_aom_sad64x64 = eb_aom_sad64x64_avx2;
        eb_aom_sad64x64x4d = eb_aom_sad64x64x4d_c;
        if (flags & HAS_AVX2) eb_aom_sad64x64x4d = eb_aom_sad64x64x4d_avx2;
        eb_aom_sad8x16 = eb_aom_sad8x16_c;
        if (flags & HAS_AVX2) eb_aom_sad8x16 = eb_aom_sad8x16_avx2;
        eb_aom_sad8x16x4d = eb_aom_sad8x16x4d_c;
        if (flags & HAS_AVX2) eb_aom_sad8x16x4d = eb_aom_sad8x16x4d_avx2;
        eb_aom_sad8x32 = eb_aom_sad8x32_c;
        if (flags & HAS_AVX2) eb_aom_sad8x32 = eb_aom_sad8x32_avx2;
        eb_aom_sad8x32x4d = eb_aom_sad8x32x4d_c;
        if (flags & HAS_AVX2) eb_aom_sad8x32x4d = eb_aom_sad8x32x4d_avx2;
        eb_aom_sad8x8 = eb_aom_sad8x8_c;
        if (flags & HAS_AVX2) eb_aom_sad8x8 = eb_aom_sad8x8_avx2;
        eb_aom_sad8x8x4d = eb_aom_sad8x8x4d_c;
        if (flags & HAS_AVX2) eb_aom_sad8x8x4d = eb_aom_sad8x8x4d_avx2;
        eb_aom_sad16x4 = eb_aom_sad16x4_c;
        if (flags & HAS_AVX2) eb_aom_sad16x4 = eb_aom_sad16x4_avx2;
        eb_aom_sad16x4x4d = eb_aom_sad16x4x4d_c;
        if (flags & HAS_AVX2) eb_aom_sad16x4x4d = eb_aom_sad16x4x4d_avx2;
        eb_aom_sad32x8 = eb_aom_sad32x8_c;
        if (flags & HAS_AVX2) eb_aom_sad32x8 = eb_aom_sad32x8_avx2;
        eb_aom_sad32x8x4d = eb_aom_sad32x8x4d_c;
        if (flags & HAS_AVX2) eb_aom_sad32x8x4d = eb_aom_sad32x8x4d_avx2;
        eb_aom_sad16x64 = eb_aom_sad16x64_c;
        if (flags & HAS_AVX2) eb_aom_sad16x64 = eb_aom_sad16x64_avx2;
        eb_aom_sad16x64x4d = eb_aom_sad16x64x4d_c;
        if (flags & HAS_AVX2) eb_aom_sad16x64x4d = eb_aom_sad16x64x4d_avx2;
        eb_aom_sad128x128 = eb_aom_sad128x128_c;
        if (flags & HAS_AVX2) eb_aom_sad128x128 = eb_aom_sad128x128_avx2;
        eb_aom_sad128x128x4d = eb_aom_sad128x128x4d_c;
        if (flags & HAS_AVX2) eb_aom_sad128x128x4d = eb_aom_sad128x128x4d_avx2;
        eb_aom_sad128x64 = eb_aom_sad128x64_c;
        if (flags & HAS_AVX2) eb_aom_sad128x64 = eb_aom_sad128x64_avx2;
        eb_aom_sad128x64x4d = eb_aom_sad128x64x4d_c;
        if (flags & HAS_AVX2) eb_aom_sad128x64x4d = eb_aom_sad128x64x4d_avx2;
        eb_aom_sad32x16 = eb_aom_sad32x16_c;
        if (flags & HAS_AVX2) eb_aom_sad32x16 = eb_aom_sad32x16_avx2;
        eb_aom_sad32x16x4d = eb_aom_sad32x16x4d_c;
        if (flags & HAS_AVX2) eb_aom_sad32x16x4d = eb_aom_sad32x16x4d_avx2;
        eb_aom_sad16x32 = eb_aom_sad16x32_c;
        if (flags & HAS_AVX2) eb_aom_sad16x32 = eb_aom_sad16x32_avx2;
        eb_aom_sad16x32x4d = eb_aom_sad16x32x4d_c;
        if (flags & HAS_AVX2) eb_aom_sad16x32x4d = eb_aom_sad16x32x4d_avx2;
        eb_aom_sad32x64 = eb_aom_sad32x64_c;
        if (flags & HAS_AVX2) eb_aom_sad32x64 = eb_aom_sad32x64_avx2;
        eb_aom_sad32x64x4d = eb_aom_sad32x64x4d_c;
        if (flags & HAS_AVX2) eb_aom_sad32x64x4d = eb_aom_sad32x64x4d_avx2;
        eb_aom_sad32x32 = eb_aom_sad32x32_c;
        if (flags & HAS_AVX2) eb_aom_sad32x32 = eb_aom_sad32x32_avx2;
        eb_aom_sad32x32x4d = eb_aom_sad32x32x4d_c;
        if (flags & HAS_AVX2) eb_aom_sad32x32x4d = eb_aom_sad32x32x4d_avx2;
        eb_aom_sad16x16 = eb_aom_sad16x16_c;
        if (flags & HAS_AVX2) eb_aom_sad16x16 = eb_aom_sad16x16_avx2;
        eb_aom_sad16x16x4d = eb_aom_sad16x16x4d_c;
        if (flags & HAS_AVX2) eb_aom_sad16x16x4d = eb_aom_sad16x16x4d_avx2;
        eb_aom_sad16x8 = eb_aom_sad16x8_c;
        if (flags & HAS_AVX2) eb_aom_sad16x8 = eb_aom_sad16x8_avx2;
        eb_aom_sad16x8x4d = eb_aom_sad16x8x4d_c;
        if (flags & HAS_AVX2) eb_aom_sad16x8x4d = eb_aom_sad16x8x4d_avx2;
        eb_aom_sad8x4 = eb_aom_sad8x4_c;
        if (flags & HAS_AVX2) eb_aom_sad8x4 = eb_aom_sad8x4_avx2;
        eb_aom_sad8x4x4d = eb_aom_sad8x4x4d_c;
        if (flags & HAS_AVX2) eb_aom_sad8x4x4d = eb_aom_sad8x4x4d_avx2;

//VARIANCE
        eb_aom_variance4x4 = eb_aom_variance4x4_c;
        if (flags & HAS_AVX2) eb_aom_variance4x4 = eb_aom_variance4x4_sse2;
        eb_aom_variance4x8 = eb_aom_variance4x8_c;
        if (flags & HAS_AVX2) eb_aom_variance4x8 = eb_aom_variance4x8_sse2;
        eb_aom_variance4x16 = eb_aom_variance4x16_c;
        if (flags & HAS_AVX2) eb_aom_variance4x16 = eb_aom_variance4x16_sse2;
        eb_aom_variance8x4 = eb_aom_variance8x4_c;
        if (flags & HAS_AVX2) eb_aom_variance8x4 = eb_aom_variance8x4_sse2;
        eb_aom_variance8x8 = eb_aom_variance8x8_c;
        if (flags & HAS_AVX2) eb_aom_variance8x8 = eb_aom_variance8x8_sse2;
        eb_aom_variance8x16 = eb_aom_variance8x16_c;
        if (flags & HAS_AVX2) eb_aom_variance8x16 = eb_aom_variance8x16_sse2;
        eb_aom_variance8x32 = eb_aom_variance8x32_c;
        if (flags & HAS_AVX2) eb_aom_variance8x32 = eb_aom_variance8x32_sse2;
        eb_aom_variance16x4 = eb_aom_variance16x4_c;
        if (flags & HAS_AVX2) eb_aom_variance16x4 = eb_aom_variance16x4_avx2;
        eb_aom_variance16x8 = eb_aom_variance16x8_c;
        if (flags & HAS_AVX2) eb_aom_variance16x8 = eb_aom_variance16x8_avx2;
        eb_aom_variance16x16 = eb_aom_variance16x16_c;
        if (flags & HAS_AVX2) eb_aom_variance16x16 = eb_aom_variance16x16_avx2;
        eb_aom_variance16x32 = eb_aom_variance16x32_c;
        if (flags & HAS_AVX2) eb_aom_variance16x32 = eb_aom_variance16x32_avx2;
        eb_aom_variance16x64 = eb_aom_variance16x64_c;
        if (flags & HAS_AVX2) eb_aom_variance16x64 = eb_aom_variance16x64_avx2;
        eb_aom_variance32x8 = eb_aom_variance32x8_c;
        if (flags & HAS_AVX2) eb_aom_variance32x8 = eb_aom_variance32x8_avx2;
        eb_aom_variance32x16 = eb_aom_variance32x16_c;
        if (flags & HAS_AVX2) eb_aom_variance32x16 = eb_aom_variance32x16_avx2;
        eb_aom_variance32x32 = eb_aom_variance32x32_c;
        if (flags & HAS_AVX2) eb_aom_variance32x32 = eb_aom_variance32x32_avx2;
        eb_aom_variance32x64 = eb_aom_variance32x64_c;
        if (flags & HAS_AVX2) eb_aom_variance32x64 = eb_aom_variance32x64_avx2;
        eb_aom_variance64x16 = eb_aom_variance64x16_c;
        if (flags & HAS_AVX2) eb_aom_variance64x16 = eb_aom_variance64x16_avx2;
        eb_aom_variance64x32 = eb_aom_variance64x32_c;
        if (flags & HAS_AVX2) eb_aom_variance64x32 = eb_aom_variance64x32_avx2;
        eb_aom_variance64x64 = eb_aom_variance64x64_c;
        if (flags & HAS_AVX2) eb_aom_variance64x64 = eb_aom_variance64x64_avx2;
        eb_aom_variance64x128 = eb_aom_variance64x128_c;
        if (flags & HAS_AVX2) eb_aom_variance64x128 = eb_aom_variance64x128_avx2;
        eb_aom_variance128x64 = eb_aom_variance128x64_c;
        if (flags & HAS_AVX2) eb_aom_variance128x64 = eb_aom_variance128x64_avx2;
        eb_aom_variance128x128 = eb_aom_variance128x128_c;
        if (flags & HAS_AVX2) eb_aom_variance128x128 = eb_aom_variance128x128_avx2;

        //QIQ
        eb_aom_quantize_b_64x64 = eb_aom_quantize_b_64x64_c_II;
        if (flags & HAS_AVX2) eb_aom_quantize_b_64x64 = eb_aom_highbd_quantize_b_64x64_avx2;

        eb_aom_highbd_quantize_b_64x64 = eb_aom_highbd_quantize_b_64x64_c;
        if (flags & HAS_AVX2) eb_aom_highbd_quantize_b_64x64 = eb_aom_highbd_quantize_b_64x64_avx2;
        // transform
        eb_av1_fwd_txfm2d_16x8 = eb_av1_fwd_txfm2d_16x8_c;
        if (flags & HAS_AVX2) eb_av1_fwd_txfm2d_16x8 = eb_av1_fwd_txfm2d_16x8_avx2;
        eb_av1_fwd_txfm2d_8x16 = eb_av1_fwd_txfm2d_8x16_c;
        if (flags & HAS_AVX2) eb_av1_fwd_txfm2d_8x16 = eb_av1_fwd_txfm2d_8x16_avx2;

        eb_av1_fwd_txfm2d_16x4 = eb_av1_fwd_txfm2d_16x4_c;
        if (flags & HAS_AVX2) eb_av1_fwd_txfm2d_16x4 = eb_av1_fwd_txfm2d_16x4_avx2;
        eb_av1_fwd_txfm2d_4x16 = eb_av1_fwd_txfm2d_4x16_c;
        if (flags & HAS_AVX2) eb_av1_fwd_txfm2d_4x16 = eb_av1_fwd_txfm2d_4x16_avx2;

        eb_av1_fwd_txfm2d_8x4 = eb_av1_fwd_txfm2d_8x4_c;
        if (flags & HAS_AVX2) eb_av1_fwd_txfm2d_8x4 = eb_av1_fwd_txfm2d_8x4_avx2;
        eb_av1_fwd_txfm2d_4x8 = eb_av1_fwd_txfm2d_4x8_c;
        if (flags & HAS_AVX2) eb_av1_fwd_txfm2d_4x8 = eb_av1_fwd_txfm2d_4x8_avx2;

        eb_av1_fwd_txfm2d_32x16 = eb_av1_fwd_txfm2d_32x16_c;
        if (flags & HAS_AVX2) eb_av1_fwd_txfm2d_32x16 = eb_av1_fwd_txfm2d_32x16_avx2;
        eb_av1_fwd_txfm2d_32x8 = eb_av1_fwd_txfm2d_32x8_c;
        if (flags & HAS_AVX2) eb_av1_fwd_txfm2d_32x8 = eb_av1_fwd_txfm2d_32x8_avx2;
        eb_av1_fwd_txfm2d_8x32 = eb_av1_fwd_txfm2d_8x32_c;
        if (flags & HAS_AVX2) eb_av1_fwd_txfm2d_8x32 = eb_av1_fwd_txfm2d_8x32_avx2;
        eb_av1_fwd_txfm2d_16x32 = eb_av1_fwd_txfm2d_16x32_c;
        if (flags & HAS_AVX2) eb_av1_fwd_txfm2d_16x32 = eb_av1_fwd_txfm2d_16x32_avx2;
        eb_av1_fwd_txfm2d_32x64 = eb_av1_fwd_txfm2d_32x64_c;
        if (flags & HAS_AVX2) eb_av1_fwd_txfm2d_32x64 = eb_av1_fwd_txfm2d_32x64_avx2;
        eb_av1_fwd_txfm2d_64x32 = eb_av1_fwd_txfm2d_64x32_c;
        if (flags & HAS_AVX2) eb_av1_fwd_txfm2d_64x32 = eb_av1_fwd_txfm2d_64x32_avx2;
        eb_av1_fwd_txfm2d_16x64 = eb_av1_fwd_txfm2d_16x64_c;
        if (flags & HAS_AVX2) eb_av1_fwd_txfm2d_16x64 = eb_av1_fwd_txfm2d_16x64_avx2;
        eb_av1_fwd_txfm2d_64x16 = eb_av1_fwd_txfm2d_64x16_c;
        if (flags & HAS_AVX2) eb_av1_fwd_txfm2d_64x16 = eb_av1_fwd_txfm2d_64x16_avx2;
        eb_av1_fwd_txfm2d_64x64 = Av1TransformTwoD_64x64_c;
        if (flags & HAS_AVX2) eb_av1_fwd_txfm2d_64x64 = eb_av1_fwd_txfm2d_64x64_avx2;
        eb_av1_fwd_txfm2d_32x32 = Av1TransformTwoD_32x32_c;
        if (flags & HAS_AVX2) eb_av1_fwd_txfm2d_32x32 = eb_av1_fwd_txfm2d_32x32_avx2;
        eb_av1_fwd_txfm2d_16x16 = Av1TransformTwoD_16x16_c;
        if (flags & HAS_AVX2) eb_av1_fwd_txfm2d_16x16 = eb_av1_fwd_txfm2d_16x16_avx2;
        eb_av1_fwd_txfm2d_8x8 = Av1TransformTwoD_8x8_c;
        if (flags & HAS_AVX2) eb_av1_fwd_txfm2d_8x8 = eb_av1_fwd_txfm2d_8x8_avx2;
        eb_av1_fwd_txfm2d_4x4 = Av1TransformTwoD_4x4_c;
        if (flags & HAS_SSE4_1) eb_av1_fwd_txfm2d_4x4 = eb_av1_fwd_txfm2d_4x4_sse4_1;

        HandleTransform16x64 = HandleTransform16x64_c;
        if (flags & HAS_AVX2) HandleTransform16x64 = HandleTransform16x64_avx2;
        HandleTransform32x64 = HandleTransform32x64_c;
        if (flags & HAS_AVX2) HandleTransform32x64 = HandleTransform32x64_avx2;
        HandleTransform64x16 = HandleTransform64x16_c;
        if (flags & HAS_AVX2) HandleTransform64x16 = HandleTransform64x16_avx2;
        HandleTransform64x32 = HandleTransform64x32_c;
        if (flags & HAS_AVX2) HandleTransform64x32 = HandleTransform64x32_avx2;
        HandleTransform64x64 = HandleTransform64x64_c;
        if (flags & HAS_AVX2) HandleTransform64x64 = HandleTransform64x64_avx2;

        // eb_aom_highbd_v_predictor
        eb_aom_highbd_v_predictor_16x16 = eb_aom_highbd_v_predictor_16x16_c;
        if (flags & HAS_AVX2) eb_aom_highbd_v_predictor_16x16 = eb_aom_highbd_v_predictor_16x16_avx2;
        eb_aom_highbd_v_predictor_16x32 = eb_aom_highbd_v_predictor_16x32_c;
        if (flags & HAS_AVX2) eb_aom_highbd_v_predictor_16x32 = eb_aom_highbd_v_predictor_16x32_avx2;
        eb_aom_highbd_v_predictor_16x4 = eb_aom_highbd_v_predictor_16x4_c;
        if (flags & HAS_AVX2) eb_aom_highbd_v_predictor_16x4 = eb_aom_highbd_v_predictor_16x4_avx2;
        eb_aom_highbd_v_predictor_16x64 = eb_aom_highbd_v_predictor_16x64_c;
        if (flags & HAS_AVX2) eb_aom_highbd_v_predictor_16x64 = eb_aom_highbd_v_predictor_16x64_avx2;
        eb_aom_highbd_v_predictor_16x8 = eb_aom_highbd_v_predictor_16x8_c;
        if (flags & HAS_AVX2) eb_aom_highbd_v_predictor_16x8 = eb_aom_highbd_v_predictor_16x8_avx2;
        eb_aom_highbd_v_predictor_2x2 = eb_aom_highbd_v_predictor_2x2_c;
        eb_aom_highbd_v_predictor_32x16 = eb_aom_highbd_v_predictor_32x16_c;
        if (flags & HAS_AVX2) eb_aom_highbd_v_predictor_32x16 = eb_aom_highbd_v_predictor_32x16_avx2;
        eb_aom_highbd_v_predictor_32x32 = eb_aom_highbd_v_predictor_32x32_c;
        if (flags & HAS_AVX2) eb_aom_highbd_v_predictor_32x32 = eb_aom_highbd_v_predictor_32x32_avx2;
        eb_aom_highbd_v_predictor_32x64 = eb_aom_highbd_v_predictor_32x64_c;
        if (flags & HAS_AVX2) eb_aom_highbd_v_predictor_32x64 = eb_aom_highbd_v_predictor_32x64_avx2;
        eb_aom_highbd_v_predictor_32x8 = eb_aom_highbd_v_predictor_32x8_c;
        if (flags & HAS_AVX2) eb_aom_highbd_v_predictor_32x8 = eb_aom_highbd_v_predictor_32x8_avx2;
        eb_aom_highbd_v_predictor_4x16 = eb_aom_highbd_v_predictor_4x16_c;
        if (flags & HAS_SSE2) eb_aom_highbd_v_predictor_4x16 = eb_aom_highbd_v_predictor_4x16_sse2;
        eb_aom_highbd_v_predictor_4x4 = eb_aom_highbd_v_predictor_4x4_c;
        if (flags & HAS_SSE2) eb_aom_highbd_v_predictor_4x4 = eb_aom_highbd_v_predictor_4x4_sse2;
        eb_aom_highbd_v_predictor_4x8 = eb_aom_highbd_v_predictor_4x8_c;
        if (flags & HAS_SSE2) eb_aom_highbd_v_predictor_4x8 = eb_aom_highbd_v_predictor_4x8_sse2;
        eb_aom_highbd_v_predictor_64x16 = eb_aom_highbd_v_predictor_64x16_c;
        if (flags & HAS_AVX2) eb_aom_highbd_v_predictor_64x16 = eb_aom_highbd_v_predictor_64x16_avx2;
        eb_aom_highbd_v_predictor_64x32 = eb_aom_highbd_v_predictor_64x32_c;
        if (flags & HAS_AVX2) eb_aom_highbd_v_predictor_64x32 = eb_aom_highbd_v_predictor_64x32_avx2;
        eb_aom_highbd_v_predictor_8x32 = eb_aom_highbd_v_predictor_8x32_c;
        if (flags & HAS_SSE2) eb_aom_highbd_v_predictor_8x32 = eb_aom_highbd_v_predictor_8x32_sse2;
        eb_aom_highbd_v_predictor_64x64 = eb_aom_highbd_v_predictor_64x64_c;
        if (flags & HAS_AVX2) eb_aom_highbd_v_predictor_64x64 = eb_aom_highbd_v_predictor_64x64_avx2;
        eb_aom_highbd_v_predictor_8x16 = eb_aom_highbd_v_predictor_8x16_c;
        if (flags & HAS_SSE2) eb_aom_highbd_v_predictor_8x16 = eb_aom_highbd_v_predictor_8x16_sse2;
        eb_aom_highbd_v_predictor_8x4 = eb_aom_highbd_v_predictor_8x4_c;
        if (flags & HAS_SSE2) eb_aom_highbd_v_predictor_8x4 = eb_aom_highbd_v_predictor_8x4_sse2;
        eb_aom_highbd_v_predictor_8x8 = eb_aom_highbd_v_predictor_8x8_c;
        if (flags & HAS_SSE2) eb_aom_highbd_v_predictor_8x8 = eb_aom_highbd_v_predictor_8x8_sse2;

        //aom_highbd_smooth_predictor
        eb_aom_highbd_smooth_predictor_16x16 = eb_aom_highbd_smooth_predictor_16x16_c;
        if (flags & HAS_AVX2) eb_aom_highbd_smooth_predictor_16x16 = eb_aom_highbd_smooth_predictor_16x16_avx2;
        eb_aom_highbd_smooth_predictor_16x32 = eb_aom_highbd_smooth_predictor_16x32_c;
        if (flags & HAS_AVX2) eb_aom_highbd_smooth_predictor_16x32 = eb_aom_highbd_smooth_predictor_16x32_avx2;
        eb_aom_highbd_smooth_predictor_16x4 = eb_aom_highbd_smooth_predictor_16x4_c;
        if (flags & HAS_AVX2) eb_aom_highbd_smooth_predictor_16x4 = eb_aom_highbd_smooth_predictor_16x4_avx2;
        eb_aom_highbd_smooth_predictor_16x64 = eb_aom_highbd_smooth_predictor_16x64_c;
        if (flags & HAS_AVX2) eb_aom_highbd_smooth_predictor_16x64 = eb_aom_highbd_smooth_predictor_16x64_avx2;
        eb_aom_highbd_smooth_predictor_16x8 = eb_aom_highbd_smooth_predictor_16x8_c;
        if (flags & HAS_AVX2) eb_aom_highbd_smooth_predictor_16x8 = eb_aom_highbd_smooth_predictor_16x8_avx2;
        eb_aom_highbd_smooth_predictor_2x2 = eb_aom_highbd_smooth_predictor_2x2_c;
        eb_aom_highbd_smooth_predictor_32x16 = eb_aom_highbd_smooth_predictor_32x16_c;
        if (flags & HAS_AVX2) eb_aom_highbd_smooth_predictor_32x16 = eb_aom_highbd_smooth_predictor_32x16_avx2;
        eb_aom_highbd_smooth_predictor_32x32 = eb_aom_highbd_smooth_predictor_32x32_c;
        if (flags & HAS_AVX2) eb_aom_highbd_smooth_predictor_32x32 = eb_aom_highbd_smooth_predictor_32x32_avx2;
        eb_aom_highbd_smooth_predictor_32x64 = eb_aom_highbd_smooth_predictor_32x64_c;
        if (flags & HAS_AVX2) eb_aom_highbd_smooth_predictor_32x64 = eb_aom_highbd_smooth_predictor_32x64_avx2;
        eb_aom_highbd_smooth_predictor_32x8 = eb_aom_highbd_smooth_predictor_32x8_c;
        if (flags & HAS_AVX2) eb_aom_highbd_smooth_predictor_32x8 = eb_aom_highbd_smooth_predictor_32x8_avx2;
        eb_aom_highbd_smooth_predictor_4x16 = eb_aom_highbd_smooth_predictor_4x16_c;
        if (flags & HAS_SSSE3) eb_aom_highbd_smooth_predictor_4x16 = eb_aom_highbd_smooth_predictor_4x16_ssse3;
        eb_aom_highbd_smooth_predictor_4x4 = eb_aom_highbd_smooth_predictor_4x4_c;
        if (flags & HAS_SSSE3) eb_aom_highbd_smooth_predictor_4x4 = eb_aom_highbd_smooth_predictor_4x4_ssse3;
        eb_aom_highbd_smooth_predictor_4x8 = eb_aom_highbd_smooth_predictor_4x8_c;
        if (flags & HAS_SSSE3) eb_aom_highbd_smooth_predictor_4x8 = eb_aom_highbd_smooth_predictor_4x8_ssse3;
        eb_aom_highbd_smooth_predictor_64x16 = eb_aom_highbd_smooth_predictor_64x16_c;
        if (flags & HAS_AVX2) eb_aom_highbd_smooth_predictor_64x16 = eb_aom_highbd_smooth_predictor_64x16_avx2;
        eb_aom_highbd_smooth_predictor_64x32 = eb_aom_highbd_smooth_predictor_64x32_c;
        if (flags & HAS_AVX2) eb_aom_highbd_smooth_predictor_64x32 = eb_aom_highbd_smooth_predictor_64x32_avx2;
        eb_aom_highbd_smooth_predictor_64x64 = eb_aom_highbd_smooth_predictor_64x64_c;
        if (flags & HAS_AVX2) eb_aom_highbd_smooth_predictor_64x64 = eb_aom_highbd_smooth_predictor_64x64_avx2;
        eb_aom_highbd_smooth_predictor_8x16 = eb_aom_highbd_smooth_predictor_8x16_c;
        if (flags & HAS_AVX2) eb_aom_highbd_smooth_predictor_8x16 = eb_aom_highbd_smooth_predictor_8x16_avx2;
        eb_aom_highbd_smooth_predictor_8x32 = eb_aom_highbd_smooth_predictor_8x32_c;
        if (flags & HAS_AVX2) eb_aom_highbd_smooth_predictor_8x32 = eb_aom_highbd_smooth_predictor_8x32_avx2;
        eb_aom_highbd_smooth_predictor_8x4 = eb_aom_highbd_smooth_predictor_8x4_c;
        if (flags & HAS_AVX2) eb_aom_highbd_smooth_predictor_8x4 = eb_aom_highbd_smooth_predictor_8x4_avx2;
        eb_aom_highbd_smooth_predictor_8x8 = eb_aom_highbd_smooth_predictor_8x8_c;
        if (flags & HAS_AVX2) eb_aom_highbd_smooth_predictor_8x8 = eb_aom_highbd_smooth_predictor_8x8_avx2;

        //aom_highbd_smooth_h_predictor
        eb_aom_highbd_smooth_h_predictor_16x16 = eb_aom_highbd_smooth_h_predictor_16x16_c;
        if (flags & HAS_AVX2) eb_aom_highbd_smooth_h_predictor_16x16 = eb_aom_highbd_smooth_h_predictor_16x16_avx2;
        eb_aom_highbd_smooth_h_predictor_16x32 = eb_aom_highbd_smooth_h_predictor_16x32_c;
        if (flags & HAS_AVX2) eb_aom_highbd_smooth_h_predictor_16x32 = eb_aom_highbd_smooth_h_predictor_16x32_avx2;
        eb_aom_highbd_smooth_h_predictor_16x4 = eb_aom_highbd_smooth_h_predictor_16x4_c;
        if (flags & HAS_AVX2) eb_aom_highbd_smooth_h_predictor_16x4 = eb_aom_highbd_smooth_h_predictor_16x4_avx2;
        eb_aom_highbd_smooth_h_predictor_16x64 = eb_aom_highbd_smooth_h_predictor_16x64_c;
        if (flags & HAS_AVX2) eb_aom_highbd_smooth_h_predictor_16x64 = eb_aom_highbd_smooth_h_predictor_16x64_avx2;
        eb_aom_highbd_smooth_h_predictor_16x8 = eb_aom_highbd_smooth_h_predictor_16x8_c;
        if (flags & HAS_AVX2) eb_aom_highbd_smooth_h_predictor_16x8 = eb_aom_highbd_smooth_h_predictor_16x8_avx2;
        eb_aom_highbd_smooth_h_predictor_2x2 = eb_aom_highbd_smooth_h_predictor_2x2_c;
        eb_aom_highbd_smooth_h_predictor_32x16 = eb_aom_highbd_smooth_h_predictor_32x16_c;
        if (flags & HAS_AVX2) eb_aom_highbd_smooth_h_predictor_32x16 = eb_aom_highbd_smooth_h_predictor_32x16_avx2;
        eb_aom_highbd_smooth_h_predictor_32x32 = eb_aom_highbd_smooth_h_predictor_32x32_c;
        if (flags & HAS_AVX2) eb_aom_highbd_smooth_h_predictor_32x32 = eb_aom_highbd_smooth_h_predictor_32x32_avx2;
        eb_aom_highbd_smooth_h_predictor_32x64 = eb_aom_highbd_smooth_h_predictor_32x64_c;
        if (flags & HAS_AVX2) eb_aom_highbd_smooth_h_predictor_32x64 = eb_aom_highbd_smooth_h_predictor_32x64_avx2;
        eb_aom_highbd_smooth_h_predictor_32x8 = eb_aom_highbd_smooth_h_predictor_32x8_c;
        if (flags & HAS_AVX2) eb_aom_highbd_smooth_h_predictor_32x8 = eb_aom_highbd_smooth_h_predictor_32x8_avx2;
        eb_aom_highbd_smooth_h_predictor_4x16 = eb_aom_highbd_smooth_h_predictor_4x16_c;
        if (flags & HAS_SSSE3) eb_aom_highbd_smooth_h_predictor_4x16 = eb_aom_highbd_smooth_h_predictor_4x16_ssse3;
        eb_aom_highbd_smooth_h_predictor_4x4 = eb_aom_highbd_smooth_h_predictor_4x4_c;
        if (flags & HAS_SSSE3) eb_aom_highbd_smooth_h_predictor_4x4 = eb_aom_highbd_smooth_h_predictor_4x4_ssse3;
        eb_aom_highbd_smooth_h_predictor_4x8 = eb_aom_highbd_smooth_h_predictor_4x8_c;
        if (flags & HAS_SSSE3) eb_aom_highbd_smooth_h_predictor_4x8 = eb_aom_highbd_smooth_h_predictor_4x8_ssse3;
        eb_aom_highbd_smooth_h_predictor_64x16 = eb_aom_highbd_smooth_h_predictor_64x16_c;
        if (flags & HAS_AVX2) eb_aom_highbd_smooth_h_predictor_64x16 = eb_aom_highbd_smooth_h_predictor_64x16_avx2;
        eb_aom_highbd_smooth_h_predictor_64x32 = eb_aom_highbd_smooth_h_predictor_64x32_c;
        if (flags & HAS_AVX2) eb_aom_highbd_smooth_h_predictor_64x32 = eb_aom_highbd_smooth_h_predictor_64x32_avx2;
        eb_aom_highbd_smooth_h_predictor_64x64 = eb_aom_highbd_smooth_h_predictor_64x64_c;
        if (flags & HAS_AVX2) eb_aom_highbd_smooth_h_predictor_64x64 = eb_aom_highbd_smooth_h_predictor_64x64_avx2;
        eb_aom_highbd_smooth_h_predictor_8x16 = eb_aom_highbd_smooth_h_predictor_8x16_c;
        if (flags & HAS_AVX2) eb_aom_highbd_smooth_h_predictor_8x16 = eb_aom_highbd_smooth_h_predictor_8x16_avx2;
        eb_aom_highbd_smooth_h_predictor_8x32 = eb_aom_highbd_smooth_h_predictor_8x32_c;
        if (flags & HAS_AVX2) eb_aom_highbd_smooth_h_predictor_8x32 = eb_aom_highbd_smooth_h_predictor_8x32_avx2;
        eb_aom_highbd_smooth_h_predictor_8x4 = eb_aom_highbd_smooth_h_predictor_8x4_c;
        if (flags & HAS_AVX2) eb_aom_highbd_smooth_h_predictor_8x4 = eb_aom_highbd_smooth_h_predictor_8x4_avx2;
        eb_aom_highbd_smooth_h_predictor_8x8 = eb_aom_highbd_smooth_h_predictor_8x8_c;
        if (flags & HAS_AVX2) eb_aom_highbd_smooth_h_predictor_8x8 = eb_aom_highbd_smooth_h_predictor_8x8_avx2;

        //aom_highbd_dc_128_predictor
        eb_aom_highbd_dc_128_predictor_16x16 = eb_aom_highbd_dc_128_predictor_16x16_c;
        if (flags & HAS_AVX2) eb_aom_highbd_dc_128_predictor_16x16 = eb_aom_highbd_dc_128_predictor_16x16_avx2;
        eb_aom_highbd_dc_128_predictor_16x32 = eb_aom_highbd_dc_128_predictor_16x32_c;
        if (flags & HAS_AVX2) eb_aom_highbd_dc_128_predictor_16x32 = eb_aom_highbd_dc_128_predictor_16x32_avx2;
        eb_aom_highbd_dc_128_predictor_16x4 = eb_aom_highbd_dc_128_predictor_16x4_c;
        if (flags & HAS_AVX2) eb_aom_highbd_dc_128_predictor_16x4 = eb_aom_highbd_dc_128_predictor_16x4_avx2;
        eb_aom_highbd_dc_128_predictor_16x64 = eb_aom_highbd_dc_128_predictor_16x64_c;
        if (flags & HAS_AVX2) eb_aom_highbd_dc_128_predictor_16x64 = eb_aom_highbd_dc_128_predictor_16x64_avx2;
        eb_aom_highbd_dc_128_predictor_16x8 = eb_aom_highbd_dc_128_predictor_16x8_c;
        if (flags & HAS_AVX2) eb_aom_highbd_dc_128_predictor_16x8 = eb_aom_highbd_dc_128_predictor_16x8_avx2;
        eb_aom_highbd_dc_128_predictor_2x2 = eb_aom_highbd_dc_128_predictor_2x2_c;
        eb_aom_highbd_dc_128_predictor_32x16 = eb_aom_highbd_dc_128_predictor_32x16_c;
        if (flags & HAS_AVX2) eb_aom_highbd_dc_128_predictor_32x16 = eb_aom_highbd_dc_128_predictor_32x16_avx2;
        eb_aom_highbd_dc_128_predictor_32x32 = eb_aom_highbd_dc_128_predictor_32x32_c;
        if (flags & HAS_AVX2) eb_aom_highbd_dc_128_predictor_32x32 = eb_aom_highbd_dc_128_predictor_32x32_avx2;
        eb_aom_highbd_dc_128_predictor_32x64 = eb_aom_highbd_dc_128_predictor_32x64_c;
        if (flags & HAS_AVX2) eb_aom_highbd_dc_128_predictor_32x64 = eb_aom_highbd_dc_128_predictor_32x64_avx2;
        eb_aom_highbd_dc_128_predictor_32x8 = eb_aom_highbd_dc_128_predictor_32x8_c;
        if (flags & HAS_AVX2) eb_aom_highbd_dc_128_predictor_32x8 = eb_aom_highbd_dc_128_predictor_32x8_avx2;
        eb_aom_highbd_dc_128_predictor_4x16 = eb_aom_highbd_dc_128_predictor_4x16_c;
        if (flags & HAS_SSE2) eb_aom_highbd_dc_128_predictor_4x16 = eb_aom_highbd_dc_128_predictor_4x16_sse2;
        eb_aom_highbd_dc_128_predictor_4x4 = eb_aom_highbd_dc_128_predictor_4x4_c;
        if (flags & HAS_SSE2) eb_aom_highbd_dc_128_predictor_4x4 = eb_aom_highbd_dc_128_predictor_4x4_sse2;
        eb_aom_highbd_dc_128_predictor_4x8 = eb_aom_highbd_dc_128_predictor_4x8_c;
        if (flags & HAS_SSE2) eb_aom_highbd_dc_128_predictor_4x8 = eb_aom_highbd_dc_128_predictor_4x8_sse2;
        eb_aom_highbd_dc_128_predictor_8x32 = eb_aom_highbd_dc_128_predictor_8x32_c;
        if (flags & HAS_SSE2) eb_aom_highbd_dc_128_predictor_8x32 = eb_aom_highbd_dc_128_predictor_8x32_sse2;
        eb_aom_highbd_dc_128_predictor_64x16 = eb_aom_highbd_dc_128_predictor_64x16_c;
        if (flags & HAS_AVX2) eb_aom_highbd_dc_128_predictor_64x16 = eb_aom_highbd_dc_128_predictor_64x16_avx2;
        eb_aom_highbd_dc_128_predictor_64x32 = eb_aom_highbd_dc_128_predictor_64x32_c;
        if (flags & HAS_AVX2) eb_aom_highbd_dc_128_predictor_64x32 = eb_aom_highbd_dc_128_predictor_64x32_avx2;
        eb_aom_highbd_dc_128_predictor_64x64 = eb_aom_highbd_dc_128_predictor_64x64_c;
        if (flags & HAS_AVX2) eb_aom_highbd_dc_128_predictor_64x64 = eb_aom_highbd_dc_128_predictor_64x64_avx2;
        eb_aom_highbd_dc_128_predictor_8x16 = eb_aom_highbd_dc_128_predictor_8x16_c;
        if (flags & HAS_SSE2) eb_aom_highbd_dc_128_predictor_8x16 = eb_aom_highbd_dc_128_predictor_8x16_sse2;
        eb_aom_highbd_dc_128_predictor_8x4 = eb_aom_highbd_dc_128_predictor_8x4_c;
        if (flags & HAS_SSE2) eb_aom_highbd_dc_128_predictor_8x4 = eb_aom_highbd_dc_128_predictor_8x4_sse2;
        eb_aom_highbd_dc_128_predictor_8x8 = eb_aom_highbd_dc_128_predictor_8x8_c;
        if (flags & HAS_SSE2) eb_aom_highbd_dc_128_predictor_8x8 = eb_aom_highbd_dc_128_predictor_8x8_sse2;

        //aom_highbd_dc_left_predictor
        eb_aom_highbd_dc_left_predictor_16x16 = eb_aom_highbd_dc_left_predictor_16x16_c;
        if (flags & HAS_AVX2) eb_aom_highbd_dc_left_predictor_16x16 = eb_aom_highbd_dc_left_predictor_16x16_avx2;
        eb_aom_highbd_dc_left_predictor_16x32 = eb_aom_highbd_dc_left_predictor_16x32_c;
        if (flags & HAS_AVX2) eb_aom_highbd_dc_left_predictor_16x32 = eb_aom_highbd_dc_left_predictor_16x32_avx2;
        eb_aom_highbd_dc_left_predictor_16x4 = eb_aom_highbd_dc_left_predictor_16x4_c;
        if (flags & HAS_AVX2) eb_aom_highbd_dc_left_predictor_16x4 = eb_aom_highbd_dc_left_predictor_16x4_avx2;
        eb_aom_highbd_dc_left_predictor_16x64 = eb_aom_highbd_dc_left_predictor_16x64_c;
        if (flags & HAS_AVX2) eb_aom_highbd_dc_left_predictor_16x64 = eb_aom_highbd_dc_left_predictor_16x64_avx2;
        eb_aom_highbd_dc_left_predictor_16x8 = eb_aom_highbd_dc_left_predictor_16x8_c;
        if (flags & HAS_AVX2) eb_aom_highbd_dc_left_predictor_16x8 = eb_aom_highbd_dc_left_predictor_16x8_avx2;
        eb_aom_highbd_dc_left_predictor_2x2 = eb_aom_highbd_dc_left_predictor_2x2_c;
        eb_aom_highbd_dc_left_predictor_32x16 = eb_aom_highbd_dc_left_predictor_32x16_c;
        eb_aom_highbd_dc_left_predictor_32x32 = eb_aom_highbd_dc_left_predictor_32x32_c;
        eb_aom_highbd_dc_left_predictor_32x64 = eb_aom_highbd_dc_left_predictor_32x64_c;
        eb_aom_highbd_dc_left_predictor_32x8 = eb_aom_highbd_dc_left_predictor_32x8_c;
        eb_aom_highbd_dc_left_predictor_4x16 = eb_aom_highbd_dc_left_predictor_4x16_c;
        if (flags & HAS_SSE2) eb_aom_highbd_dc_left_predictor_4x16 = eb_aom_highbd_dc_left_predictor_4x16_sse2;
        eb_aom_highbd_dc_left_predictor_4x4 = eb_aom_highbd_dc_left_predictor_4x4_c;
        if (flags & HAS_SSE2) eb_aom_highbd_dc_left_predictor_4x4 = eb_aom_highbd_dc_left_predictor_4x4_sse2;
        eb_aom_highbd_dc_left_predictor_4x8 = eb_aom_highbd_dc_left_predictor_4x8_c;
        if (flags & HAS_SSE2) eb_aom_highbd_dc_left_predictor_4x8 = eb_aom_highbd_dc_left_predictor_4x8_sse2;
        eb_aom_highbd_dc_left_predictor_8x32 = eb_aom_highbd_dc_left_predictor_8x32_c;
        if (flags & HAS_SSE2) eb_aom_highbd_dc_left_predictor_8x32 = eb_aom_highbd_dc_left_predictor_8x32_sse2;
        eb_aom_highbd_dc_left_predictor_64x16 = eb_aom_highbd_dc_left_predictor_64x16_c;
        eb_aom_highbd_dc_left_predictor_64x32 = eb_aom_highbd_dc_left_predictor_64x32_c;
        eb_aom_highbd_dc_left_predictor_64x64 = eb_aom_highbd_dc_left_predictor_64x64_c;
        eb_aom_highbd_dc_left_predictor_8x16 = eb_aom_highbd_dc_left_predictor_8x16_c;
        if (flags & HAS_SSE2) eb_aom_highbd_dc_left_predictor_8x16 = eb_aom_highbd_dc_left_predictor_8x16_sse2;
        eb_aom_highbd_dc_left_predictor_8x4 = eb_aom_highbd_dc_left_predictor_8x4_c;
        if (flags & HAS_SSE2) eb_aom_highbd_dc_left_predictor_8x4 = eb_aom_highbd_dc_left_predictor_8x4_sse2;
        eb_aom_highbd_dc_left_predictor_8x8 = eb_aom_highbd_dc_left_predictor_8x8_c;
        if (flags & HAS_SSE2) eb_aom_highbd_dc_left_predictor_8x8 = eb_aom_highbd_dc_left_predictor_8x8_sse2;

#ifndef NON_AVX512_SUPPORT
        if (flags & HAS_AVX2) eb_aom_highbd_dc_left_predictor_32x8 = aom_highbd_dc_left_predictor_32x8_avx512;
        if (flags & HAS_AVX2) eb_aom_highbd_dc_left_predictor_32x16 = aom_highbd_dc_left_predictor_32x16_avx512;
        if (flags & HAS_AVX2) eb_aom_highbd_dc_left_predictor_32x32 = aom_highbd_dc_left_predictor_32x32_avx512;
        if (flags & HAS_AVX2) eb_aom_highbd_dc_left_predictor_32x64 = aom_highbd_dc_left_predictor_32x64_avx512;
        if (flags & HAS_AVX2) eb_aom_highbd_dc_left_predictor_64x16 = aom_highbd_dc_left_predictor_64x16_avx512;
        if (flags & HAS_AVX2) eb_aom_highbd_dc_left_predictor_64x32 = aom_highbd_dc_left_predictor_64x32_avx512;
        if (flags & HAS_AVX2) eb_aom_highbd_dc_left_predictor_64x64 = aom_highbd_dc_left_predictor_64x64_avx512;
#else
        if (flags & HAS_AVX2) eb_aom_highbd_dc_left_predictor_32x8 = eb_aom_highbd_dc_left_predictor_32x8_avx2;
        if (flags & HAS_AVX2) eb_aom_highbd_dc_left_predictor_32x16 = eb_aom_highbd_dc_left_predictor_32x16_avx2;
        if (flags & HAS_AVX2) eb_aom_highbd_dc_left_predictor_32x32 = eb_aom_highbd_dc_left_predictor_32x32_avx2;
        if (flags & HAS_AVX2) eb_aom_highbd_dc_left_predictor_32x64 = eb_aom_highbd_dc_left_predictor_32x64_avx2;
        if (flags & HAS_AVX2) eb_aom_highbd_dc_left_predictor_64x16 = eb_aom_highbd_dc_left_predictor_64x16_avx2;
        if (flags & HAS_AVX2) eb_aom_highbd_dc_left_predictor_64x32 = eb_aom_highbd_dc_left_predictor_64x32_avx2;
        if (flags & HAS_AVX2) eb_aom_highbd_dc_left_predictor_64x64 = eb_aom_highbd_dc_left_predictor_64x64_avx2;
#endif // !NON_AVX512_SUPPORT

        eb_aom_highbd_dc_predictor_16x16 = eb_aom_highbd_dc_predictor_16x16_c;
        if (flags & HAS_AVX2) eb_aom_highbd_dc_predictor_16x16 = eb_aom_highbd_dc_predictor_16x16_avx2;
        eb_aom_highbd_dc_predictor_16x32 = eb_aom_highbd_dc_predictor_16x32_c;
        if (flags & HAS_AVX2) eb_aom_highbd_dc_predictor_16x32 = eb_aom_highbd_dc_predictor_16x32_avx2;
        eb_aom_highbd_dc_predictor_16x4 = eb_aom_highbd_dc_predictor_16x4_c;
        if (flags & HAS_AVX2) eb_aom_highbd_dc_predictor_16x4 = eb_aom_highbd_dc_predictor_16x4_avx2;
        eb_aom_highbd_dc_predictor_16x64 = eb_aom_highbd_dc_predictor_16x64_c;
        if (flags & HAS_AVX2) eb_aom_highbd_dc_predictor_16x64 = eb_aom_highbd_dc_predictor_16x64_avx2;
        eb_aom_highbd_dc_predictor_16x8 = eb_aom_highbd_dc_predictor_16x8_c;
        if (flags & HAS_AVX2) eb_aom_highbd_dc_predictor_16x8 = eb_aom_highbd_dc_predictor_16x8_avx2;
        eb_aom_highbd_dc_predictor_2x2 = eb_aom_highbd_dc_predictor_2x2_c;
        eb_aom_highbd_dc_predictor_32x16 = eb_aom_highbd_dc_predictor_32x16_c;
        eb_aom_highbd_dc_predictor_32x32 = eb_aom_highbd_dc_predictor_32x32_c;
        eb_aom_highbd_dc_predictor_32x64 = eb_aom_highbd_dc_predictor_32x64_c;
        eb_aom_highbd_dc_predictor_32x8 = eb_aom_highbd_dc_predictor_32x8_c;
        eb_aom_highbd_dc_predictor_4x16 = eb_aom_highbd_dc_predictor_4x16_c;
        if (flags & HAS_SSE2) eb_aom_highbd_dc_predictor_4x16 = eb_aom_highbd_dc_predictor_4x16_sse2;
        eb_aom_highbd_dc_predictor_4x4 = eb_aom_highbd_dc_predictor_4x4_c;
        if (flags & HAS_SSE2) eb_aom_highbd_dc_predictor_4x4 = eb_aom_highbd_dc_predictor_4x4_sse2;
        eb_aom_highbd_dc_predictor_4x8 = eb_aom_highbd_dc_predictor_4x8_c;
        if (flags & HAS_SSE2) eb_aom_highbd_dc_predictor_4x8 = eb_aom_highbd_dc_predictor_4x8_sse2;
        eb_aom_highbd_dc_predictor_64x16 = eb_aom_highbd_dc_predictor_64x16_c;
        eb_aom_highbd_dc_predictor_64x32 = eb_aom_highbd_dc_predictor_64x32_c;
        eb_aom_highbd_dc_predictor_64x64 = eb_aom_highbd_dc_predictor_64x64_c;
        eb_aom_highbd_dc_predictor_8x16 = eb_aom_highbd_dc_predictor_8x16_c;
        if (flags & HAS_SSE2) eb_aom_highbd_dc_predictor_8x16 = eb_aom_highbd_dc_predictor_8x16_sse2;
        eb_aom_highbd_dc_predictor_8x4 = eb_aom_highbd_dc_predictor_8x4_c;
        if (flags & HAS_SSE2) eb_aom_highbd_dc_predictor_8x4 = eb_aom_highbd_dc_predictor_8x4_sse2;
        eb_aom_highbd_dc_predictor_8x8 = eb_aom_highbd_dc_predictor_8x8_c;
        if (flags & HAS_SSE2) eb_aom_highbd_dc_predictor_8x8 = eb_aom_highbd_dc_predictor_8x8_sse2;
        eb_aom_highbd_dc_predictor_8x32 = eb_aom_highbd_dc_predictor_8x32_c;
        if (flags & HAS_SSE2) eb_aom_highbd_dc_predictor_8x32 = eb_aom_highbd_dc_predictor_8x32_sse2;

#ifndef NON_AVX512_SUPPORT
        if (flags & HAS_AVX2) eb_aom_highbd_dc_predictor_32x8 = aom_highbd_dc_predictor_32x8_avx512;
        if (flags & HAS_AVX2) eb_aom_highbd_dc_predictor_32x16 = aom_highbd_dc_predictor_32x16_avx512;
        if (flags & HAS_AVX2) eb_aom_highbd_dc_predictor_32x32 = aom_highbd_dc_predictor_32x32_avx512;
        if (flags & HAS_AVX2) eb_aom_highbd_dc_predictor_32x64 = aom_highbd_dc_predictor_32x64_avx512;
        if (flags & HAS_AVX2) eb_aom_highbd_dc_predictor_64x16 = aom_highbd_dc_predictor_64x16_avx512;
        if (flags & HAS_AVX2) eb_aom_highbd_dc_predictor_64x32 = aom_highbd_dc_predictor_64x32_avx512;
        if (flags & HAS_AVX2) eb_aom_highbd_dc_predictor_64x64 = aom_highbd_dc_predictor_64x64_avx512;
#else
        if (flags & HAS_AVX2) eb_aom_highbd_dc_predictor_32x8 = eb_aom_highbd_dc_predictor_32x8_avx2;
        if (flags & HAS_AVX2) eb_aom_highbd_dc_predictor_32x16 = eb_aom_highbd_dc_predictor_32x16_avx2;
        if (flags & HAS_AVX2) eb_aom_highbd_dc_predictor_32x32 = eb_aom_highbd_dc_predictor_32x32_avx2;
        if (flags & HAS_AVX2) eb_aom_highbd_dc_predictor_32x64 = eb_aom_highbd_dc_predictor_32x64_avx2;
        if (flags & HAS_AVX2) eb_aom_highbd_dc_predictor_64x16 = eb_aom_highbd_dc_predictor_64x16_avx2;
        if (flags & HAS_AVX2) eb_aom_highbd_dc_predictor_64x32 = eb_aom_highbd_dc_predictor_64x32_avx2;
        if (flags & HAS_AVX2) eb_aom_highbd_dc_predictor_64x64 = eb_aom_highbd_dc_predictor_64x64_avx2;
#endif // !NON_AVX512_SUPPORT
        //aom_highbd_dc_top_predictor
        eb_aom_highbd_dc_top_predictor_16x16 = eb_aom_highbd_dc_top_predictor_16x16_c;
        if (flags & HAS_AVX2) eb_aom_highbd_dc_top_predictor_16x16 = eb_aom_highbd_dc_top_predictor_16x16_avx2;
        eb_aom_highbd_dc_top_predictor_16x32 = eb_aom_highbd_dc_top_predictor_16x32_c;
        if (flags & HAS_AVX2) eb_aom_highbd_dc_top_predictor_16x32 = eb_aom_highbd_dc_top_predictor_16x32_avx2;
        eb_aom_highbd_dc_top_predictor_16x4 = eb_aom_highbd_dc_top_predictor_16x4_c;
        if (flags & HAS_AVX2) eb_aom_highbd_dc_top_predictor_16x4 = eb_aom_highbd_dc_top_predictor_16x4_avx2;
        eb_aom_highbd_dc_top_predictor_16x64 = eb_aom_highbd_dc_top_predictor_16x64_c;
        if (flags & HAS_AVX2) eb_aom_highbd_dc_top_predictor_16x64 = eb_aom_highbd_dc_top_predictor_16x64_avx2;
        eb_aom_highbd_dc_top_predictor_16x8 = eb_aom_highbd_dc_top_predictor_16x8_c;
        if (flags & HAS_AVX2) eb_aom_highbd_dc_top_predictor_16x8 = eb_aom_highbd_dc_top_predictor_16x8_avx2;
        eb_aom_highbd_dc_top_predictor_2x2 = eb_aom_highbd_dc_top_predictor_2x2_c;
        eb_aom_highbd_dc_top_predictor_32x16 = eb_aom_highbd_dc_top_predictor_32x16_c;
        eb_aom_highbd_dc_top_predictor_32x32 = eb_aom_highbd_dc_top_predictor_32x32_c;
        eb_aom_highbd_dc_top_predictor_32x64 = eb_aom_highbd_dc_top_predictor_32x64_c;
        eb_aom_highbd_dc_top_predictor_32x8 = eb_aom_highbd_dc_top_predictor_32x8_c;
        eb_aom_highbd_dc_top_predictor_4x16 = eb_aom_highbd_dc_top_predictor_4x16_c;
        if (flags & HAS_SSE2) eb_aom_highbd_dc_top_predictor_4x16 = eb_aom_highbd_dc_top_predictor_4x16_sse2;
        eb_aom_highbd_dc_top_predictor_4x4 = eb_aom_highbd_dc_top_predictor_4x4_c;
        if (flags & HAS_SSE2) eb_aom_highbd_dc_top_predictor_4x4 = eb_aom_highbd_dc_top_predictor_4x4_sse2;
        eb_aom_highbd_dc_top_predictor_4x8 = eb_aom_highbd_dc_top_predictor_4x8_c;
        if (flags & HAS_SSE2) eb_aom_highbd_dc_top_predictor_4x8 = eb_aom_highbd_dc_top_predictor_4x8_sse2;
        eb_aom_highbd_dc_top_predictor_64x16 = eb_aom_highbd_dc_top_predictor_64x16_c;
        eb_aom_highbd_dc_top_predictor_64x32 = eb_aom_highbd_dc_top_predictor_64x32_c;
        eb_aom_highbd_dc_top_predictor_64x64 = eb_aom_highbd_dc_top_predictor_64x64_c;
        eb_aom_highbd_dc_top_predictor_8x16 = eb_aom_highbd_dc_top_predictor_8x16_c;
        if (flags & HAS_SSE2) eb_aom_highbd_dc_top_predictor_8x16 = eb_aom_highbd_dc_top_predictor_8x16_sse2;
        if (flags & HAS_SSE2) eb_aom_highbd_dc_top_predictor_8x32 = eb_aom_highbd_dc_top_predictor_8x32_c;
        eb_aom_highbd_dc_top_predictor_8x4 = eb_aom_highbd_dc_top_predictor_8x4_c;
        if (flags & HAS_SSE2) eb_aom_highbd_dc_top_predictor_8x4 = eb_aom_highbd_dc_top_predictor_8x4_sse2;
        eb_aom_highbd_dc_top_predictor_8x8 = eb_aom_highbd_dc_top_predictor_8x8_c;
        if (flags & HAS_SSE2) eb_aom_highbd_dc_top_predictor_8x8 = eb_aom_highbd_dc_top_predictor_8x8_sse2;

#ifndef NON_AVX512_SUPPORT
        if (flags & HAS_AVX2) eb_aom_highbd_dc_top_predictor_32x8 = aom_highbd_dc_top_predictor_32x8_avx512;
        if (flags & HAS_AVX2) eb_aom_highbd_dc_top_predictor_32x16 = aom_highbd_dc_top_predictor_32x16_avx512;
        if (flags & HAS_AVX2) eb_aom_highbd_dc_top_predictor_32x32 = aom_highbd_dc_top_predictor_32x32_avx512;
        if (flags & HAS_AVX2) eb_aom_highbd_dc_top_predictor_32x64 = aom_highbd_dc_top_predictor_32x64_avx512;
        if (flags & HAS_AVX2) eb_aom_highbd_dc_top_predictor_64x16 = aom_highbd_dc_top_predictor_64x16_avx512;
        if (flags & HAS_AVX2) eb_aom_highbd_dc_top_predictor_64x32 = aom_highbd_dc_top_predictor_64x32_avx512;
        if (flags & HAS_AVX2) eb_aom_highbd_dc_top_predictor_64x64 = aom_highbd_dc_top_predictor_64x64_avx512;
#else
        if (flags & HAS_AVX2) eb_aom_highbd_dc_top_predictor_32x8 = eb_aom_highbd_dc_top_predictor_32x8_avx2;
        if (flags & HAS_AVX2) eb_aom_highbd_dc_top_predictor_32x16 = eb_aom_highbd_dc_top_predictor_32x16_avx2;
        if (flags & HAS_AVX2) eb_aom_highbd_dc_top_predictor_32x32 = eb_aom_highbd_dc_top_predictor_32x32_avx2;
        if (flags & HAS_AVX2) eb_aom_highbd_dc_top_predictor_32x64 = eb_aom_highbd_dc_top_predictor_32x64_avx2;
        if (flags & HAS_AVX2) eb_aom_highbd_dc_top_predictor_64x16 = eb_aom_highbd_dc_top_predictor_64x16_avx2;
        if (flags & HAS_AVX2) eb_aom_highbd_dc_top_predictor_64x32 = eb_aom_highbd_dc_top_predictor_64x32_avx2;
        if (flags & HAS_AVX2) eb_aom_highbd_dc_top_predictor_64x64 = eb_aom_highbd_dc_top_predictor_64x64_avx2;
#endif
        // eb_aom_highbd_h_predictor
        eb_aom_highbd_h_predictor_16x4 = eb_aom_highbd_h_predictor_16x4_c;
        if (flags & HAS_AVX2) eb_aom_highbd_h_predictor_16x4 = eb_aom_highbd_h_predictor_16x4_avx2;
        eb_aom_highbd_h_predictor_16x64 = eb_aom_highbd_h_predictor_16x64_c;
        if (flags & HAS_AVX2) eb_aom_highbd_h_predictor_16x64 = eb_aom_highbd_h_predictor_16x64_avx2;
        eb_aom_highbd_h_predictor_16x8 = eb_aom_highbd_h_predictor_16x8_c;
        if (flags & HAS_SSE2) eb_aom_highbd_h_predictor_16x8 = eb_aom_highbd_h_predictor_16x8_sse2;
        eb_aom_highbd_h_predictor_2x2 = eb_aom_highbd_h_predictor_2x2_c;
        eb_aom_highbd_h_predictor_32x16 = eb_aom_highbd_h_predictor_32x16_c;
        if (flags & HAS_SSE2) eb_aom_highbd_h_predictor_32x16 = eb_aom_highbd_h_predictor_32x16_sse2;
        eb_aom_highbd_h_predictor_32x32 = eb_aom_highbd_h_predictor_32x32_c;
        if (flags & HAS_SSE2) eb_aom_highbd_h_predictor_32x32 = eb_aom_highbd_h_predictor_32x32_sse2;
        eb_aom_highbd_h_predictor_32x64 = eb_aom_highbd_h_predictor_32x64_c;
        if (flags & HAS_AVX2) eb_aom_highbd_h_predictor_32x64 = eb_aom_highbd_h_predictor_32x64_avx2;
        eb_aom_highbd_h_predictor_32x8 = eb_aom_highbd_h_predictor_32x8_c;
        if (flags & HAS_AVX2) eb_aom_highbd_h_predictor_32x8 = eb_aom_highbd_h_predictor_32x8_avx2;
        eb_aom_highbd_h_predictor_4x16 = eb_aom_highbd_h_predictor_4x16_c;
        if (flags & HAS_SSE2) eb_aom_highbd_h_predictor_4x16 = eb_aom_highbd_h_predictor_4x16_sse2;
        eb_aom_highbd_h_predictor_4x4 = eb_aom_highbd_h_predictor_4x4_c;
        if (flags & HAS_SSE2) eb_aom_highbd_h_predictor_4x4 = eb_aom_highbd_h_predictor_4x4_sse2;
        eb_aom_highbd_h_predictor_4x8 = eb_aom_highbd_h_predictor_4x8_c;
        if (flags & HAS_SSE2) eb_aom_highbd_h_predictor_4x8 = eb_aom_highbd_h_predictor_4x8_sse2;
        eb_aom_highbd_h_predictor_64x16 = eb_aom_highbd_h_predictor_64x16_c;
        if (flags & HAS_AVX2) eb_aom_highbd_h_predictor_64x16 = eb_aom_highbd_h_predictor_64x16_avx2;
        eb_aom_highbd_h_predictor_64x32 = eb_aom_highbd_h_predictor_64x32_c;
        if (flags & HAS_AVX2) eb_aom_highbd_h_predictor_64x32 = eb_aom_highbd_h_predictor_64x32_avx2;
        eb_aom_highbd_h_predictor_8x32 = eb_aom_highbd_h_predictor_8x32_c;
        if (flags & HAS_SSE2) eb_aom_highbd_h_predictor_8x32 = eb_aom_highbd_h_predictor_8x32_sse2;
        eb_aom_highbd_h_predictor_64x64 = eb_aom_highbd_h_predictor_64x64_c;
        if (flags & HAS_AVX2) eb_aom_highbd_h_predictor_64x64 = eb_aom_highbd_h_predictor_64x64_avx2;
        eb_aom_highbd_h_predictor_8x16 = eb_aom_highbd_h_predictor_8x16_c;
        if (flags & HAS_SSE2) eb_aom_highbd_h_predictor_8x16 = eb_aom_highbd_h_predictor_8x16_sse2;
        eb_aom_highbd_h_predictor_8x4 = eb_aom_highbd_h_predictor_8x4_c;
        if (flags & HAS_SSE2) eb_aom_highbd_h_predictor_8x4 = eb_aom_highbd_h_predictor_8x4_sse2;
        eb_aom_highbd_h_predictor_8x8 = eb_aom_highbd_h_predictor_8x8_c;
        if (flags & HAS_SSE2) eb_aom_highbd_h_predictor_8x8 = eb_aom_highbd_h_predictor_8x8_sse2;
        eb_aom_highbd_h_predictor_16x16 = eb_aom_highbd_h_predictor_16x16_c;
        if (flags & HAS_SSE2) eb_aom_highbd_h_predictor_16x16 = eb_aom_highbd_h_predictor_16x16_sse2;
        eb_aom_highbd_h_predictor_16x32 = eb_aom_highbd_h_predictor_16x32_c;
        if (flags & HAS_SSE2) eb_aom_highbd_h_predictor_16x32 = eb_aom_highbd_h_predictor_16x32_sse2;

        eb_aom_fft2x2_float = eb_aom_fft2x2_float_c;
        eb_aom_fft4x4_float = eb_aom_fft4x4_float_c;
        if (flags & HAS_SSE2) eb_aom_fft4x4_float = eb_aom_fft4x4_float_sse2;
        eb_aom_fft16x16_float = eb_aom_fft16x16_float_c;
        if (flags & HAS_AVX2) eb_aom_fft16x16_float = eb_aom_fft16x16_float_avx2;
        eb_aom_fft32x32_float = eb_aom_fft32x32_float_c;
        if (flags & HAS_AVX2) eb_aom_fft32x32_float = eb_aom_fft32x32_float_avx2;
        eb_aom_fft8x8_float = eb_aom_fft8x8_float_c;
        if (flags & HAS_AVX2) eb_aom_fft8x8_float = eb_aom_fft8x8_float_avx2;

        /*if (flags & HAS_AVX2)*/ eb_aom_ifft16x16_float = eb_aom_ifft16x16_float_avx2;
        /*if (flags & HAS_AVX2)*/ eb_aom_ifft32x32_float = eb_aom_ifft32x32_float_avx2;
        /*if (flags & HAS_AVX2)*/ eb_aom_ifft8x8_float = eb_aom_ifft8x8_float_avx2;
        eb_aom_ifft2x2_float = eb_aom_ifft2x2_float_c;
        /*if (flags & HAS_SSE2)*/ eb_aom_ifft4x4_float = eb_aom_ifft4x4_float_sse2;
    }
#endif

#ifdef __cplusplus
}  // extern "C"
#endif

#endif
