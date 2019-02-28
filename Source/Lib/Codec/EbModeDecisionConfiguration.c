/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/


#include "EbModeDecisionConfiguration.h"
#include "EbLambdaRateTables.h"
#include "EbUtility.h"
#include "EbModeDecisionProcess.h"

static const uint32_t me2Nx2NOffset[4] = { 0, 1, 5, 21 };

/********************************************
 * Constants
 ********************************************/
#define ADD_CU_STOP_SPLIT             0   // Take into account & Stop Splitting
#define ADD_CU_CONTINUE_SPLIT         1   // Take into account & Continue Splitting
#define DO_NOT_ADD_CU_CONTINUE_SPLIT  2   // Do not take into account & Continue Splitting

#define DEPTH_64                      0   // Depth corresponding to the CU size
#define DEPTH_32                      1   // Depth corresponding to the CU size
#define DEPTH_16                      2   // Depth corresponding to the CU size
#define DEPTH_8                       3   // Depth corresponding to the CU size

static const uint8_t parentCuIndex[85] =
{
    0,
    0, 0, 0, 1, 2, 3, 5, 0, 1, 2, 3, 10, 0, 1, 2, 3, 15, 0, 1, 2, 3,
    21, 0, 0, 1, 2, 3, 5, 0, 1, 2, 3, 10, 0, 1, 2, 3, 15, 0, 1, 2, 3,
    42, 0, 0, 1, 2, 3, 5, 0, 1, 2, 3, 10, 0, 1, 2, 3, 15, 0, 1, 2, 3,
    36, 0, 0, 1, 2, 3, 5, 0, 1, 2, 3, 10, 0, 1, 2, 3, 15, 0, 1, 2, 3,
};

const uint8_t incrementalCount[85] = {

    //64x64
    0,
    //32x32
    4, 4,
    4, 4,
    //16x16
    0, 0, 0, 0,
    0, 4, 0, 4,
    0, 0, 0, 0,
    0, 4, 0, 4,
    //8x8
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 4, 0, 0, 0, 4,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 4, 0, 0, 0, 4

};

/*******************************************
mdcSetDepth : set depth to be tested
*******************************************/
#define REFINEMENT_P        0x01
#define REFINEMENT_Pp1      0x02
#define REFINEMENT_Pp2      0x04
#define REFINEMENT_Pp3      0x08
#define REFINEMENT_Pm1      0x10
#define REFINEMENT_Pm2      0x20
#define REFINEMENT_Pm3      0x40



EbErrorType MdcRefinement(
    MdcpLocalCodingUnit                   *local_cu_array,
    uint32_t                                  cu_index,
    uint32_t                                  depth,
    uint8_t                                   refinementLevel,
    uint8_t                                   lowestLevel)
{
    EbErrorType return_error = EB_ErrorNone;


    if (refinementLevel & REFINEMENT_P) {
        if (lowestLevel == REFINEMENT_P) {
            local_cu_array[cu_index].stopSplit = EB_TRUE;
        }

    }
    else {
        local_cu_array[cu_index].slectedCu = EB_FALSE;
    }

    if (refinementLevel & REFINEMENT_Pp1) {

        if (depth < 3 && cu_index < 81) {
            local_cu_array[cu_index + 1].slectedCu = EB_TRUE;
            local_cu_array[cu_index + 1 + depth_offset[depth + 1]].slectedCu = EB_TRUE;
            local_cu_array[cu_index + 1 + 2 * depth_offset[depth + 1]].slectedCu = EB_TRUE;
            local_cu_array[cu_index + 1 + 3 * depth_offset[depth + 1]].slectedCu = EB_TRUE;
        }
        if (lowestLevel == REFINEMENT_Pp1) {
            if (depth < 3 && cu_index < 81) {
                local_cu_array[cu_index + 1].stopSplit = EB_TRUE;
                local_cu_array[cu_index + 1 + depth_offset[depth + 1]].stopSplit = EB_TRUE;
                local_cu_array[cu_index + 1 + 2 * depth_offset[depth + 1]].stopSplit = EB_TRUE;
                local_cu_array[cu_index + 1 + 3 * depth_offset[depth + 1]].stopSplit = EB_TRUE;
            }
        }
    }

    if (refinementLevel & REFINEMENT_Pp2) {
        if (depth < 2 && cu_index < 65) {
            local_cu_array[cu_index + 1 + 1].slectedCu = EB_TRUE;
            local_cu_array[cu_index + 1 + 1 + depth_offset[depth + 2]].slectedCu = EB_TRUE;
            local_cu_array[cu_index + 1 + 1 + 2 * depth_offset[depth + 2]].slectedCu = EB_TRUE;
            local_cu_array[cu_index + 1 + 1 + 3 * depth_offset[depth + 2]].slectedCu = EB_TRUE;

            local_cu_array[cu_index + 1 + depth_offset[depth + 1] + 1].slectedCu = EB_TRUE;
            local_cu_array[cu_index + 1 + depth_offset[depth + 1] + 1 + depth_offset[depth + 2]].slectedCu = EB_TRUE;
            local_cu_array[cu_index + 1 + depth_offset[depth + 1] + 1 + 2 * depth_offset[depth + 2]].slectedCu = EB_TRUE;
            local_cu_array[cu_index + 1 + depth_offset[depth + 1] + 1 + 3 * depth_offset[depth + 2]].slectedCu = EB_TRUE;

            local_cu_array[cu_index + 1 + 2 * depth_offset[depth + 1] + 1].slectedCu = EB_TRUE;
            local_cu_array[cu_index + 1 + 2 * depth_offset[depth + 1] + 1 + depth_offset[depth + 2]].slectedCu = EB_TRUE;
            local_cu_array[cu_index + 1 + 2 * depth_offset[depth + 1] + 1 + 2 * depth_offset[depth + 2]].slectedCu = EB_TRUE;
            local_cu_array[cu_index + 1 + 2 * depth_offset[depth + 1] + 1 + 3 * depth_offset[depth + 2]].slectedCu = EB_TRUE;

            local_cu_array[cu_index + 1 + 3 * depth_offset[depth + 1] + 1].slectedCu = EB_TRUE;
            local_cu_array[cu_index + 1 + 3 * depth_offset[depth + 1] + 1 + depth_offset[depth + 2]].slectedCu = EB_TRUE;
            local_cu_array[cu_index + 1 + 3 * depth_offset[depth + 1] + 1 + 2 * depth_offset[depth + 2]].slectedCu = EB_TRUE;
            local_cu_array[cu_index + 1 + 3 * depth_offset[depth + 1] + 1 + 3 * depth_offset[depth + 2]].slectedCu = EB_TRUE;
        }
        if (lowestLevel == REFINEMENT_Pp2) {
            if (depth < 2 && cu_index < 65) {
                local_cu_array[cu_index + 1 + 1].stopSplit = EB_TRUE;
                local_cu_array[cu_index + 1 + 1 + depth_offset[depth + 2]].stopSplit = EB_TRUE;
                local_cu_array[cu_index + 1 + 1 + 2 * depth_offset[depth + 2]].stopSplit = EB_TRUE;
                local_cu_array[cu_index + 1 + 1 + 3 * depth_offset[depth + 2]].stopSplit = EB_TRUE;

                local_cu_array[cu_index + 1 + depth_offset[depth + 1] + 1].stopSplit = EB_TRUE;
                local_cu_array[cu_index + 1 + depth_offset[depth + 1] + 1 + depth_offset[depth + 2]].stopSplit = EB_TRUE;
                local_cu_array[cu_index + 1 + depth_offset[depth + 1] + 1 + 2 * depth_offset[depth + 2]].stopSplit = EB_TRUE;
                local_cu_array[cu_index + 1 + depth_offset[depth + 1] + 1 + 3 * depth_offset[depth + 2]].stopSplit = EB_TRUE;

                local_cu_array[cu_index + 1 + 2 * depth_offset[depth + 1] + 1].stopSplit = EB_TRUE;
                local_cu_array[cu_index + 1 + 2 * depth_offset[depth + 1] + 1 + depth_offset[depth + 2]].stopSplit = EB_TRUE;
                local_cu_array[cu_index + 1 + 2 * depth_offset[depth + 1] + 1 + 2 * depth_offset[depth + 2]].stopSplit = EB_TRUE;
                local_cu_array[cu_index + 1 + 2 * depth_offset[depth + 1] + 1 + 3 * depth_offset[depth + 2]].stopSplit = EB_TRUE;

                local_cu_array[cu_index + 1 + 3 * depth_offset[depth + 1] + 1].stopSplit = EB_TRUE;
                local_cu_array[cu_index + 1 + 3 * depth_offset[depth + 1] + 1 + depth_offset[depth + 2]].stopSplit = EB_TRUE;
                local_cu_array[cu_index + 1 + 3 * depth_offset[depth + 1] + 1 + 2 * depth_offset[depth + 2]].stopSplit = EB_TRUE;
                local_cu_array[cu_index + 1 + 3 * depth_offset[depth + 1] + 1 + 3 * depth_offset[depth + 2]].stopSplit = EB_TRUE;
            }
        }
    }

    if (refinementLevel & REFINEMENT_Pp3) {
        uint8_t inLoop;
        uint8_t outLoop;
        uint8_t cu_index = 2;
        if (depth == 0) {

            for (outLoop = 0; outLoop < 16; ++outLoop) {
                for (inLoop = 0; inLoop < 4; ++inLoop) {
                    local_cu_array[++cu_index].slectedCu = EB_TRUE;

                }
                cu_index += cu_index == 21 ? 2 : cu_index == 42 ? 2 : cu_index == 63 ? 2 : 1;

            }
            if (lowestLevel == REFINEMENT_Pp3) {
                cu_index = 2;
                for (outLoop = 0; outLoop < 16; ++outLoop) {
                    for (inLoop = 0; inLoop < 4; ++inLoop) {
                        local_cu_array[++cu_index].stopSplit = EB_TRUE;
                    }
                    cu_index += cu_index == 21 ? 2 : cu_index == 42 ? 2 : cu_index == 63 ? 2 : 1;
                }
            }
        }

    }

    if (refinementLevel & REFINEMENT_Pm1) {
        if (depth > 0) {
            local_cu_array[cu_index - 1 - parentCuIndex[cu_index]].slectedCu = EB_TRUE;
        }
        if (lowestLevel == REFINEMENT_Pm1) {
            if (depth > 0) {
                local_cu_array[cu_index - 1 - parentCuIndex[cu_index]].stopSplit = EB_TRUE;
            }
        }
    }

    if (refinementLevel & REFINEMENT_Pm2) {
        if (depth == 2) {
            local_cu_array[0].slectedCu = EB_TRUE;
        }
        if (depth == 3) {
            local_cu_array[1].slectedCu = EB_TRUE;
            local_cu_array[22].slectedCu = EB_TRUE;
            local_cu_array[43].slectedCu = EB_TRUE;
            local_cu_array[64].slectedCu = EB_TRUE;
        }
        if (lowestLevel == REFINEMENT_Pm2) {
            if (depth == 2) {
                local_cu_array[0].stopSplit = EB_TRUE;
            }
            if (depth == 3) {
                local_cu_array[1].stopSplit = EB_TRUE;
                local_cu_array[22].stopSplit = EB_TRUE;
                local_cu_array[43].stopSplit = EB_TRUE;
                local_cu_array[64].stopSplit = EB_TRUE;
            }
        }
    }

    if (refinementLevel & REFINEMENT_Pm3) {
        if (depth == 3) {
            local_cu_array[0].slectedCu = EB_TRUE;
        }
        if (lowestLevel == REFINEMENT_Pm2) {
            if (depth == 3) {
                local_cu_array[0].stopSplit = EB_TRUE;
            }
        }
    }

    return return_error;
}
/*******************************************
Cost Computation for intra CU
*******************************************/
// TO be updated with AV1 functions
EbErrorType MdcIntraCuRate(
    uint32_t                                  cu_depth,
    MdRateEstimationContext               *md_rate_estimation_ptr,
    PictureControlSet                    *picture_control_set_ptr,
    uint64_t                                  *mdcIntraRate)
{
    EbErrorType return_error = EB_ErrorNone;

    (void)mdcIntraRate;
    (void)picture_control_set_ptr;
    (void)md_rate_estimation_ptr;
    (void)cu_depth;
    //uint32_t chroma_mode = EB_INTRA_CHROMA_DM;
    //EbPartMode partitionMode = SIZE_2Nx2N;
    //int32_t prediction_index = -1;

    //// Number of bits for each synatax element
    //uint64_t partSizeIntraBitsNum = 0;
    //uint64_t intraLumaModeBitsNum = 0;

    //// Luma and chroma rate
    //uint64_t lumaRate;
    //uint64_t chromaRate;

    //EncodeContext *encode_context_ptr = ((SequenceControlSet*)(picture_control_set_ptr->sequence_control_set_wrapper_ptr->object_ptr))->encode_context_ptr;

    //CHECK_REPORT_ERROR(
    //    (partitionMode == SIZE_2Nx2N),
    //    encode_context_ptr->app_callback_ptr,
    //    EB_ENC_RD_COST_ERROR3);

    //// Estimate Chroma Mode Bits
    //chromaRate = md_rate_estimation_ptr->intraChromaBits[chroma_mode];

    //// Estimate Partition Size Bits :
    //// *Note - Intra is implicitly 2Nx2N
    //partSizeIntraBitsNum = ((uint8_t)cu_depth == (((SequenceControlSet *)picture_control_set_ptr->sequence_control_set_wrapper_ptr->object_ptr)->max_sb_depth - 1)) ?
    //    md_rate_estimation_ptr->intraPartSizeBits[partitionMode] :
    //    ZERO_COST;

    //// Estimate Luma Mode Bits for Intra
    ///*intra_luma_mode_context(
    //cu_ptr,
    //luma_mode,
    //&prediction_index);*/

    //intraLumaModeBitsNum = (prediction_index != -1) ?
    //    md_rate_estimation_ptr->intraLumaBits[prediction_index] :
    //    md_rate_estimation_ptr->intraLumaBits[3];

    //// Rate of the candiadate mode is equal to the sum of the rate of independent syntax element
    //lumaRate =
    //    partSizeIntraBitsNum +
    //    intraLumaModeBitsNum;

    //// Assign fast cost
    //*(mdcIntraRate) = chromaRate + lumaRate;

    return return_error;
}

/*******************************************
Cost Computation for inter CU
*******************************************/
// TO be updated with AV1 functions
uint64_t MdcInterCuRate(
    uint32_t                                   direction,
    int16_t                                   x_mv_l0,
    int16_t                                   y_mv_l0,
    int16_t                                   x_mv_l1,
    int16_t                                   y_mv_l1)
{

    int32_t MVs_0, MVs_1, MVs_2, MVs_3;

    uint64_t rate = 86440;


    switch (direction)
    {
    default:
    case UNI_PRED_LIST_0:

        // Estimate the Motion Vector Prediction Index Bits
        rate += 23196;

        // Estimate the Motion Vector Difference Bits
        MVs_0 = ABS(x_mv_l0);
        MVs_1 = ABS(y_mv_l0);
        MVs_0 = MVs_0 > 499 ? 499 : MVs_0;
        MVs_1 = MVs_1 > 499 ? 499 : MVs_1;
        rate += mv_bit_table[MVs_0][MVs_1];
        break;

    case UNI_PRED_LIST_1:

        // Estimate the Motion Vector Prediction Index Bits
        rate += 23196;

        // Estimate the Motion Vector Difference Bits

        MVs_2 = ABS(x_mv_l1);
        MVs_3 = ABS(y_mv_l1);
        MVs_2 = MVs_2 > 499 ? 499 : MVs_2;
        MVs_3 = MVs_3 > 499 ? 499 : MVs_3;


        rate += mv_bit_table[MVs_2][MVs_3];




        break;

    case BI_PRED:
        //LIST 0 RATE ESTIMATION

        // Estimate the Motion Vector Prediction Index Bits

        rate += 46392;

        // Estimate the Motion Vector Difference Bits

        MVs_0 = ABS(x_mv_l0);
        MVs_1 = ABS(y_mv_l0);
        MVs_0 = MVs_0 > 499 ? 499 : MVs_0;
        MVs_1 = MVs_1 > 499 ? 499 : MVs_1;


        rate += mv_bit_table[MVs_0][MVs_1];

        // Estimate the Motion Vector Difference Bits
        MVs_2 = ABS(x_mv_l1);
        MVs_3 = ABS(y_mv_l1);
        MVs_2 = MVs_2 > 499 ? 499 : MVs_2;
        MVs_3 = MVs_3 > 499 ? 499 : MVs_3;


        rate += mv_bit_table[MVs_2][MVs_3];

        break;

    }

    return rate;
}


/*******************************************
Derive the contouring class
If (AC energy < 32 * 32) then apply aggressive action (Class 1),
else if (AC energy < 32 * 32 * 1.6) OR (32 * 32 * 3.5 < AC energy < 32 * 32 * 4.5 AND non-8x8) then moderate action (Class 2),
else no action
*******************************************/
uint8_t derive_contouring_class(
    PictureParentControlSet   *parent_pcs_ptr,
    uint16_t                       sb_index,
    uint8_t                        leaf_index)
{
    uint8_t contouringClass = 0;

    SequenceControlSet *sequence_control_set_ptr = (SequenceControlSet*)parent_pcs_ptr->sequence_control_set_wrapper_ptr->object_ptr;

    if (parent_pcs_ptr->is_sb_homogeneous_over_time[sb_index]) {
        if (leaf_index > 0) {
            SbParams            *sb_params = &sequence_control_set_ptr->sb_params_array[sb_index];
            if (sb_params->is_edge_sb) {

                if (parent_pcs_ptr->sb_y_src_energy_cu_array[sb_index][(leaf_index - 1) / 21 + 1] < ANTI_CONTOURING_TH_1) {
                    contouringClass = 2;
                }
                else if (parent_pcs_ptr->sb_y_src_energy_cu_array[sb_index][(leaf_index - 1) / 21 + 1] < ANTI_CONTOURING_TH_2) {
                    contouringClass = 3;
                }
                else if (parent_pcs_ptr->sb_y_src_energy_cu_array[sb_index][(leaf_index - 1) / 21 + 1] < (ANTI_CONTOURING_TH_1 + ANTI_CONTOURING_TH_2)) {
                    contouringClass = 3;
                }
            }
            else {
                if (parent_pcs_ptr->sb_y_src_energy_cu_array[sb_index][(leaf_index - 1) / 21 + 1] < ANTI_CONTOURING_TH_0) {
                    contouringClass = 1;
                }
                else if (parent_pcs_ptr->sb_y_src_energy_cu_array[sb_index][(leaf_index - 1) / 21 + 1] < ANTI_CONTOURING_TH_1) {
                    contouringClass = 2;
                }
                else if (parent_pcs_ptr->sb_y_src_energy_cu_array[sb_index][(leaf_index - 1) / 21 + 1] < ANTI_CONTOURING_TH_2) {
                    contouringClass = 3;
                }
            }
        }
    }
    return(contouringClass);
}


void RefinementPredictionLoop(
    SequenceControlSet                   *sequence_control_set_ptr,
    PictureControlSet                    *picture_control_set_ptr,
    LargestCodingUnit                    *sb_ptr,
    uint32_t                                  sb_index,
    ModeDecisionConfigurationContext     *context_ptr)
{

    MdcpLocalCodingUnit  *local_cu_array = context_ptr->local_cu_array;
    SbParams            *sb_params = &sequence_control_set_ptr->sb_params_array[sb_index];
    uint32_t                  temporal_layer_index = picture_control_set_ptr->temporal_layer_index;
    uint32_t                  cu_index = 0;
    uint8_t                   stationary_edge_over_time_flag = (&(picture_control_set_ptr->parent_pcs_ptr->sb_stat_array[sb_index]))->stationary_edge_over_time_flag;

    sb_ptr->pred64 = EB_FALSE;
    while (cu_index < CU_MAX_COUNT)
    {
        if (sb_params->raster_scan_cu_validity[md_scan_to_raster_scan[cu_index]] && (local_cu_array[cu_index].earlySplitFlag == EB_FALSE))
        {
            local_cu_array[cu_index].slectedCu = EB_TRUE;
            sb_ptr->pred64 = (cu_index == 0) ? EB_TRUE : sb_ptr->pred64;
            uint32_t depth = get_coded_unit_stats(cu_index)->depth;
            uint8_t refinementLevel;

            if (sb_ptr->picture_control_set_ptr->slice_type == I_SLICE) {

                {
                    uint8_t lowestLevel = 0x00;

                    if (sequence_control_set_ptr->input_resolution == INPUT_SIZE_4K_RANGE)
                        refinementLevel = ndp_refinement_control_islice_m4[0][depth]; // HG: why always 0
                    else
                        refinementLevel = ndp_refinement_control_islice_1080_p_m4[0][depth]; // HG: why always 0

                    if (depth <= 1 && stationary_edge_over_time_flag > 0) {
                        if (depth == 0)
                            refinementLevel = Predp1 + Predp2 + Predp3;
                        else
                            refinementLevel = Pred + Predp1 + Predp2;

                    }
                    lowestLevel = (refinementLevel & REFINEMENT_Pp3) ? REFINEMENT_Pp3 : (refinementLevel & REFINEMENT_Pp2) ? REFINEMENT_Pp2 : (refinementLevel & REFINEMENT_Pp1) ? REFINEMENT_Pp1 :
                        (refinementLevel & REFINEMENT_P) ? REFINEMENT_P :
                        (refinementLevel & REFINEMENT_Pm1) ? REFINEMENT_Pm1 : (refinementLevel & REFINEMENT_Pm2) ? REFINEMENT_Pm2 : (refinementLevel & REFINEMENT_Pm3) ? REFINEMENT_Pm3 : 0x00;

                    MdcRefinement(
                        &(*context_ptr->local_cu_array),
                        cu_index,
                        depth,
                        refinementLevel,
                        lowestLevel);
                }
            }
            else {
                if (picture_control_set_ptr->parent_pcs_ptr->pic_depth_mode == PIC_SB_SWITCH_DEPTH_MODE && (picture_control_set_ptr->parent_pcs_ptr->sb_md_mode_array[sb_index] == LCU_PRED_OPEN_LOOP_DEPTH_MODE || picture_control_set_ptr->parent_pcs_ptr->sb_md_mode_array[sb_index] == LCU_PRED_OPEN_LOOP_1_NFL_DEPTH_MODE)) {
                    refinementLevel = Pred;
                }
                else


                    if (picture_control_set_ptr->parent_pcs_ptr->pic_depth_mode == PIC_OPEN_LOOP_DEPTH_MODE ||
                        (picture_control_set_ptr->parent_pcs_ptr->pic_depth_mode == PIC_SB_SWITCH_DEPTH_MODE && (picture_control_set_ptr->parent_pcs_ptr->sb_md_mode_array[sb_index] == LCU_OPEN_LOOP_DEPTH_MODE || picture_control_set_ptr->parent_pcs_ptr->sb_md_mode_array[sb_index] == LCU_AVC_DEPTH_MODE)))

                        refinementLevel = ndp_refinement_control_nref[temporal_layer_index][depth];
                    else
                        refinementLevel = ndp_refinement_control_fast[temporal_layer_index][depth];

                if (picture_control_set_ptr->parent_pcs_ptr->cu8x8_mode == CU_8x8_MODE_1) {
                    refinementLevel = ((refinementLevel & REFINEMENT_Pp1) && depth == 2) ? refinementLevel - REFINEMENT_Pp1 :
                        ((refinementLevel & REFINEMENT_Pp2) && depth == 1) ? refinementLevel - REFINEMENT_Pp2 :
                        ((refinementLevel & REFINEMENT_Pp3) && depth == 0) ? refinementLevel - REFINEMENT_Pp3 : refinementLevel;
                }

                uint8_t lowestLevel = 0x00;

                lowestLevel = (refinementLevel & REFINEMENT_Pp3) ? REFINEMENT_Pp3 : (refinementLevel & REFINEMENT_Pp2) ? REFINEMENT_Pp2 : (refinementLevel & REFINEMENT_Pp1) ? REFINEMENT_Pp1 :
                    (refinementLevel & REFINEMENT_P) ? REFINEMENT_P :
                    (refinementLevel & REFINEMENT_Pm1) ? REFINEMENT_Pm1 : (refinementLevel & REFINEMENT_Pm2) ? REFINEMENT_Pm2 : (refinementLevel & REFINEMENT_Pm3) ? REFINEMENT_Pm3 : 0x00;

                MdcRefinement(
                    &(*context_ptr->local_cu_array),
                    cu_index,
                    depth,
                    refinementLevel,
                    lowestLevel);
            }

            cu_index += depth_offset[depth];

        }
        else {

            cu_index++;
        }
    } // End while 1 CU Loop
}


void PrePredictionRefinement(
    SequenceControlSet                   *sequence_control_set_ptr,
    PictureControlSet                    *picture_control_set_ptr,
    LargestCodingUnit                    *sb_ptr,
    uint32_t                                  sb_index,
    uint32_t                                 *startDepth,
    uint32_t                                 *endDepth
)
{
    SbParams    *sb_params = &sequence_control_set_ptr->sb_params_array[sb_index];

    EB_SLICE        slice_type = picture_control_set_ptr->slice_type;

    uint8_t           edge_block_num = picture_control_set_ptr->parent_pcs_ptr->edge_results_ptr[sb_index].edge_block_num;

    SbStat      *sb_stat_ptr = &(picture_control_set_ptr->parent_pcs_ptr->sb_stat_array[sb_index]);
    uint8_t           stationary_edge_over_time_flag = sb_stat_ptr->stationary_edge_over_time_flag;

    uint8_t           aura_status_iii = sb_ptr->aura_status_iii;






    if (picture_control_set_ptr->parent_pcs_ptr->high_dark_low_light_area_density_flag && picture_control_set_ptr->parent_pcs_ptr->temporal_layer_index > 0 && picture_control_set_ptr->parent_pcs_ptr->sharp_edge_sb_flag[sb_index] && !picture_control_set_ptr->parent_pcs_ptr->similar_colocated_sb_array_ii[sb_index]) {
        *startDepth = DEPTH_16;
    }

    if ((slice_type != I_SLICE && picture_control_set_ptr->high_intra_slection == 0) && (sb_params->is_complete_sb)) {

        if (picture_control_set_ptr->scene_caracteristic_id == EB_FRAME_CARAC_0) {

            if (picture_control_set_ptr->parent_pcs_ptr->grass_percentage_in_picture > 60 && aura_status_iii) {
                *startDepth = DEPTH_16;
            }
        }
    }

    if (picture_control_set_ptr->parent_pcs_ptr->logo_pic_flag && edge_block_num)
    {
        *startDepth = DEPTH_16;
    }


    // S-LOGO

    if (stationary_edge_over_time_flag > 0) {

        *startDepth = DEPTH_16;
        *endDepth = DEPTH_16;

    }

    if (picture_control_set_ptr->parent_pcs_ptr->complex_sb_array[sb_ptr->index] == SB_COMPLEXITY_STATUS_2) {
        *startDepth = DEPTH_16;
    }
}


void ForwardCuToModeDecision(
    SequenceControlSet                   *sequence_control_set_ptr,
    PictureControlSet                    *picture_control_set_ptr,

    uint32_t                                  sb_index,
    ModeDecisionConfigurationContext     *context_ptr
)
{

    uint8_t                   cu_index = 0;
    uint32_t                  cuClass = DO_NOT_ADD_CU_CONTINUE_SPLIT;
    EbBool                 split_flag = EB_TRUE;
    MdcLcuData           *resultsPtr = &picture_control_set_ptr->mdc_sb_array[sb_index];
    SbParams            *sb_params = &sequence_control_set_ptr->sb_params_array[sb_index];
    MdcpLocalCodingUnit  *local_cu_array = context_ptr->local_cu_array;
    EB_SLICE                slice_type = picture_control_set_ptr->slice_type;


    // CU Loop
    const CodedUnitStats *cuStatsPtr = get_coded_unit_stats(0);

    SbStat *sb_stat_ptr = &(picture_control_set_ptr->parent_pcs_ptr->sb_stat_array[sb_index]);

    EbBool    testAllDepthIntraSliceFlag = EB_FALSE;


    testAllDepthIntraSliceFlag = slice_type == I_SLICE &&
        (sb_stat_ptr->stationary_edge_over_time_flag || picture_control_set_ptr->parent_pcs_ptr->logo_pic_flag ||
        (picture_control_set_ptr->parent_pcs_ptr->very_low_var_pic_flag && picture_control_set_ptr->parent_pcs_ptr->low_motion_content_flag)) ?
        EB_TRUE : testAllDepthIntraSliceFlag;


    resultsPtr->leaf_count = 0;

    cu_index = 0;

    while (cu_index < CU_MAX_COUNT)
    {
        split_flag = EB_TRUE;
        if (sb_params->raster_scan_cu_validity[md_scan_to_raster_scan[cu_index]])
        {
            cuStatsPtr = get_coded_unit_stats(cu_index);

            switch (cuStatsPtr->depth) {

            case 0:
            case 1:
            case 2:

                cuClass = DO_NOT_ADD_CU_CONTINUE_SPLIT;


                if (slice_type == I_SLICE) {
                    if (testAllDepthIntraSliceFlag) {
                        cuClass = ADD_CU_CONTINUE_SPLIT;
                    }
                    else {
                        cuClass = local_cu_array[cu_index].slectedCu == EB_TRUE ? ADD_CU_CONTINUE_SPLIT : cuClass;
                        cuClass = local_cu_array[cu_index].stopSplit == EB_TRUE ? ADD_CU_STOP_SPLIT : cuClass;
                    }
                }
                else {
                    cuClass = local_cu_array[cu_index].slectedCu == EB_TRUE ? ADD_CU_CONTINUE_SPLIT : cuClass;
                    cuClass = local_cu_array[cu_index].stopSplit == EB_TRUE ? ADD_CU_STOP_SPLIT : cuClass;

                }

                // Take into account MAX CU size & MAX intra size (from the API)
                cuClass = (cuStatsPtr->size > sequence_control_set_ptr->max_cu_size || (slice_type == I_SLICE && cuStatsPtr->size > sequence_control_set_ptr->max_intra_size)) ?
                    DO_NOT_ADD_CU_CONTINUE_SPLIT :
                    cuClass;

                // Take into account MIN CU size & Min intra size(from the API)
                cuClass = (cuStatsPtr->size == sequence_control_set_ptr->min_cu_size || (slice_type == I_SLICE && cuStatsPtr->size == sequence_control_set_ptr->min_intra_size)) ?
                    ADD_CU_STOP_SPLIT :
                    cuClass;

                switch (cuClass) {

                case ADD_CU_STOP_SPLIT:
                    // Stop
                    resultsPtr->leaf_data_array[resultsPtr->leaf_count].leaf_index = cu_index;
                    resultsPtr->leaf_data_array[resultsPtr->leaf_count++].split_flag = split_flag = EB_FALSE;

                    break;

                case ADD_CU_CONTINUE_SPLIT:
                    // Go Down + consider the current CU as candidate
                    resultsPtr->leaf_data_array[resultsPtr->leaf_count].leaf_index = cu_index;
                    resultsPtr->leaf_data_array[resultsPtr->leaf_count++].split_flag = split_flag = EB_TRUE;

                    break;

                case DO_NOT_ADD_CU_CONTINUE_SPLIT:
                    // Go Down + do not consider the current CU as candidate
                    split_flag = EB_TRUE;

                    break;

                default:
                    resultsPtr->leaf_data_array[resultsPtr->leaf_count].leaf_index = cu_index;
                    resultsPtr->leaf_data_array[resultsPtr->leaf_count++].split_flag = split_flag = EB_TRUE;

                    break;
                }

                break;
            case 3:

                resultsPtr->leaf_data_array[resultsPtr->leaf_count].leaf_index = cu_index;
                resultsPtr->leaf_data_array[resultsPtr->leaf_count++].split_flag = split_flag = EB_FALSE;

                break;

            default:
                resultsPtr->leaf_data_array[resultsPtr->leaf_count].leaf_index = cu_index;
                resultsPtr->leaf_data_array[resultsPtr->leaf_count++].split_flag = split_flag = EB_TRUE;
                break;
            }
        }

        cu_index += (split_flag == EB_TRUE) ? 1 : depth_offset[cuStatsPtr->depth];

    } // End CU Loop

}



void MdcInterDepthDecision(
    ModeDecisionConfigurationContext     *context_ptr,
    uint32_t                                 origin_x,
    uint32_t                                 origin_y,
    uint32_t                                 endDepth,
    uint32_t                                 splitFlagBits0,
    uint32_t                                 splitFlagBits1,
    uint64_t                                 lambda,

    uint32_t                                 cu_index
)
{

    uint32_t               leftCuIndex;
    uint32_t               topCuIndex;
    uint32_t               topLeftCuIndex;
    uint32_t               depthZeroCandidateCuIndex;
    uint32_t               depthOneCandidateCuIndex = cu_index;
    uint32_t               depthTwoCandidateCuIndex = cu_index;
    uint64_t               depthNRate = 0;
    uint64_t               depthNPlusOneRate = 0;
    uint64_t               depthNCost = 0;
    uint64_t               depthNPlusOneCost = 0;
    MdcpLocalCodingUnit *local_cu_array = context_ptr->local_cu_array;
    /*** Stage 0: Inter depth decision: depth 2 vs depth 3 ***/
    // Walks to the last coded 8x8 block for merging
    uint8_t  group_of8x8_blocks_count = context_ptr->group_of8x8_blocks_count;
    uint8_t  group_of16x16_blocks_count = context_ptr->group_of16x16_blocks_count;
    if ((GROUP_OF_4_8x8_BLOCKS(origin_x, origin_y))) {

        group_of8x8_blocks_count++;

        // From the last coded cu index, get the indices of the left, top, and top left cus
        leftCuIndex = cu_index - DEPTH_THREE_STEP;
        topCuIndex = leftCuIndex - DEPTH_THREE_STEP;
        topLeftCuIndex = topCuIndex - DEPTH_THREE_STEP;

        // From the top left index, get the index of the candidate pu for merging
        depthTwoCandidateCuIndex = topLeftCuIndex - 1;

        // Compute depth N cost
        local_cu_array[depthTwoCandidateCuIndex].splitContext = 0;

        depthNRate = (((lambda  * splitFlagBits0) + MD_OFFSET) >> MD_SHIFT);

        depthNCost = (local_cu_array[depthTwoCandidateCuIndex]).earlyCost + depthNRate;

        if (endDepth < 3) {

            (local_cu_array[depthTwoCandidateCuIndex]).earlySplitFlag = EB_FALSE;
            (local_cu_array[depthTwoCandidateCuIndex]).earlyCost = depthNCost;

        }
        else {

            // Assign rate
            depthNPlusOneRate = (((lambda  * splitFlagBits1) + MD_OFFSET) >> MD_SHIFT);

            depthNPlusOneCost = (local_cu_array[cu_index]).earlyCost + (local_cu_array[leftCuIndex]).earlyCost + (local_cu_array[topCuIndex]).earlyCost + (local_cu_array[topLeftCuIndex]).earlyCost + depthNPlusOneRate;

            if (depthNCost <= depthNPlusOneCost) {

                // If the cost is low enough to warrant not spliting further:
                // 1. set the split flag of the candidate pu for merging to false
                // 2. update the last pu index
                (local_cu_array[depthTwoCandidateCuIndex]).earlySplitFlag = EB_FALSE;
                (local_cu_array[depthTwoCandidateCuIndex]).earlyCost = depthNCost;

            }
            else {
                // If the cost is not low enough:
                // update the cost of the candidate pu for merging
                // this update is required for the next inter depth decision
                (&local_cu_array[depthTwoCandidateCuIndex])->earlyCost = depthNPlusOneCost;
            }

        }
    }

    // Walks to the last coded 16x16 block for merging
    if (GROUP_OF_4_16x16_BLOCKS(get_coded_unit_stats(depthTwoCandidateCuIndex)->origin_x, get_coded_unit_stats(depthTwoCandidateCuIndex)->origin_y) &&
        (group_of8x8_blocks_count == 4)) {

        group_of8x8_blocks_count = 0;
        group_of16x16_blocks_count++;

        // From the last coded pu index, get the indices of the left, top, and top left pus
        leftCuIndex = depthTwoCandidateCuIndex - DEPTH_TWO_STEP;
        topCuIndex = leftCuIndex - DEPTH_TWO_STEP;
        topLeftCuIndex = topCuIndex - DEPTH_TWO_STEP;

        // From the top left index, get the index of the candidate pu for merging
        depthOneCandidateCuIndex = topLeftCuIndex - 1;

        if (get_coded_unit_stats(depthOneCandidateCuIndex)->depth == 1) {
            depthNRate = (((lambda  *splitFlagBits0) + MD_OFFSET) >> MD_SHIFT);

            depthNCost = local_cu_array[depthOneCandidateCuIndex].earlyCost + depthNRate;
            if (endDepth < 2) {

                local_cu_array[depthOneCandidateCuIndex].earlySplitFlag = EB_FALSE;
                local_cu_array[depthOneCandidateCuIndex].earlyCost = depthNCost;

            }
            else {
                // Compute depth N+1 cost
                // Assign rate
                depthNPlusOneRate = (((lambda  *splitFlagBits1) + MD_OFFSET) >> MD_SHIFT);

                depthNPlusOneCost = local_cu_array[depthTwoCandidateCuIndex].earlyCost +
                    local_cu_array[leftCuIndex].earlyCost +
                    local_cu_array[topCuIndex].earlyCost +
                    local_cu_array[topLeftCuIndex].earlyCost +
                    depthNPlusOneRate;

                // Inter depth comparison: depth 1 vs depth 2
                if (depthNCost <= depthNPlusOneCost) {
                    // If the cost is low enough to warrant not spliting further:
                    // 1. set the split flag of the candidate pu for merging to false
                    // 2. update the last pu index
                    local_cu_array[depthOneCandidateCuIndex].earlySplitFlag = EB_FALSE;
                    local_cu_array[depthOneCandidateCuIndex].earlyCost = depthNCost;
                }
                else {
                    // If the cost is not low enough:
                    // update the cost of the candidate pu for merging
                    // this update is required for the next inter depth decision
                    local_cu_array[depthOneCandidateCuIndex].earlyCost = depthNPlusOneCost;
                }
            }
        }
    }

    // Stage 2: Inter depth decision: depth 0 vs depth 1

    // Walks to the last coded 32x32 block for merging
    // Stage 2 isn't performed in I slices since the abcense of 64x64 candidates
    if (GROUP_OF_4_32x32_BLOCKS(get_coded_unit_stats(depthOneCandidateCuIndex)->origin_x, get_coded_unit_stats(depthOneCandidateCuIndex)->origin_y) &&
        (group_of16x16_blocks_count == 4)) {

        group_of16x16_blocks_count = 0;

        // From the last coded pu index, get the indices of the left, top, and top left pus
        leftCuIndex = depthOneCandidateCuIndex - DEPTH_ONE_STEP;
        topCuIndex = leftCuIndex - DEPTH_ONE_STEP;
        topLeftCuIndex = topCuIndex - DEPTH_ONE_STEP;

        // From the top left index, get the index of the candidate pu for merging
        depthZeroCandidateCuIndex = topLeftCuIndex - 1;

        if (get_coded_unit_stats(depthZeroCandidateCuIndex)->depth == 0) {

            // Compute depth N cost
            depthNRate = (((lambda  *splitFlagBits0) + MD_OFFSET) >> MD_SHIFT);

            depthNCost = (&local_cu_array[depthZeroCandidateCuIndex])->earlyCost + depthNRate;
            if (endDepth < 1) {

                (&local_cu_array[depthZeroCandidateCuIndex])->earlySplitFlag = EB_FALSE;

            }
            else {
                // Compute depth N+1 cost

                // Assign rate
                depthNPlusOneRate = (((lambda  *splitFlagBits1) + MD_OFFSET) >> MD_SHIFT);

                depthNPlusOneCost = local_cu_array[depthOneCandidateCuIndex].earlyCost +
                    local_cu_array[leftCuIndex].earlyCost +
                    local_cu_array[topCuIndex].earlyCost +
                    local_cu_array[topLeftCuIndex].earlyCost +
                    depthNPlusOneRate;

                // Inter depth comparison: depth 0 vs depth 1
                if (depthNCost <= depthNPlusOneCost) {
                    // If the cost is low enough to warrant not spliting further:
                    // 1. set the split flag of the candidate pu for merging to false
                    // 2. update the last pu index
                    (&local_cu_array[depthZeroCandidateCuIndex])->earlySplitFlag = EB_FALSE;
                }
            }
        }
    }

    context_ptr->group_of8x8_blocks_count = group_of8x8_blocks_count;
    context_ptr->group_of16x16_blocks_count = group_of16x16_blocks_count;
}

void PredictionPartitionLoop(
    SequenceControlSet                   *sequence_control_set_ptr,
    PictureControlSet                    *picture_control_set_ptr,
    uint32_t                                  sb_index,
    uint32_t                                  tbOriginX,
    uint32_t                                  tbOriginY,
    uint32_t                                  startDepth,
    uint32_t                                  endDepth,
    ModeDecisionConfigurationContext     *context_ptr
)
{

    MdRateEstimationContext *md_rate_estimation_ptr = context_ptr->md_rate_estimation_ptr;
    MdcpLocalCodingUnit *local_cu_array = context_ptr->local_cu_array;
    MdcpLocalCodingUnit   *cu_ptr;

    SbParams *sb_params = &sequence_control_set_ptr->sb_params_array[sb_index];
    uint32_t      cuInterSad = 0;
    uint64_t      cuInterRate = 0;
    uint32_t      cuIntraSad = 0;
    uint64_t      cuIntraRate = 0;
    uint64_t      cuIntraCost = 0;
    uint32_t     cuIndexInRaterScan;
    uint64_t      cuInterCost = 0;
    uint32_t cu_index = 0;
    uint32_t start_index = 0;

    (void)tbOriginX;
    (void)tbOriginY;

    const CodedUnitStats *cuStatsPtr;

    for (cu_index = start_index; cu_index < CU_MAX_COUNT; ++cu_index)

    {

        local_cu_array[cu_index].slectedCu = EB_FALSE;
        local_cu_array[cu_index].stopSplit = EB_FALSE;

        cu_ptr = &local_cu_array[cu_index];
        cuIndexInRaterScan = md_scan_to_raster_scan[cu_index];
        if (sb_params->raster_scan_cu_validity[cuIndexInRaterScan])
        {
            uint32_t depth;
            uint32_t size;
            cuStatsPtr = get_coded_unit_stats(cu_index);

            depth = cuStatsPtr->depth;
            size = cuStatsPtr->size;
            cu_ptr->earlySplitFlag = (depth < endDepth) ? EB_TRUE : EB_FALSE;

            if (depth >= startDepth && depth <= endDepth) {
                //reset the flags here:   all CU splitFalg=TRUE. default: we always split. interDepthDecision will select where  to stop splitting(ie setting the flag to False)

                {
                    MdcIntraCuRate(
                        depth,
                        md_rate_estimation_ptr,
                        picture_control_set_ptr,
                        &cuIntraRate);


                    OisCu32Cu16Results            *oisCu32Cu16ResultsPtr = picture_control_set_ptr->parent_pcs_ptr->ois_cu32_cu16_results[sb_index];
                    OisCu8Results                   *oisCu8ResultsPtr = picture_control_set_ptr->parent_pcs_ptr->ois_cu8_results[sb_index];
                    OisCandidate * OisCuPtr = cuIndexInRaterScan < RASTER_SCAN_CU_INDEX_8x8_0 ?
                        oisCu32Cu16ResultsPtr->sorted_ois_candidate[cuIndexInRaterScan] : oisCu8ResultsPtr->sorted_ois_candidate[cuIndexInRaterScan - RASTER_SCAN_CU_INDEX_8x8_0];


                    if (size > 32) {
                        cuIntraSad =
                            oisCu32Cu16ResultsPtr->sorted_ois_candidate[1][0].distortion +
                            oisCu32Cu16ResultsPtr->sorted_ois_candidate[2][0].distortion +
                            oisCu32Cu16ResultsPtr->sorted_ois_candidate[3][0].distortion +
                            oisCu32Cu16ResultsPtr->sorted_ois_candidate[4][0].distortion;
                    }
                    else if (size == 32) {
                        cuIntraSad = oisCu32Cu16ResultsPtr->sorted_ois_candidate[cuIndexInRaterScan][0].distortion;
                    }
                    else {
                        if (size > 8) {
                            cuIntraSad = oisCu32Cu16ResultsPtr->sorted_ois_candidate[cuIndexInRaterScan][0].distortion;
                        }
                        else {
                            if (OisCuPtr[0].valid_distortion) {
                                cuIntraSad = OisCuPtr[0].distortion;
                            }
                            else {

                                const CodedUnitStats  *cu_stats = get_coded_unit_stats(parent_block_index[cu_index]);
                                const uint32_t me2Nx2NTableOffset = cu_stats->cu_num_in_depth + me2Nx2NOffset[cu_stats->depth];


                                if (oisCu8ResultsPtr->sorted_ois_candidate[me2Nx2NTableOffset][0].valid_distortion) {
                                    cuIntraSad = oisCu8ResultsPtr->sorted_ois_candidate[me2Nx2NTableOffset][0].distortion;
                                }
                                else {
                                    cuIntraSad = 0;
                                }
                            }
                        }
                    }


                    cuIntraCost = (cuIntraSad << COST_PRECISION) + ((context_ptr->lambda * cuIntraRate + MD_OFFSET) >> MD_SHIFT);
                    cu_ptr->earlyCost = cuIntraCost;

                }

                if (picture_control_set_ptr->slice_type != I_SLICE) {



                    MeCuResults * mePuResult = &picture_control_set_ptr->parent_pcs_ptr->me_results[sb_index][cuIndexInRaterScan];
                    cuInterRate = MdcInterCuRate(
                        mePuResult->distortion_direction[0].direction,
                        mePuResult->x_mv_l0,
                        mePuResult->y_mv_l0,
                        mePuResult->x_mv_l1,
                        mePuResult->y_mv_l1);
                    cuInterSad = mePuResult->distortion_direction[0].distortion;


                    cuInterCost = (cuInterSad << COST_PRECISION) + ((context_ptr->lambda * cuInterRate + MD_OFFSET) >> MD_SHIFT);
                    cu_ptr->earlyCost = cuInterCost;

                }
#if ENCODER_MODE_CLEANUP
                if (1){
#else
                if (picture_control_set_ptr->enc_mode <= ENC_M3) {
#endif

                    cu_ptr->earlyCost = picture_control_set_ptr->slice_type == I_SLICE ? cuIntraCost : MIN(cuInterCost, cuIntraCost);
                }
                else {
                    cu_ptr->earlyCost = picture_control_set_ptr->slice_type == I_SLICE ? cuIntraCost : cuInterCost;
                }

                if (endDepth == 2) {
                    context_ptr->group_of8x8_blocks_count = depth == 2 ? incrementalCount[cuIndexInRaterScan] : 0;
                }

                if (endDepth == 1) {
                    context_ptr->group_of16x16_blocks_count = depth == 1 ? incrementalCount[cuIndexInRaterScan] : 0;
                }
                MdcInterDepthDecision(
                    context_ptr,
                    cuStatsPtr->origin_x,
                    cuStatsPtr->origin_y,
                    endDepth,
                    md_rate_estimation_ptr->splitFlagBits[0],
                    md_rate_estimation_ptr->splitFlagBits[3],
                    context_ptr->lambda,
                    cu_index);
            }
            else {
                cu_ptr->earlyCost = ~0u;
            }

        }

    }// End CU Loop

}

EbErrorType early_mode_decision_lcu(
    SequenceControlSet                   *sequence_control_set_ptr,
    PictureControlSet                    *picture_control_set_ptr,
    LargestCodingUnit                    *sb_ptr,
    uint32_t                                  sb_index,
    ModeDecisionConfigurationContext     *context_ptr)
{

    EbErrorType    return_error = EB_ErrorNone;

    uint32_t          tbOriginX = sb_ptr->origin_x;
    uint32_t          tbOriginY = sb_ptr->origin_y;
    EB_SLICE        slice_type = picture_control_set_ptr->slice_type;

    uint32_t      startDepth = (picture_control_set_ptr->parent_pcs_ptr->pic_depth_mode == PIC_SB_SWITCH_DEPTH_MODE && picture_control_set_ptr->parent_pcs_ptr->sb_md_mode_array[sb_index] == LCU_AVC_DEPTH_MODE) ?
        DEPTH_16 :
        DEPTH_64;


    uint32_t      endDepth = (slice_type == I_SLICE) ? DEPTH_8 : DEPTH_16;

    context_ptr->group_of8x8_blocks_count = 0;
    context_ptr->group_of16x16_blocks_count = 0;


    // The MDC refinements had been taken into account at the budgeting algorithm & therefore could be skipped for the ADP case
    if (picture_control_set_ptr->parent_pcs_ptr->pic_depth_mode != PIC_SB_SWITCH_DEPTH_MODE) {
        PrePredictionRefinement(
            sequence_control_set_ptr,
            picture_control_set_ptr,
            sb_ptr,
            sb_index,
            &startDepth,
            &endDepth);
    }

    PredictionPartitionLoop(
        sequence_control_set_ptr,
        picture_control_set_ptr,
        sb_index,
        tbOriginX,
        tbOriginY,
        startDepth,
        endDepth,
        context_ptr
    );

    RefinementPredictionLoop(
        sequence_control_set_ptr,
        picture_control_set_ptr,
        sb_ptr,
        sb_index,
        context_ptr);

    ForwardCuToModeDecision(
        sequence_control_set_ptr,
        picture_control_set_ptr,

        sb_index,
        context_ptr);

    return return_error;

}



