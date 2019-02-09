// INTEL CONFIDENTIAL
// Copyright © 2018 Intel Corporation.
//
// This software and the related documents are Intel copyrighted materials,
// and your use of them is governed by the express license under which they were provided to you.
// Unless the License provides otherwise, you may not use, modify, copy, publish, distribute, disclose or transmit
// this software or the related documents without Intel's prior written permission.
// This software and the related documents are provided as is, with no express or implied warranties,
// other than those that are expressly stated in the License.

#ifndef EbCdefProcess_h
#define EbCdefProcess_h

#include "EbDefinitions.h"

#include "EbSystemResourceManager.h"
#include "EbPictureBufferDesc.h"
#include "EbSequenceControlSet.h"
#include "EbUtility.h"
#include "EbPsnr.h"
#include "EbPictureControlSet.h"

/**************************************
 * Cdef Context
 **************************************/
typedef struct CdefContext_s
{
    EbFifo_t                       *cdef_input_fifo_ptr;
    EbFifo_t                       *cdef_output_fifo_ptr;



} CdefContext_t;

/**************************************
 * Extern Function Declarations
 **************************************/
extern EbErrorType cdef_context_ctor(
    CdefContext_t **context_dbl_ptr,
    EbFifo_t                       *cdef_input_fifo_ptr,
    EbFifo_t                       *cdef_output_fifo_ptr,
    EbBool                  is16bit,
    uint32_t                max_input_luma_width,
    uint32_t                max_input_luma_height
   );

extern void* CdefKernel(void *input_ptr);

#endif