/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#include <stdlib.h>
#include "EbPictureDecisionQueue.h"

static void pa_reference_queue_entry_dctor(EbPtr p) {
    PaReferenceQueueEntry* obj = (PaReferenceQueueEntry*)p;
    EB_FREE(obj->list0.list);
    EB_FREE(obj->list1.list);
}

EbErrorType pa_reference_queue_entry_ctor(PaReferenceQueueEntry* entry_ptr) {
    entry_ptr->dctor = pa_reference_queue_entry_dctor;
    EB_MALLOC(entry_ptr->list0.list, sizeof(int32_t) * (1 << MAX_TEMPORAL_LAYERS));
    EB_MALLOC(entry_ptr->list1.list, sizeof(int32_t) * (1 << MAX_TEMPORAL_LAYERS));

    return EB_ErrorNone;
}
