/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#include <stdlib.h>
#include "EbPictureManagerQueue.h"

EbErrorType input_queue_entry_ctor(InputQueueEntry *entry_ptr) {
    (void)entry_ptr;
    return EB_ErrorNone;
}

static void reference_queue_entry_dctor(EbPtr p) {
    ReferenceQueueEntry *obj = (ReferenceQueueEntry *)p;
    EB_FREE(obj->list0.list);
    EB_FREE(obj->list1.list);
}

EbErrorType reference_queue_entry_ctor(ReferenceQueueEntry *entry_ptr) {
    entry_ptr->dctor                = reference_queue_entry_dctor;
    entry_ptr->reference_object_ptr = (EbObjectWrapper *)EB_NULL;
    entry_ptr->picture_number       = ~0u;
    entry_ptr->dependent_count      = 0;
    entry_ptr->reference_available  = EB_FALSE;

    EB_MALLOC(entry_ptr->list0.list, sizeof(int32_t) * (1 << MAX_TEMPORAL_LAYERS));

    EB_MALLOC(entry_ptr->list1.list, sizeof(int32_t) * (1 << MAX_TEMPORAL_LAYERS));

    return EB_ErrorNone;
}
