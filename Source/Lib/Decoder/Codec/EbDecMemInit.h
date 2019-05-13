/*
* Copyright(c) 2019 Netflix, Inc.
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbDecMemInit_h
#define EbDecMemInit_h

#ifdef __cplusplus
extern "C" {
#endif

extern EbMemoryMapEntry                 *svt_dec_memory_map;
extern uint32_t                         *svt_dec_memory_map_index;
extern uint64_t                         *svt_dec_total_lib_memory;
extern uint32_t                         svt_dec_lib_malloc_count;

#ifdef _MSC_VER
#define EB_ALLIGN_MALLOC_DEC(type, pointer, n_elements, pointer_class) \
pointer = (type) _aligned_malloc(n_elements,ALVALUE); \
if (pointer == (type)EB_NULL) { \
    return EB_ErrorInsufficientResources; \
    } \
    else { \
    svt_dec_memory_map[*(svt_dec_memory_map_index)].ptr_type = pointer_class; \
    svt_dec_memory_map[(*(svt_dec_memory_map_index))++].ptr = pointer; \
    if (n_elements % 8 == 0) { \
        *svt_dec_total_lib_memory += (n_elements); \
    } \
    else { \
        *svt_dec_total_lib_memory += ((n_elements) + (8 - ((n_elements) % 8))); \
    } \
} \
if (*(svt_dec_memory_map_index) >= MAX_NUM_PTR) { \
    return EB_ErrorInsufficientResources; \
} \
svt_dec_lib_malloc_count++;

#else
#define EB_ALLIGN_MALLOC(type, pointer, n_elements, pointer_class) \
if (posix_memalign((void**)(&(pointer)), ALVALUE, n_elements) != 0) { \
    return EB_ErrorInsufficientResources; \
        } \
            else { \
    pointer = (type) pointer;  \
    svt_dec_memory_map[*(svt_dec_memory_map_index)].ptr_type = pointer_class; \
    svt_dec_memory_map[(*(svt_dec_memory_map_index))++].ptr = pointer; \
    if (n_elements % 8 == 0) { \
        *svt_dec_total_lib_memory += (n_elements); \
            } \
            else { \
        *svt_dec_total_lib_memory += ((n_elements) + (8 - ((n_elements) % 8))); \
    } \
} \
if (*(svt_dec_memory_map_index) >= MAX_NUM_PTR) { \
    return EB_ErrorInsufficientResources; \
    } \
svt_dec_lib_malloc_count++;
#endif

#define EB_MALLOC_DEC(type, pointer, n_elements, pointer_class) \
pointer = (type) malloc(n_elements); \
if (pointer == (type)EB_NULL) { \
    return EB_ErrorInsufficientResources; \
    } \
    else { \
    svt_dec_memory_map[*(svt_dec_memory_map_index)].ptr_type = pointer_class; \
    svt_dec_memory_map[(*(svt_dec_memory_map_index))++].ptr = pointer; \
    if (n_elements % 8 == 0) { \
        *svt_dec_total_lib_memory += (n_elements); \
    } \
    else { \
        *svt_dec_total_lib_memory += ((n_elements) + (8 - ((n_elements) % 8))); \
    } \
} \
if (*(svt_dec_memory_map_index) >= MAX_NUM_PTR) { \
    return EB_ErrorInsufficientResources; \
} \
svt_dec_lib_malloc_count++;

EbErrorType dec_mem_init(EbDecHandle  *dec_handle_ptr);

#ifdef __cplusplus
    }
#endif
#endif // EbDecMemInit_h