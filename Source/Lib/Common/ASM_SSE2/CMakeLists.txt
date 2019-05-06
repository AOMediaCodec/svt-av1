#
# Copyright(c) 2019 Intel Corporation
# SPDX - License - Identifier: BSD - 2 - Clause - Patent
#

# ASM_SSE2 Directory CMakeLists.txt

# Include Encoder Subdirectories
include_directories(${PROJECT_SOURCE_DIR}/Source/API/)
include_directories(${PROJECT_SOURCE_DIR}/Source/Lib/Common/Codec/)
include_directories(${PROJECT_SOURCE_DIR}/Source/Lib/Common/C_DEFAULT/)
include_directories(${PROJECT_SOURCE_DIR}/Source/Lib/Common/ASM_SSE2/)

if(NOT (MSVC OR (WIN32 AND CMAKE_C_COMPILER_ID STREQUAL "INTEL")))
    if(CMAKE_C_COMPILER_ID STREQUAL "Intel")
        set(CMAKE_C_FLAGS "-static-intel -w")
    else()
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -msse2")
    endif()
endif()

file(GLOB all_files
    "*.h"
    "*.asm"
    "*.c")

add_library(COMMON_ASM_SSE2
    ${all_files})