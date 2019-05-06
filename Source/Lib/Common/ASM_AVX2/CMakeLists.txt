#
# Copyright(c) 2019 Intel Corporation
# SPDX - License - Identifier: BSD - 2 - Clause - Patent
#

# Common/ASM_AVX2 Directory CMakeLists.txt

# Include Encoder Subdirectories
include_directories(${PROJECT_SOURCE_DIR}/Source/API/)
include_directories(${PROJECT_SOURCE_DIR}/Source/Lib/Common/Codec/)
include_directories(${PROJECT_SOURCE_DIR}/Source/Lib/Common/C_DEFAULT/)
include_directories(${PROJECT_SOURCE_DIR}/Source/Lib/Common/ASM_SSE2/)
include_directories(${PROJECT_SOURCE_DIR}/Source/Lib/Common/ASM_SSSE3/)
include_directories(${PROJECT_SOURCE_DIR}/Source/Lib/Common/ASM_SSE4_1/)
include_directories(${PROJECT_SOURCE_DIR}/Source/Lib/Common/ASM_AVX2/)
link_directories(${PROJECT_SOURCE_DIR}/Source/Lib/Common/ASM_SSSE3/)

if(MSVC OR (WIN32 AND CMAKE_C_COMPILER_ID STREQUAL "INTEL"))
    # Intel Windows
    #(*Note - The Warning level /W0 should be made to /W4 at some point)
    if(CMAKE_C_COMPILER_ID MATCHES "Intel")
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /W0 /Qwd10148 /Qwd10010 /Qwd10157")
    else()
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /arch:AVX2")
    endif()
else()
    if(CMAKE_C_COMPILER_ID STREQUAL "Intel")
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -static-intel -w")
    else()
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -mavx2")
    endif()
endif()

file(GLOB all_files
    "*.h"
    "*.asm"
    "*.c")

add_library(COMMON_ASM_AVX2
    ${all_files})