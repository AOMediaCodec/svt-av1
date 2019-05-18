#!/bin/bash
#
# Copyright(c) 2019 Intel Corporation
# SPDX - License - Identifier: BSD - 2 - Clause - Patent
#

function cmake_ify {
  
    PATH=$PATH:/usr/local/bin/
    cmake ..                                        \
        -DCMAKE_ASM_NASM_COMPILER=$CMAKE_ASSEMBLER  \
        -G Xcode                                    \
        
}

# Defines
CMAKE_ASSEMBLER=yasm

mkdir _build && cd _build

cmake_ify

exit
