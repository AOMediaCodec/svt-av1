#
# Copyright(c) 2018 Intel Corporation
# SPDX - License - Identifier: BSD - 2 - Clause - Patent
#

cmake_minimum_required(VERSION 3.5)
set(project_name "svt-av1")

project(${project_name} C CXX ASM_NASM)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${CMAKE_BUILD_TYPE}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${CMAKE_BUILD_TYPE}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${CMAKE_BUILD_TYPE}/bin)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_C_STANDARD 99)
set(CMAKE_CXX_STANDARD 11)
set(CAN_USE_ASSEMBLER TRUE)
set(CMAKE_NASM_LINK_EXECUTABLE "ld -o <TARGET> <OBJECTS>")
set(CMAKE_ASM_NASM_SOURCE_FILE_EXTENSIONS asm)
find_program(YASM_EXE yasm)
if(YASM_EXE)
    set(CMAKE_ASM_NASM_COMPILER YASM_EXE CACHE PATH "yasm path")
endif()

if(UNIX)
    set(CMAKE_ASM_NASM_FLAGS "-DUNIX64")
else()
    set(CMAKE_ASM_NASM_FLAGS "-DWIN64")
endif()

if(APPLE)
    set(CMAKE_MACOSX_RPATH 1)
    set(CMAKE_ASM_NASM_OBJECT_FORMAT macho64)
    set(CMAKE_C_ARCHIVE_CREATE "<CMAKE_AR> Scr <TARGET> <LINK_FLAGS> <OBJECTS>")
    set(CMAKE_CXX_ARCHIVE_CREATE "<CMAKE_AR> Scr <TARGET> <LINK_FLAGS> <OBJECTS>")
    set(CMAKE_C_ARCHIVE_FINISH "<CMAKE_RANLIB> -no_warning_for_no_symbols -c <TARGET>")
    set(CMAKE_CXX_ARCHIVE_FINISH "<CMAKE_RANLIB> -no_warning_for_no_symbols -c <TARGET>")
elseif(UNIX)
    set(CMAKE_ASM_NASM_OBJECT_FORMAT elf64)
else()
    set(CMAKE_ASM_NASM_OBJECT_FORMAT win64)
endif()

if(MSVC OR (WIN32 AND CMAKE_C_COMPILER_ID STREQUAL "INTEL"))
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /MP")
    set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} /Od")
    set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} /O2")
    set(CMAKE_ASM_NASM_FLAGS_DEBUG "${CMAKE_ASM_NASM_FLAGS_DEBUG} -g -F cv8")
    set(CMAKE_EXE_LINKER_FLAGS_DEBUG "${CMAKE_EXE_LINKER_FLAGS_DEBUG} /DEBUG")
else()
    set(CMAKE_C_FLAGS
        "${CMAKE_C_FLAGS} -Wall -Wextra -Wformat -Wformat-security -mavx -fstack-protector-strong")
    if(NOT MINGW)
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC -fPIE")
    endif()
    set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -O2 -D_FORTIFY_SOURCE=2")
    set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -O0 -g")
    if(NOT MINGW)
        set(CMAKE_EXE_LINKER_FLAGS
            "${CMAKE_EXE_LINKER_FLAGS} -pie")
    endif()
    if(UNIX AND NOT APPLE)
        set(CMAKE_EXE_LINKER_FLAGS
            "${CMAKE_EXE_LINKER_FLAGS} -pie -z noexecstack -z relro -z now")
    endif()
    set(CMAKE_SHARED_LINKER_FLAGS
        "${CMAKE_SHARED_LINKER_FLAGS} ${CMAKE_EXE_LINKER_FLAGS}")
endif()

option(BUILD_TESTING "" OFF)
include(CTest)

set(COVERAGE false CACHE BOOL "For use with coveralls and GTest")

# Prepare for Coveralls
if(COVERAGE)
    if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
        set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} --coverage")
        set(CMAKE_EXE_LINKER_FLAGS_DEBUG "${CMAKE_EXE_LINKER_FLAGS_DEBUG} --coverage")
    endif()
endif()

#set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/CMakeModules)
#if(USE_CODE_COVERAGE AND CMAKE_COMPILER_IS_GNUCXX)
#    include(CodeCoverage)
#    setup_target_for_coverage(${PROJECT_NAME}_coverage ${PROJECT_TEST_NAME} coverage)
#endif()

# Add Subdirectories
add_subdirectory(Source/Lib/Common)
add_subdirectory(Source/Lib/Encoder)
add_subdirectory(Source/App/EncApp)
add_subdirectory(Source/Lib/Decoder)
add_subdirectory(Source/App/DecApp)
if(BUILD_TESTING)
    add_subdirectory(test)
    add_subdirectory(third_party/googletest)
endif()
