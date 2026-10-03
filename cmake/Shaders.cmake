# HLSL -> DXIL (Direct3D 12) and SPIR-V (Vulkan) with dxc, embedded into the
# executable as C++ headers.
#
# File naming: <name>.<stage>.hlsl, where stage is vert | frag | comp.
# Entry point is always `main`.
# Each shader produces the header <binary_dir>/generated/shaders/<name>.<stage>.h
# defining `inline constexpr unsigned char k<Name><Stage>Dxil[]` and
# `k<Name><Stage>Spirv[]`.
#
# dxc is Microsoft's official release, pinned below and fetched at configure
# time (the Windows SDK's dxc cannot emit SPIR-V). Set ASTRAXIS_DXC to use
# another dxc with SPIR-V support (e.g. on hosts without a prebuilt release).

set(ASTRAXIS_DXC "" CACHE FILEPATH "dxc with SPIR-V support; empty = fetch the pinned release")

if(ASTRAXIS_DXC)
    set(ASTRAXIS_DXC_EXECUTABLE "${ASTRAXIS_DXC}")
else()
    # DirectXShaderCompiler v1.9.2609 (2026-09-29), SHA-256 from the GitHub release.
    set(_dxc_base "https://github.com/microsoft/DirectXShaderCompiler/releases/download/v1.9.2609")
    if(CMAKE_HOST_WIN32)
        set(_dxc_url "${_dxc_base}/dxc_2026_09_29.zip")
        set(_dxc_sha256 ad31b1fc8443175d204f77a611fdb3ef2ec42759bdc2f1167368de24a4a7e7f1)
        set(_dxc_exe bin/x64/dxc.exe)
    elseif(CMAKE_HOST_SYSTEM_NAME STREQUAL "Linux" AND CMAKE_HOST_SYSTEM_PROCESSOR MATCHES "^(x86_64|AMD64)$")
        set(_dxc_url "${_dxc_base}/linux_dxc_2026_09_28.x86_x64.tar.gz")
        set(_dxc_sha256 96faadc7f5c282d2ffda49804beb4c3ee38127bc252b723234e3c5cdf7aa39a1)
        set(_dxc_exe bin/dxc) # finds lib/libdxcompiler.so through its $ORIGIN/../lib rpath
    else()
        message(FATAL_ERROR "No prebuilt dxc for ${CMAKE_HOST_SYSTEM_NAME}/${CMAKE_HOST_SYSTEM_PROCESSOR}; "
                            "set ASTRAXIS_DXC to a dxc with SPIR-V support")
    endif()

    include(FetchContent)
    FetchContent_Declare(dxc
        URL      "${_dxc_url}"
        URL_HASH SHA256=${_dxc_sha256}
        DOWNLOAD_EXTRACT_TIMESTAMP TRUE
    )
    FetchContent_MakeAvailable(dxc) # no CMakeLists.txt inside: only downloads and extracts
    set(ASTRAXIS_DXC_EXECUTABLE "${dxc_SOURCE_DIR}/${_dxc_exe}")
endif()

if(NOT EXISTS "${ASTRAXIS_DXC_EXECUTABLE}")
    message(FATAL_ERROR "dxc not found at ${ASTRAXIS_DXC_EXECUTABLE}")
endif()
message(STATUS "dxc: ${ASTRAXIS_DXC_EXECUTABLE}")

set(_ASTRAXIS_EMBED_SCRIPT "${CMAKE_CURRENT_LIST_DIR}/EmbedBinary.cmake")

function(astraxis_add_shaders target)
    set(gen_dir "${CMAKE_CURRENT_BINARY_DIR}/generated")
    set(headers "")
    # Shared include files: any change recompiles every shader.
    file(GLOB includes CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/shaders/*.hlsli")

    foreach(src IN LISTS ARGN)
        get_filename_component(file_name "${src}" NAME)
        if(NOT file_name MATCHES "^([A-Za-z0-9_]+)\\.(vert|frag|comp)\\.hlsl$")
            message(FATAL_ERROR "Shader '${src}' must be named <name>.<vert|frag|comp>.hlsl")
        endif()
        set(name "${CMAKE_MATCH_1}")
        set(stage "${CMAKE_MATCH_2}")

        if(stage STREQUAL "vert")
            set(profile vs_6_0)
        elseif(stage STREQUAL "frag")
            set(profile ps_6_0)
        else()
            set(profile cs_6_0)
        endif()

        # Symbol prefix: k + PascalCase(name) + PascalCase(stage); suffixed Dxil / Spirv.
        string(REPLACE "_" ";" parts "${name};${stage}")
        set(symbol "k")
        foreach(part IN LISTS parts)
            string(SUBSTRING "${part}" 0 1 head)
            string(SUBSTRING "${part}" 1 -1 tail)
            string(TOUPPER "${head}" head)
            string(APPEND symbol "${head}${tail}")
        endforeach()

        set(src_abs "${CMAKE_CURRENT_SOURCE_DIR}/${src}")
        set(dxil "${gen_dir}/shaders/${name}.${stage}.dxil")
        set(spirv "${gen_dir}/shaders/${name}.${stage}.spv")
        set(header "${gen_dir}/shaders/${name}.${stage}.h")
        set(common_args -nologo -T ${profile} -E main -O3 -I "${CMAKE_CURRENT_SOURCE_DIR}/shaders")

        add_custom_command(
            OUTPUT "${header}"
            COMMAND "${CMAKE_COMMAND}" -E make_directory "${gen_dir}/shaders"
            COMMAND "${ASTRAXIS_DXC_EXECUTABLE}" ${common_args} -Fo "${dxil}" "${src_abs}"
            COMMAND "${ASTRAXIS_DXC_EXECUTABLE}" ${common_args} -spirv -fspv-target-env=vulkan1.0 -Fo "${spirv}" "${src_abs}"
            COMMAND "${CMAKE_COMMAND}" -DDXIL=${dxil} -DSPIRV=${spirv} -DOUTPUT=${header} -DSYMBOL=${symbol}
                    -P "${_ASTRAXIS_EMBED_SCRIPT}"
            DEPENDS "${src_abs}" ${includes} "${_ASTRAXIS_EMBED_SCRIPT}"
            COMMENT "dxc ${file_name}"
            VERBATIM
        )
        list(APPEND headers "${header}")
    endforeach()

    # Listed for IDE visibility only; stop Visual Studio from running its own FXC on them.
    set_source_files_properties(${ARGN} ${includes} PROPERTIES VS_TOOL_OVERRIDE "None")
    add_custom_target(${target}_build DEPENDS ${headers} SOURCES ${ARGN} ${includes})
    add_library(${target} INTERFACE)
    add_dependencies(${target} ${target}_build)
    target_include_directories(${target} INTERFACE "${gen_dir}")
endfunction()
