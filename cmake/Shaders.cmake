# HLSL -> DXIL compilation with dxc, embedded into the executable as C++ headers.
#
# File naming: <name>.<stage>.hlsl, where stage is vert | frag | comp.
# Entry point is always `main`.
# Each shader produces the header <binary_dir>/generated/shaders/<name>.<stage>.h
# defining `inline constexpr unsigned char k<Name><Stage>Dxil[]`.

# Locate dxc: prefer the Vulkan SDK (if any), else the newest Windows SDK.
set(_dxc_hints "")
if(DEFINED ENV{VULKAN_SDK})
    list(APPEND _dxc_hints "$ENV{VULKAN_SDK}/Bin")
endif()
file(GLOB _winsdk_dxc_dirs "C:/Program Files (x86)/Windows Kits/10/bin/10.*/x64")
list(SORT _winsdk_dxc_dirs COMPARE NATURAL ORDER DESCENDING)
list(APPEND _dxc_hints ${_winsdk_dxc_dirs})

find_program(DXC_EXECUTABLE dxc HINTS ${_dxc_hints} REQUIRED)
message(STATUS "dxc: ${DXC_EXECUTABLE}")

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

        # Symbol: k + PascalCase(name) + PascalCase(stage) + Dxil
        string(REPLACE "_" ";" parts "${name};${stage}")
        set(symbol "k")
        foreach(part IN LISTS parts)
            string(SUBSTRING "${part}" 0 1 head)
            string(SUBSTRING "${part}" 1 -1 tail)
            string(TOUPPER "${head}" head)
            string(APPEND symbol "${head}${tail}")
        endforeach()
        string(APPEND symbol "Dxil")

        set(src_abs "${CMAKE_CURRENT_SOURCE_DIR}/${src}")
        set(dxil "${gen_dir}/shaders/${name}.${stage}.dxil")
        set(header "${gen_dir}/shaders/${name}.${stage}.h")

        add_custom_command(
            OUTPUT "${header}"
            COMMAND "${CMAKE_COMMAND}" -E make_directory "${gen_dir}/shaders"
            COMMAND "${DXC_EXECUTABLE}" -nologo -T ${profile} -E main -O3
                    -I "${CMAKE_CURRENT_SOURCE_DIR}/shaders"
                    -Fo "${dxil}" "${src_abs}"
            COMMAND "${CMAKE_COMMAND}" -DINPUT=${dxil} -DOUTPUT=${header} -DSYMBOL=${symbol}
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
