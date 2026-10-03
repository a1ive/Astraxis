# Script mode: cmake -DDXIL=<file> -DSPIRV=<file> -DOUTPUT=<header> -DSYMBOL=<prefix> -P EmbedBinary.cmake
# Writes a header with `inline constexpr unsigned char <prefix>Dxil[]` and
# `<prefix>Spirv[]` holding the bytes of the two shader binaries.

function(embed_array input name out_var)
    file(READ "${input}" hex HEX)
    string(LENGTH "${hex}" hex_len)
    math(EXPR size "${hex_len} / 2")
    string(REGEX REPLACE "([0-9a-f][0-9a-f])" "0x\\1," bytes "${hex}")
    # Break into lines of 16 bytes for readability.
    string(REGEX REPLACE "((0x[0-9a-f][0-9a-f],){16})" "\\1\n    " bytes "${bytes}")
    set(${out_var} "inline constexpr unsigned char ${name}[${size}] = {\n    ${bytes}\n};\n" PARENT_SCOPE)
endfunction()

embed_array("${DXIL}" "${SYMBOL}Dxil" dxil_array)
embed_array("${SPIRV}" "${SYMBOL}Spirv" spirv_array)

file(WRITE "${OUTPUT}"
"// Generated from ${DXIL} and ${SPIRV}. Do not edit.
#pragma once

${dxil_array}
${spirv_array}")
