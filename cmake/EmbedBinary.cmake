# Script mode: cmake -DINPUT=<file> -DOUTPUT=<header> -DSYMBOL=<name> -P EmbedBinary.cmake
# Writes a header with `inline constexpr unsigned char SYMBOL[]` holding the file bytes.

file(READ "${INPUT}" hex HEX)
string(LENGTH "${hex}" hex_len)
math(EXPR size "${hex_len} / 2")

string(REGEX REPLACE "([0-9a-f][0-9a-f])" "0x\\1," bytes "${hex}")
# Break into lines of 16 bytes for readability.
string(REGEX REPLACE "((0x[0-9a-f][0-9a-f],){16})" "\\1\n    " bytes "${bytes}")

file(WRITE "${OUTPUT}"
"// Generated from ${INPUT}. Do not edit.
#pragma once

inline constexpr unsigned char ${SYMBOL}[${size}] = {
    ${bytes}
};
")
