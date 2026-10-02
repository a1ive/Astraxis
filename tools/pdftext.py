"""Crude PDF text extraction with the standard library only (Flate streams, Tj/TJ operators)."""
import re
import sys
import zlib


def streams(data):
    for m in re.finditer(rb'stream\r?\n', data):
        start = m.end()
        end = data.find(b'endstream', start)
        if end < 0:
            continue
        raw = data[start:end]
        try:
            yield zlib.decompress(raw)
        except zlib.error:
            try:
                yield zlib.decompressobj().decompress(raw)
            except zlib.error:
                pass


def unescape(s):
    out = []
    i = 0
    while i < len(s):
        c = s[i]
        if c == '\\' and i + 1 < len(s):
            n = s[i + 1]
            if n in '()\\':
                out.append(n)
                i += 2
                continue
            if n.isdigit():
                j = i + 1
                while j < len(s) and j < i + 4 and s[j].isdigit():
                    j += 1
                out.append(chr(int(s[i + 1:j], 8)))
                i = j
                continue
            out.append(n)
            i += 2
            continue
        out.append(c)
        i += 1
    return ''.join(out)


def text_of(content):
    s = content.decode('latin-1')
    lines = []
    for block in re.findall(r'BT(.*?)ET', s, re.S):
        parts = []
        for tj in re.finditer(r'\[(.*?)\]\s*TJ|\((.*?)(?<!\\)\)\s*Tj|(T\*|Td|TD|Tm)', block, re.S):
            if tj.group(1) is not None:
                seg = ''
                for item in re.finditer(r'\((.*?)(?<!\\)\)|(-?\d+\.?\d*)', tj.group(1), re.S):
                    if item.group(1) is not None:
                        seg += unescape(item.group(1))
                    elif item.group(2) and float(item.group(2)) < -200:
                        seg += ' '
                parts.append(seg)
            elif tj.group(2) is not None:
                parts.append(unescape(tj.group(2)))
            else:
                parts.append(' ')
        lines.append(''.join(parts))
    return '\n'.join(lines)


data = open(sys.argv[1], 'rb').read()
out = []
for st in streams(data):
    if b'BT' in st and (b'TJ' in st or b'Tj' in st):
        out.append(text_of(st))
text = '\n'.join(out)
sys.stdout.reconfigure(encoding='utf-8')
print(text)
