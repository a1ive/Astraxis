"""Minimal reader for NAIF DAF files (the container of SPK and CK kernels).

Only the Python standard library is used.
"""

import struct


class Daf:
    """The comment area and the segments of a DAF file: for each segment its
    double and integer summary components and its data words."""

    def __init__(self, path, idword):
        with open(path, 'rb') as f:
            data = f.read()
        if not data.startswith(idword):
            raise ValueError('%s: not a %s file' % (path, idword.decode('ascii').strip()))
        end = '<' if data[88:96].startswith(b'LTL') else '>'
        nd, ni = struct.unpack(end + 'ii', data[8:16])
        fward = struct.unpack(end + 'i', data[76:80])[0]
        # The comment area: records 2 .. fward - 1, lines ended by NUL.
        self.comments = data[1024:(fward - 1) * 1024].replace(b'\0', b'\n').decode('latin-1')
        self.segments = []  # (doubles, integers, words), in file order
        record = fward
        summary_size = nd + (ni + 1) // 2
        while record:
            base = (record - 1) * 1024
            nxt, _, count = struct.unpack(end + 'ddd', data[base:base + 24])
            for i in range(int(count)):
                off = base + 24 + i * summary_size * 8
                dc = struct.unpack(end + '%dd' % nd, data[off:off + nd * 8])
                ic = struct.unpack(end + '%di' % ni, data[off + nd * 8:off + nd * 8 + ni * 4])
                first, last = ic[-2:]
                words = struct.unpack(end + '%dd' % (last - first + 1), data[(first - 1) * 8:last * 8])
                self.segments.append((dc, ic, words))
            record = int(nxt)
