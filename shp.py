"""Lecteur minimal de shapefile (polylignes, polygones, points) + DBF, sans dépendance externe."""
import struct


def read_dbf(path):
    with open(path, "rb") as f:
        data = f.read()
    nrec = struct.unpack("<I", data[4:8])[0]
    hlen, rlen = struct.unpack("<HH", data[8:12])
    fields = []
    pos = 32
    while data[pos] != 0x0D:
        name = data[pos:pos + 11].split(b"\x00")[0].decode("ascii")
        ftype = chr(data[pos + 11])
        flen = data[pos + 16]
        fields.append((name, ftype, flen))
        pos += 32
    rows = []
    for i in range(nrec):
        rec = data[hlen + i * rlen: hlen + (i + 1) * rlen]
        if not rec or rec[0] == 0x2A:  # enregistrement supprimé
            rows.append(None)
            continue
        off = 1
        row = {}
        for name, ftype, flen in fields:
            raw = rec[off:off + flen]
            off += flen
            s = raw.decode("utf-8", errors="replace").strip()
            if ftype in "NF":
                try:
                    row[name] = int(s) if ftype == "N" and "." not in s else float(s)
                except ValueError:
                    row[name] = None
            else:
                row[name] = s
        rows.append(row)
    return rows


def read_shp(path):
    """Retourne une liste de géométries : dict(type, parts=[[(x,y),...],...])."""
    with open(path, "rb") as f:
        data = f.read()
    stype = struct.unpack("<i", data[32:36])[0]
    pos = 100
    geoms = []
    while pos + 8 <= len(data):
        _, clen = struct.unpack(">ii", data[pos:pos + 8])
        pos += 8
        rec = data[pos:pos + clen * 2]
        pos += clen * 2
        t = struct.unpack("<i", rec[:4])[0]
        if t == 0:
            geoms.append(None)
            continue
        if t in (3, 5):  # polyligne / polygone
            npart, npts = struct.unpack("<ii", rec[36:44])
            parts = list(struct.unpack("<%di" % npart, rec[44:44 + 4 * npart]))
            ptoff = 44 + 4 * npart
            pts = struct.unpack("<%dd" % (2 * npts), rec[ptoff:ptoff + 16 * npts])
            xy = [(pts[2 * i], pts[2 * i + 1]) for i in range(npts)]
            parts.append(npts)
            geoms.append({"type": "line" if t == 3 else "polygon",
                          "parts": [xy[parts[k]:parts[k + 1]] for k in range(npart)]})
        elif t == 1:
            x, y = struct.unpack("<dd", rec[4:20])
            geoms.append({"type": "point", "parts": [[(x, y)]]})
        else:
            geoms.append(None)
    return geoms
