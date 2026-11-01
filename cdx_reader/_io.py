"""Narrow binary compatibility repairs before the unchanged pycdxml decoder."""
from io import BytesIO
from pathlib import Path
from pycdxml import cdxml_converter


def normalize_integer_widths(blob):
    """Only lossless widths of Charge and graphic enums; never drop a property.

    CDX is a tag/length stream, not an XML-like text. Walk actual records so a
    byte sequence inside text/image payload cannot be mistaken for a property.
    Unknown tags are preserved; the normal decoder remains responsible for them.
    """
    if len(blob) < 28 or not blob.startswith(b'VjCD0100') or blob[22:24] != b'\x00\x80':
        return blob, []
    out = bytearray(blob[:22]); pos = 22; stack = []; repairs = []
    while pos < len(blob):
        if pos + 2 > len(blob):
            raise ValueError('truncated CDX tag')
        start = pos; tag = int.from_bytes(blob[pos:pos+2],'little'); pos += 2
        if tag == 0:
            if stack: stack.pop()
            out.extend(blob[start:pos]); continue
        if tag & 0x8000:
            if pos + 4 > len(blob):
                raise ValueError('truncated CDX object ID')
            stack.append((tag,int.from_bytes(blob[pos:pos+4],'little')))
            pos += 4; out.extend(blob[start:pos]); continue
        if pos + 2 > len(blob):
            raise ValueError('truncated CDX property size')
        size = int.from_bytes(blob[pos:pos+2],'little'); pos += 2
        if size == 0xffff:
            if pos + 4 > len(blob):
                raise ValueError('truncated extended CDX property size')
            size = int.from_bytes(blob[pos:pos+4],'little'); pos += 4
        end = pos + size
        if end > len(blob):
            raise ValueError('truncated CDX property payload')
        payload = blob[pos:end]
        graphic = {0x0a00:'GraphicType',0x0a01:'LineType',0x0a02:'ArrowType',
                   0x0a03:'RectangleType',0x0a06:'BracketType',0x0a07:'SymbolType'}
        width = 1 if tag == 0x0421 and size in (2,4) else 2 if tag in graphic and size in (1,4) else None
        if width:
            signed = tag == 0x0421
            value = int.from_bytes(payload,'little',signed=signed)
            try:
                canonical = value.to_bytes(width,'little',signed=signed)
            except OverflowError:
                canonical = None
            if canonical is not None:
                out.extend(tag.to_bytes(2,'little') + width.to_bytes(2,'little') + canonical)
                repairs.append(dict(property='Charge' if signed else graphic[tag],
                    object=stack[-1][1] if stack else None,offset=start,old_width=size,width=width,value=value))
            else:
                out.extend(blob[start:end])
        else:
            out.extend(blob[start:end])
        pos = end
    return bytes(out), repairs


def read_document(source):
    """Read a path/bytes/BytesIO without mutating pycdxml's global type registry."""
    blob = source if isinstance(source,bytes) else source.getvalue() if isinstance(source,BytesIO) else Path(source).read_bytes()
    try:
        return cdxml_converter.read_cdx(BytesIO(blob))
    except ValueError as exc:
        if str(exc) not in ('INT8 should consist of exactly 1 byte.',
                            'CDXGraphicType should consist of exactly 2 bytes.',
                            'CDXLineType should consist 2 bytes.',
                            'CDXArrowType should consist of exactly 2 bytes.',
                            'CDXRectangleType should consist 2 bytes.',
                            'CDXBracketType should consist of exactly 2 byte.',
                            'CDXSymbolType should consist of exactly 2 bytes.'):
            raise
        normalized, repairs = normalize_integer_widths(blob)
        if not repairs:
            raise
        doc = cdxml_converter.read_cdx(BytesIO(normalized))
        doc._cdx_binary_repairs = repairs
        return doc
