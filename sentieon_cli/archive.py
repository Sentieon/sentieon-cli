import os
import sys

def ar_load(path):
    """
    Functionality for reading GNU ar archives.
    With an archive file, list all archive members.
    With an archive/member, extract and return the member as bytes
    """

    base, names = None, []
    while True:
        try:
            fp = open(path, 'rb')
        except NotADirectoryError:
            path, base = os.path.split(path)
        except:
            raise
        else:
            break

    magic = fp.read(8)
    if magic != b'!<arch>\n':
        raise RuntimeError('not an ar file')
    strtab = None
    hdr = fp.read(60)
    while len(hdr) > 0:
        if len(hdr) != 60:
            raise RuntimeError('ar file truncated')
        if hdr[58:] != b'`\n':
            raise RuntimeError('ar file corrupted')
        name = hdr[:16].rstrip(b' ')
        size = int(hdr[48:58].decode())
        pad = size % 2
        if name[:1] == b'/':
            # sysv special files
            if name == b'/':
                # symbol table
                name = None
            elif name == b'//':
                # string table
                strtab = fp.read(size)
                size = 0
                name = None
            else:
                # long name in string table
                idx = int(name[1:].decode())
                end = strtab.find(b'/\n', idx)
                if end == -1:
                    raise RuntimeError('invalid entry')
                name = strtab[idx:end]
        elif name[-1:] == b'/':
            # sysv style
            name = name[:-1]
        elif name[:3] == b'#1/':
            # bsd style long name
            nlen = int(name[3:].decode())
            name = fp.read(nlen)
            size -= nlen
        else:
            # bsd style
            pass
        if name and base == name.decode():
            return fp.read(size)
        if name and base is None:
            names.append(name.decode())
        fp.seek(size + pad, 1)
        hdr = fp.read(60)
    return names

if __name__ == '__main__':
    print(ar_load(sys.argv[1]))