import re
import os
import sys

kvpat = re.compile(r'(.*?)=(".*?"|.*?)(?:,|$)')
def parse_line(line):
    s = line.index('<')
    e = line.index('>')
    return dict(kvpat.findall(line[s+1:e]))

def parse_kv(kv):
    kv = kv.split('=',1)
    return len(kv) == 2 and kv or [kv[0], True]

infos, fmts = {}, {}
for line in sys.stdin:
    if line[0] == '#':
        if line.startswith('##INFO'):
            d = parse_line(line)
            if d['Number'] in 'ARG':
                infos[d['ID']] = d['Number']
        if line.startswith('##FORMAT'):
            d = parse_line(line)
            if d['Number'] in 'ARG':
                fmts[d['ID']] = d['Number']
        sys.stdout.write(line)
        continue

    vals = line.rstrip().split('\t')
    if len(vals) < 9 or vals[8] == '.':
        # pop vcf only site, skip
        continue

    ka, kr, kg, info = None, None, None, None
    if len(vals) >= 8 and vals[7] != '.':
        info = list(map(parse_kv, vals[7].split(';')))
        for k,v in info:
            if k == 'AF':
                ka = [x != '.' for x in v.split(',')]
                kr = [True] + ka
                break

    if ka is None:
        # no-call, drop
        continue

    if all(ka):
        sys.stdout.write(line)
        continue

    if vals[4] != '.':
        v = vals[4].split(',')
        v = [x for x,y in zip(v, ka) if y]
        r = vals[3]
        #n = len(os.path.commonprefix([r[:0:-1]] + [a[:0:-1] for a in v]))
        n = min(len(r), min(map(len, v))) - 1
        if n > 0:
            r, v = r[:-n], [a[:-n] for a in v]
        vals[3] = r
        vals[4] = ','.join(v)

    if info:
        for i,(k,v) in enumerate(info):
            n = infos.get(k) if v != '.' else None
            if n == 'A':
                v = v.split(',')
                info[i] = k + '=' + ','.join(x for x,y in zip(v, ka) if y)
            elif n == 'R':
                v = v.split(',')
                info[i] = k + '=' + ','.join(x for x,y in zip(v, kr) if y)
            elif v is True:
                info[i] = k
            else:
                info[i] = k + '=' + v
        vals[7] = ';'.join(info)

    if len(vals) >= 10 and vals[9] != '.':
        fmt = vals[8].split(':')
        kg = [x and y for i,x in enumerate(kr) for y in kr[:i+1]]
        for j,val in enumerate(vals[9:]):
            val = val.split(':')
            for i,v in enumerate(val):
                n = fmts.get(fmt[i]) if v != '.' else None
                if n == 'A':
                    v = v.split(',')
                    val[i] = ','.join(x for x,y in zip(v, ka) if y)
                elif n == 'R':
                    v = v.split(',')
                    val[i] = ','.join(x for x,y in zip(v, kr) if y)
                elif n == 'G':
                    v = v.split(',')
                    val[i] = ','.join(x for x,y in zip(v, kg) if y)
                # GT not affected, new alt alleles are added at the end
            vals[9+j] = ':'.join(val)

    line = '\t'.join(vals) + '\n'
    sys.stdout.write(line)
