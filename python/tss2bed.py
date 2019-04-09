import os, sys, re

for fn in sys.argv[1:]:
    if fn.endswith('.fa') or fn.endswith('.mfa'):
        filename_out = fn[0:fn.rfind('.')] + '.bed'
    print(filename_out)
    with open(fn) as fi, open(filename_out, 'w') as fo:
        for line in fi:
            if line.startswith('>'):
                m = re.search('(chr\\w+):(.):(\\d+)\\D(\\d+)', line)
                if m:
                    fo.write('{}\t{}\t{}\n'.format(m.group(1), m.group(3), m.group(4)))
