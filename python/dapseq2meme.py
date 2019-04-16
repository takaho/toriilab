import os, sys, re

with open(sys.argv[2], 'w') as fo:
    fo.write("""MEME version 4.4

ALPHABET= ACGT

strands: + -

Background letter frequencies (from /mnt/thumper-e1/home/shhuang/projects/dap/analysis.v4/gem07_rep_memechip03/FHA_tnt/FHA2_col_b/background):
A 0.30800 C 0.19200 G 0.19200 T 0.30800

""")
    for path, dns, fns in os.walk(sys.argv[1]):
        for fn in fns:
            if fn == 'meme_m1.txt':
                sys.stderr.write(path + '\n')
                with open(os.path.join(path, fn)) as fi:
                    contents = fi.readlines()
                    for i, l in enumerate(contents):
                        if l.startswith('MOTIF'):
                            for ln in contents[i:]:
                                fo.write(ln)
                            fo.write('\n')

                        