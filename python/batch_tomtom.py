import os, sys, re

"""
dreme.txt
Background letter frequencies (from dataset):
A 0.341 C 0.159 G 0.159 T 0.340


MOTIF GGCCCA DREME-1

#             Word    RC Word        Pos        Neg    P-value    E-value
# BEST      GGCCCA     TGGGCC       2874       2797   1.9e-549   7.2e-543
#           GGCCCA     TGGGCC       2874       2797   1.9e-549   7.2e-543

letter-probability matrix: alength= 4 w= 6 nsites= 4420 E= 7.2e-543
0.000000 0.000000 1.000000 0.000000
0.000000 0.000000 1.000000 0.000000
0.000000 1.000000 0.000000 0.000000
0.000000 1.000000 0.000000 0.000000
0.000000 1.000000 0.000000 0.000000
1.000000 0.000000 0.000000 0.000000
"""

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-i', nargs='+')
parser.add_argument('-o', default='out')
parser.add_argument('--verbose', action='store_true')
args = parser.parse_args()

class MEMEMotif(object):
    def __init__(self, name, length, n_sites=0, score=0, matrix=None):
        self.name = name
        self.length = length
        self.n_sites = n_sites
        self.score = score
        self.set_matrix(matrix)
        # self.matrix = matrix
    def set_matrix(self, mat):
        if mat is None:
            return
        if len(mat) != self.length:
            raise Exception('matrix shape should be [{},4]'.format(self.length))
        for r in mat:
            if len(r) != 4:
                raise Exception('matrix shape should be [{},4]'.format(self.length))
        self.matrix = mat
    def __repr__(self):
        ostr = '{} Score:{:.2e} Sites:{}\n'.format(self.name, self.score, self.n_sites)
        for i in range(4):
            line = ''
            for j in range(self.length):
                if j > 0: line += ' '
                line += '{:.3f}'.format(self.matrix[j][i])
            ostr += line + '\n'
        return ostr
            
for fn in args.i:
    if os.path.isdir(fn):
        name = os.path.basename(fn)
        fn = os.path.join(fn, 'dreme.txt')
    elif os.path.isfile(fn):
        name = os.path.basename(os.path.dirname(fn))
    else:
        sys.stderr('Skip {} : directory or dreme.txt is acceptable\n'.format(fn))
        continue
    if os.path.exists(fn) is False:
        sys.stderr.write('no {}\n'.format(fn))
        continue
    with open(fn) as fi:
        contents = fi.readlines()
    # contents = open(fn).readlines().split('\t')
    num_lines = len(contents)
    i = 0
    motifs = []
    background = [.25, .25, .25, .25]
    while i < num_lines:
        line = contents[i]
        # print(line)
        if line.startswith('Background letter frequencies'):
            m = re.match('A\\s+([\\d\\.]+)\\s+C\\s+([\\d\\.]+)\\s+G\\s+([\\d\\.]+)\\s+T\\s+([\\d\\.]+)\\s+', contents[i+1])
            if m:
                acgt = [float(x_) for x_ in m.groups()]
                background = [float(x_) / sum(acgt) for x in acgt]
            i += 1
        elif line.startswith('MOTIF'):
            m = re.match('MOTIF\\s+(\\w+)', line)
            name = m.group(1)
            j = i + 1
            while j < num_lines:
                m = re.search('alength=\\s*(\\d+) w=\\s*(\\d+) nsites=\\s*(\\d+) E=\\s*([\\d\\.e\\-]+)', contents[j]) #7.2e-543', info)
                if m:
                    l = int(m.group(1))
                    w = int(m.group(2))
                    n = int(m.group(3))
                    score = float(m.group(4))
                    j = j + 1
                    break
                j += 1
            # print(name)
            # for l in contents[i:i+7]:
            #     print(l)
            # info = contents[i + 6]
            # m = re.search('alength=\\s*(\\d+) w=\\s*(\\d+) nsites=\\s*(\\d+) E=\\s*([\\d\\.e\\-]+)', info) #7.2e-543', info)
            # l = int(m.group(1))
            # w = int(m.group(2))
            # n = int(m.group(3))
            # score = float(m.group(4))
            # j = i + 7
            # print(l, w, n)
            # print(contents[j:j+w+1])
            matrix = []
            for k in range(w):
            # while j < num_lines:#min(i + 7 + w, num_lines):
                # print(contents[j+k])
                freq = [float(x_) for x_ in re.split('\\s+', contents[j+k].strip())]
                matrix.append(freq)
                k += 1
            motif = MEMEMotif(name, w, n, score, matrix)
            print(str(motif))
            motifs.append(motif)
            i = j + 1
        else:
            i += 1