"""
Aggregate peaks detected using HOMER and MEME (DREME)
"""
import os, sys, re
import pathlib, argparse
import numpy as np
import collections

class MEMEMotif(object):
    def __init__(self, name, length, n_sites=0, score=0, matrix=None):
        self.name = name
        self.length = length
        self.n_sites = n_sites
        self.score = score
        self.set_matrix(matrix)
        self.background = [.25] * 4
        # self.matrix = matrix
    def set_background(self, bg):
        if len(bg) == 4:
            self.background = bg
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
    def output_meme(self, ost):
        a, c, g, t = self.background
        ost.write('MEME version 5\n\nALPHABET= ACGT\n\nstrands: + -\n\nBackground letter frequencies (from uniform background):\nA {:.4f} C {:.4f} G {:.4f} T {:.4f}\n\n'.format(a, c, g, t))
        ost.write('MOTIF {}\n\nletter-probability matrix: alength= 4 w= {} nsites= {} E={:.3e}\n\n'.format(self.name, self.length, self.n_sites, self.score))
        for i in range(4):
            line = ''
            for j in range(self.length):
                if j > 0: line += ' '
                line += '{:.3f}'.format(self.matrix[j][i])
            ost.write(line)
            ost.write('\n')
    @classmethod
    def load_dreme(cls, fn):
        if os.path.isdir(fn):
            title = os.path.basename(fn)
            fn = os.path.join(fn, 'dreme.txt')
        elif os.path.isfile(fn):
            title = os.path.basename(os.path.dirname(fn))
        else:
            raise Exception('Skip {} : directory or dreme.txt is acceptable\n'.format(fn))

        if os.path.exists(fn) is False:
            raise Exception('no {}\n'.format(fn))

        with open(fn) as fi:
            contents = fi.readlines()
        num_lines = len(contents)
        i = 0
        motifs = []
        background = [.25, .25, .25, .25]
        while i < num_lines:
            line = contents[i]
            if line.startswith('Background letter frequencies'):
                m = re.match('A\\s+([\\d\\.]+)\\s+C\\s+([\\d\\.]+)\\s+G\\s+([\\d\\.]+)\\s+T\\s+([\\d\\.]+)\\s+', contents[i+1])
                if m:
                    acgt = [float(x_) for x_ in m.groups()]
                    if len(acgt) == 4:
                        background = [float(x_) / sum(acgt) for x_ in acgt]
                i += 1
            elif line.startswith('MOTIF'):
                m = re.match('MOTIF\\s+(\\w+)', line)
                name = m.group(1)
                start = -1
                j = i + 1
                while j < num_lines:
                    if contents[j].startswith('letter-probability'):
                        start = j
                        break
                    j += 1
                if start < 0:
                    i += 1
                    continue
                info = contents[start]
                m = re.search('alength=\\s*(\\d+) w=\\s*(\\d+) nsites=\\s*(\\d+) E=\\s*([\\d\\.e\\-]+)', info) #7.2e-543', info)
                l = int(m.group(1))
                w = int(m.group(2))
                n = int(m.group(3))
                score = float(m.group(4))
                j = start + 1
                matrix = []
                while j < min(start + 1 + w, num_lines):
                    freq = [float(x_) for x_ in re.split('\\s+', contents[j].strip())]
                    matrix.append(freq)
                    j += 1
                motif = MEMEMotif(name, w, n, score, matrix)
                # motif.output_meme(sys.stdout)
                motifs.append(motif)
                i = j + 1
            else:
                i += 1
        return motifs

class HomerMotif(object):
    def __init__(self, motif, name, matrix, info):
        self.motif = motif
        self.name = name
        self.dbname = None
        self.matrix = matrix
        self.info = info
    def set_knwon_name(self, name):
        if self.dbname is None:
            self.dbname = []
        self.dbname.append(name)

    def __get_enrichment(self):
        if 'T' in self.info and 'B' in self.info:
            t = self.info['T'][1]
            b = self.info['B'][1]
            if t > 0 and b > 0:
                return t / b
        return 1.0

    def __get_hits(self):
        if 'T' in self.info:
            return self.info['T'][0]
        else:
            return 0

    def __get_background(self):
        if 'B' in self.info:
            return self.info['B'][0]
        else:
            return 0

    def __get_known_name(self):
        priority_names = 'SEP1', 'SOC1', 'PI', 'AGL15', 'ATHB-5', 'SEP3', 'TBP3', 'SEP4', 'SEP5',  
        if self.dbname is None:
            return '.'
        for n in priority_names:
            if n in self.dbname:
                return n
        return self.dbname[0]

    enrichment = property(__get_enrichment)
    hits = property(__get_hits)
    background = property(__get_background)
    known = property(__get_known_name)

    def __repr__(self):
        ostr = '{} {}\n'.format(self.motif, self.name)
        for i in range(4):
            line = ''
            for j in range(self.length):
                if j > 0: line += ' '
                line += '{:.3f}'.format(self.matrix[j][i])
            ostr += line + '\n'
        return ostr
    @classmethod
    def convert_to_meme(cls, motifs, filename_template=None):
        if filename_template is None:
            filename_template = os.path.join(os.path.dirname(__file__), 'meme_template.txt')
        with open(filename_template) as fi:
            contents = fi.read()
        ostr = contents
        for i, motif in enumerate(motifs):
            ostr += 'MOTIF {} HOMER-{}\n\n'.format(motif.motif, i + 1)#motif.name)
            ostr += 'letter-probability matrix: alength= 4 w= {} nsites= {} E= {:.4e}\n'.format(motif.length, motif.info.get('T', [0,0])[0], motif.info.get('P', 0))
            for row in motif.matrix:
                ostr += '{:.6f}\t{:.6f}\t{:.6f}\t{:.6f}\n'.format(row[0], row[1], row[2], row[3])
            ostr += '\n\n'
        return ostr

    length = property(lambda s:len(s.matrix))
    @classmethod
    def load_motifs(cls, srcdir):
        motifs = {}
        for fn in os.listdir(srcdir):
            m = re.match('motif(\\d+)\\.motif$', fn)
            if m is None:
                continue
            number = int(m.group(1))
            with open(os.path.join(srcdir, fn)) as fi:
                # items = re.split('\\s+', fi.readline()[1:-1])
                items = fi.readline()[1:-1].split('\t')
                # print(items)
                motif = items[0]
                name = items[1]
                evalue = float(items[2])
                score = float(items[3])
                dig = float(items[4])
                info = {}
                for elem in items[5].split(','):
                    key, val = elem.split(':', 1)
                    if key == 'T' or key == 'B':
                        m = re.match('([\\d\\.]+)\\(([\\d\\.]+)%\\)', val)
                        if m:
                            count, percentage = float(m.group(1)), float(m.group(2))
                            info[key] = count, percentage
                    elif key == 'P':
                        pvalue = float(val)
                info['E'] = evalue
                matrix = []
                for line in fi:
                    items = re.split('\\s+', line.strip())
                    matrix.append([float(x) for x in items])
                motifs[number] = HomerMotif(motif, name, matrix, info)
        return [motifs[n] for n in sorted(motifs.keys())]

def check_match(seq1, seq2):
    import Bio.Seq
    f = seq1
    r = Bio.Seq.reverse_complement(seq1)
    if f.find(seq2) >= 0 or r.find(seq2) >= 0 or seq2.find(f) >= 0 or seq2.find(r) >= 0:
        return True
    return False

parser = argparse.ArgumentParser()
parser.add_argument('-t', nargs='+', default=['AtML1', 'SPCH', 'MUTE', 'FAMA', 'GC1'])
parser.add_argument('-o', default='out')
parser.add_argument('-i', help='source directory')
parser.add_argument('--minimum-hits', type=int, default=100, help='minimum target sites')

args = parser.parse_args()
minimum_hits = args.minimum_hits
reproduced_motifs = {}

# motifs = MEMEMotif.load_dreme('500/GC1.dreme')
# seq = 'TGGGCC'
# for key in motifs:
#     print('{}\t{}\t{}'.format(seq, key.name, check_match(key.name, seq)))
# exit()
columns = []
rows = set()
results = {}
srcdir = args.i
db = 'motifs/jaspar_arabid.meme'

def tomtom(filename_input, filename_output, db, params={}):
    tomtom = params.get('tomtom', 'tomtom')
    dstdir = params.get('oc', 'out')
    dist = params.get('dist', 'pearson')
    cmd = tomtom, '-oc', dstdir, '-dist', dist, '-evalue', '-min-overlap', '5', '-thresh', '10', fileanme_intput, db
    import subprocess
    sys.stderr.write(' '.join(cmd) + '\n')
    subprocess.Popen(cmd).wait()

for tf in args.t:
    sys.stderr.write('#{}\n'.format(tf))
    src_dreme = os.path.join(srcdir, tf + '.dreme/dreme.txt')
    src_tomtom = os.path.join(srcdir, tf + '.dreme/tomtom/tomtom.txt')
    # src_homer = os.path.join(srcdir, tf + '.homer/homerMotifs.all.motifs')
    srcdir_homer = os.path.join(srcdir, tf + '.homer/homerResults/')#.all.motifs')
    src_homer_tomtom = os.path.join(srcdir, tf + '.homer/tomtom/tomtom.txt')
    homer_tomtom_input = os.path.join(srcdir, tf + '.homer/tomtom.input')
    tomtom = collections.OrderedDict()
    homer = collections.OrderedDict()
    dreme = collections.OrderedDict()
    for mot in MEMEMotif.load_dreme(src_dreme):
        dreme[mot.name] = mot

    #HOMER
    motif_stat = {}
    homer_motifs = HomerMotif.load_motifs(srcdir_homer)
    if os.path.exists(src_homer_tomtom) is False:
        if os.path.exists(homer_tomtom_input) is False or os.path.getsize(homer_tomtom_input) < 1000:
            with open(homer_tomtom_input, 'w') as fo:
                fo.write(HomerMotif.convert_to_meme(homer_motifs))
        sys.stderr.write('{} is not prepared tomtom, execute tomtom with {}\n'.format(tf, homer_tomtom_input))
        tomtom(homer_tomtom_input, src_homer_tomtom, db)
    else:
        with open(src_homer_tomtom) as fi:
            for line in fi:
                items = line.strip().split('\t')
                if line.startswith('#') or len(items) < 5: continue
                # print(items)
                seq = items[0]
                dbname = items[1]
                for m in homer_motifs:
                    if m.motif == seq or m.name == seq or m.name == dbname or m.motif == dbname:
                        m.set_knwon_name(dbname)
                        if dbname not in homer:homer[dbname] = []
                        homer[dbname].append(seq)
                        # print('{} {} //  {} {}'.format(seq, dbname, m.motif, m.name))
                        break
    columns.append(tf)
    #DREME-TOMTOM
    if os.path.exists(src_tomtom) is False or os.path.getsize(src_tomtom) < 1000:
        tomtom(src_dreme, src_tomtom, db)
    with open(src_tomtom) as fi:
        for line in fi:
            items = line.strip().split('\t')
            if line.startswith('#') or len(items) < 5: continue
            name = items[0]
            dbname = items[1]
            if dbname not in tomtom: tomtom[dbname] = []
            tomtom[dbname].append(name)
            # tomtom[name].append(dbname)

    # print(len(tomtom), len(homer))
    results[tf] = {}
    for i, motif in enumerate(homer_motifs):
        if motif.hits < minimum_hits:
            continue
        # search database name
        touched = set()
        reproduced = False
        # print(motif.motif)
        for dbname, seqs in homer.items():
            for seq in seqs:
                # print(seq, motif.motif)
                if motif.motif == seq:
                    # print("NAMED")#motif.motif, seq)
                    if dbname in tomtom: # reproduced
                        reproduced = True
                        break
        if reproduced:
            name = motif.known
            rows.add(name)
            print('{}\t{}\t{:.3f}'.format(i + 1, name, motif.enrichment))#str(motif.dbname), motif.name))
            results[tf][motif.known] = motif

# output results
indexes = list(sorted(rows))
header = 'NAME\tMotif\t' + '\t'.join(args.t)
filename_output = args.o
cluster_script = 'python/cluster_motifs.py'
with open(filename_output, 'w') as fo:
    fo.write(header + '\n')
    for name in rows:
        seq = None
        values = []
        # print(name)
        for column in columns:
            result = results[column]
            if name in result:
                detected = result[name]
                if seq is None: seq = detected.motif
                values.append(detected.enrichment)
            else:
                values.append(0)
        if seq is not None:
            ostr = '{}\t{}'.format(name, seq)
            for v in values: ostr += '\t{:.4f}'.format(v)
            fo.write(ostr + '\n')

if os.path.exists(cluster_script):
    import subprocess
    filename_cluster = filename_output + '.cluster'
    cmd = 'python', cluster_script, '-i', filename_output, '-o', filename_cluster
    sys.stderr.write(' '.join(cmd) + '\n')
    subprocess.Popen(cmd).wait()