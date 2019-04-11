"""
Aggregate peaks detected using HOMER and MEME (DREME)
"""
import os, sys, re
import pathlib, argparse
import numpy as np
# import pandas as pd
import collections
from motifdetection import MEMEMotif, HOMERMotif, ATACMotif, load_consensus

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
parser.add_argument('-n', type=int, default=10)
parser.add_argument('--motif', default=None)
parser.add_argument('--minimum-hits', type=int, default=100, help='minimum target sites')

args = parser.parse_args()
filename_output = args.o
minimum_hits = args.minimum_hits
reproduced_motifs = {}
hit_threshold = args.n

if args.motif is not None:
    motif_db = load_consensus(args.motif)
else:
    motif_db = {}
for key, val in motif_db.items():
    print('{}\t{}'.format(key, val))
exit()    
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

celltypes = args.t

for celltype in celltypes:
    sys.stderr.write('#{}\n'.format(celltype))
    src_dreme = os.path.join(srcdir, celltype + '.dreme/dreme.txt')
    src_tomtom = os.path.join(srcdir, celltype + '.dreme/tomtom/tomtom.txt')
    # src_homer = os.path.join(srcdir, celltype + '.homer/homerMotifs.all.motifs')
    srcdir_homer = os.path.join(srcdir, celltype + '.homer/homerResults/')#.all.motifs')
    src_homer_tomtom = os.path.join(srcdir, celltype + '.homer/tomtom/tomtom.txt')
    homer_tomtom_input = os.path.join(srcdir, celltype + '.homer/tomtom.input')
    tomtom = collections.OrderedDict()
    homer = collections.OrderedDict()
    dreme = collections.OrderedDict()

    # DREME
    dreme_motifs = MEMEMotif.load_motifs(src_dreme)
    ATACMotif.annotate_motifs(dreme_motifs, src_tomtom)
    for mot in dreme_motifs:
        dreme[mot.name] = mot

    #HOMER
    motif_stat = {}
    homer_motifs = HOMERMotif.load_motifs(srcdir_homer)
    if os.path.exists(src_homer_tomtom) is False: # execute if no tomtom outputs found
        if os.path.exists(homer_tomtom_input) is False or os.path.getsize(homer_tomtom_input) < 1000:
            with open(homer_tomtom_input, 'w') as fo:
                fo.write(HOMERMotif.convert_to_meme(homer_motifs))
        sys.stderr.write('{} is not prepared tomtom, execute tomtom with {}\n'.format(celltype, homer_tomtom_input))
        tomtom(homer_tomtom_input, src_homer_tomtom, db)
    ATACMotif.annotate_motifs(homer_motifs, src_homer_tomtom)#tomtom)

    # Examine reproducibility
    shared_names = set()
    reproduced_motifs = []
    print('MEME {}, HOMER {}'.format(len(dreme_motifs), len(homer_motifs)))
    # for m in homer_motifs:
    #     print(m.get_known())
    # print()
    # for m in dreme_motifs:
    #     print(m.get_known())
    # exit()
    for motif in homer_motifs:
        reproduced = False
        for m in dreme_motifs:
            # print(motif.consensus, m.consensus, m.get_known())
            if motif.consensus in m.consensus or m.consensus in motif.consensus:
                reproduced = True
            for c in motif.get_shared_consensus(m):
                shared_names.add(c)
                reproduced = True
        if reproduced:
            reproduced_motifs.append(motif)
    arranged = ATACMotif.classify_and_select_best(reproduced_motifs)
    # print(len(arranged))
    results[celltype] = arranged
    for name, m in arranged.items():
        print('{}\t{}\t{:.2f}\t{:.0f}\t{}'.format(celltype, m.consensus, m.enrichment, m.n_sites, name))

normalized_names = {}
enrichment_matrix = {}
all_motifs = set()
for celltype, motifs in results.items():
    results = {}
    used_motif = set()
    for joined_name, motif in motifs.items():
        names = joined_name.split('//')
        joined = None
        for name in names:
            if name in results:
                if results[name].n_sites < motif.n_sites:
                    results[name] = motif
            else:
                results[name] = motif
                all_motifs.add(name)
    enrichment_matrix[celltype] = results

def name_with_priority(name):
    if name in ('SEP1', 'SOC1', 'PI', 'AGL15', 'ATHB-5', 'SEP3', 'TBP3', 'SEP4', 'SEP5',  'TCP15'):
        return ' ' + name
    return name

used_ids = set()
unique_names = []
for m in sorted(all_motifs, key=name_with_priority):
    ids = []
    novel = False
    for enr in enrichment_matrix.values():
        if m in enr:
            idnum = id(enr[m])
            if idnum not in used_ids:
                novel = True
                used_ids.add(idnum)
            ids.append(idnum)
    if novel:
        unique_names.append(m)

#valuematrix = np.zeros((len(unique_names), len(celltypes)))
#valuematrix = []
with open(filename_output, 'w') as fo:
    header = 'NAME\tMotif\t' + '\t'.join(celltypes) + '\n'
    fo.write(header)
    for name in unique_names:
        row = '{}\t{}'.format(name, motif_db.get(name, name))
        for celltype in celltypes:
            enr = enrichment_matrix[celltype]
            if name in enr and enr[name].n_sites > hit_threshold:
                value = enr[name].enrichment
                visible = True
            else:
                value = 1.0
            row += '\t{:.6f}'.format(value)
        if visible:
            fo.write(row + '\n')#\t{:.6f}'.format(value))

#df.to_csv(filename_output, sep='\t')
exit()
#             if name in normalized_names:
#                 joined = normalized_names[name]
#                 for n in names:
#                     joined.add(n)
#                 for n in names:
#                     normalized_names[n] = joined
#                 break
#         if joined is None:
#             joined = set(names)
#             for n in names:
#                 normalized_names[n] = joined
# touched = set()
# for n, alias in normalized_names.items():
#     if id(alias) in touched:
#         continue
#     print(','.join(sorted(alias)))
#     touched.add(id(alias))
# print(len(touched))
# exit()
                




# exit()

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
                values.append(1)
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