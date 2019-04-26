import os, sys, re, collections
import sqlite3, subprocess
import argparse

def setup_database(bamfile, table=None):
    if table is None:
        table = re.sub('\\W', '', os.path.basename(bamfile).split('.')[0])
    dbdir = os.path.expanduser('~/Research/Torii/db/')
    db = os.path.join(dbdir, 'ataccount.db')
    with sqlite3.connect(db) as cnx:
        cur = cnx.cursor()
        cur.execute('select name from sqlite_master where name="{}"'.format(table))#counts"')
        if cur.fetchone() is None:
            cur.execute('CREATE TABLE {} (chromosome not null, start int4 not null, stop int4 not null, mapped int4 not null)'.format(table))
            cur.execute('CREATE INDEX index_{0} ON {0}(chromosome, start, stop)'.format(table))
        cur.execute('SELECT NAME FROM sqlite_master WHERE NAME="total"')
        r = cur.fetchone()
        if r is None:
            cur.execute('CREATE TABLE total (name NOT NULL, count INT4 NOT NULL)')
    return db

def get_total_count(bamfile):
    table = re.sub('\\W', '', os.path.basename(bamfile).split('.')[0])
    db = setup_database(bamfile, table)
    with sqlite3.connect(db) as cnx:
        cur = cnx.cursor()
        cur.execute('SELECT count FROM total WHERE name=?', (table, ))
        r = cur.fetchone()
        num_mapped = 0
        if r is None or r[0] == 0:
            cmd = 'samtools', 'flagstat', bamfile
            sys.stderr.write(' '.join(cmd) + '\n')
            proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)
            for line in proc.stdout:
                sys.stderr.write(line.decode('utf8'))
                m = re.match('(\\d+)\\s+\\+\\s+(\\d+)\\s+mapped', line.decode('utf8'))
                if m:
                    num_mapped = int(m.group(1))
            proc.stdout.close()
            proc.wait()
            if num_mapped > 0:
                cur.execute('INSERT INTO total values(?, ?)', (table, num_mapped))
                cnx.commit()
        else:
            num_mapped = r[0]
        return num_mapped

def count_reads(bamfile, regions):
    bai = bamfile + '.bai'
    counts = [0] * len(regions)
    if os.path.exists(bai) is False:
        cmd = 'samtools', 'index', bam
        subprocess.Popen(cmd).wait()
    table = re.sub('\\W', '', os.path.basename(bamfile).split('.')[0])
    # print(table)
    db = setup_database(bamfile, table)
    # print(db)
    with sqlite3.connect(db) as cnx:
        cur = cnx.cursor()
        # cur.execute('SELECT count FROM total WHERE name=?', (table, ))
        # r = cur.fetchone()
        # if r is None:
        #     cmd = 'samtools', 'flagstat', bamfile
        #     sys.stderr.write(' '.join(cmd) + '\n')
        #     proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)
        #     num_mapped = 0
        #     for line in proc.stdout:
        #         sys.stderr.write(line.decode('utf8'))
        #         m = re.match('(\\d+)\\s+\\+\\s+(\\d+)\\s+mapped', line.decode('utf8'))
        #         if m:
        #             num_mapped = int(m.group(1))
        #     proc.stdout.close()
        #     proc.wait()
        #     if num_mapped > 0:
        #         cur.execute('INSERT INTO total values(?, ?)', (table, num_mapped))
        cur.execute('SELECT chromosome, start, stop, mapped FROM ' + table)
        saved = {}
        n = 0
        for r in cur.fetchall():
            if r[0] not in saved:
                saved[r[0]] = {}
            saved[r[0]][r[1]] = r
            n += 1
        # sys.stderr.write('{} counts saved  \n'.format(n))
        remnants = []
        for i, r in enumerate(regions):
            if len(r) >= 3:
                if r[0] in saved and r[1] in saved[r[0]]:
                    counts[i] = saved[r[0]][r[1]][3]
                    # print(saved[r[0]][r[1]])
                else:
                    remnants.append(i)
        chunk_size = 200
        i = 0
        while i < len(remnants):
            index = remnants[i]
            cmd = ['samtools', 'view', bamfile]
            stop = min(i + chunk_size, len(remnants))
            # print(index, stop, stop - index)
            for j in range(i, stop):
                c_, s_, st_ = regions[remnants[j]][0:3]
                cmd.append('{}:{}-{}'.format(c_, s_, st_))
            # sys.stderr.write(' '.join(cmd) + '\n')
            proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)
            sys.stderr.write('{:.2f}% {}           \r'.format(100.0 * index / len(regions), bamfile))
            # sys.stderr.write(' '.join(cmd) + '\n')
            touched = set()
            for line in proc.stdout:
                items = line.decode('utf-8').strip().split('\t')
                if items[0] in touched: continue
                touched.add(items[0])
                chrom = items[2]
                # print(items)
                if items[6] == '=':
                    st = int(items[3])
                    dist = int(items[8])
                    if abs(dist) < 10000:
                        length = len(items[9])
                        center = st + (dist + length) // 2
                        for j in range(i, stop):
                            c_, s_, e_ = regions[remnants[j]][0:3]
                            if items[2] == c_ and s_ <= center <= e_:
                                counts[remnants[j]] += 1
            i += chunk_size
            proc.stdout.close()
            proc.wait()
            del proc

        for index in remnants:
            chrom, start, stop = regions[index][0:3]
            # print(i, counts[i])
            cur.execute('INSERT INTO ' + table + ' VALUES(?, ?, ?, ?)', (chrom, start, stop, counts[index]))
        cnx.commit()
    return counts

# load peask
def load_summits(filename):
    summits = []
    with open(filename) as fi:
        for line in fi:
            items = line.strip().split('\t')
            chrom = items[0]
            if chrom[3:].isdigit():
                start = int(items[1])
                if len(items) > 4:
                    enrichment = float(items[4])
                else:
                    enrichment = 1
                summits.append((chrom, start, enrichment))
    return summits

def merge_peaks(peaks, distance=50):
    def get_locationkey(pos):
        c = int(pos[0][3:]) * 1000000000 + (pos[1] + pos[2]) // 2
        return c
    peaks = sorted([p for p in peaks if p[0][3:].isdigit()], key=get_locationkey)
    merged = []
    i = 0
    while i < len(peaks):
        p = peaks[i]
        # print(i, p)
        chrom, start, stop = p[0:3]
        p0 = p[1]
        j = i + 1
        while j < len(peaks):
            q = peaks[j]
            p1 = q[1]
            if p[0] != q[0] or p0 + distance < p1:
                break
            start += p1
            stop += q[2]
            j += 1
        n = j - i
        if n > 1:
            # print(peaks[i:j])
            merged.append((chrom, start // n, stop / n))
        else:
            merged.append(peaks[i])
        i = j
    return merged

def load_bed(filename):
    peaks = []
    with open(filename) as fi:
        for line in fi:
            items = line.strip().split('\t')
            if len(items) >= 3 and items[1].isdigit() and items[2].isdigit():
                chrom = items[0]
                start = int(items[1])
                stop = int(items[2])
                name = items[3]
                peaks.append((chrom, start, stop, name))
    return peaks


parser = argparse.ArgumentParser()
parser.add_argument('-i', nargs='+', default=None)
parser.add_argument('-d', type=int, default=100, metavar='number', help='tolerance')
parser.add_argument('-o', default=None, metavar='filename', help='output filename')
parser.add_argument('--insert-others', action='store_true')
args = parser.parse_args()
filename_output = args.o
peaks = []
tolerance = args.d
outputs_total = args.insert_others
merge_request = False
for b in args.i:
    if len(peaks) > 0:
        merge_request = True
    if b.endswith('_summits.bed'):
        for p in load_summits(b):
            peaks.append((p[0], p[1] - tolerance, p[1] + tolerance))
    elif b.endswith('.bed') or b.endswith('.xls'):
        for p in load_bed(b):
            if p[0][3:].isdigit():
                peaks.append(p)
if merge_request:
    peaks = merge_peaks(peaks)
# print(len(peaks))
# print(len(peaks))

data = collections.OrderedDict()
bamdir = 'bamfiles'
countdata = {}
samples = []
celltypes = {}
totalreads = {}
with open('description.txt') as fi:
    for line in fi:
        items = line.strip().split('\t')
        if len(items) > 1:
            state = items[0]
            names = [x_.strip() for x_ in items[1].split(',')]
            bamfiles = []
            for bam in names:
                fn = os.path.join(bamdir, bam + '.bam')
                print(fn)
                if os.path.exists(fn) and os.path.getsize(fn) > 100000:
                    bamfiles.append(fn)
                    sys.stderr.write('loading counts from {}:{}\n'.format(state, bam))
                    countdata[fn] = count_reads(fn, peaks)
                    if outputs_total: 
                        totalreads[fn] = get_total_count(fn)
            if len(bamfiles) > 0:
                    # celltypes[state] = celltypes.get(state, 0) + 1
                if state in data:
                    data[state].append(bamfiles)
                else:
                    data[state] = [bamfiles, ]
                samples.append('{}_rep{}'.format(state, len(data[state])))
# print(samples)

ostr = sys.stdout
if args.o is not None:
    ostr = open(args.o, 'w')

ostr.write('Peak\t' + '\t'.join(samples) + '\n')
index = 0
assigned = {}
for key in data.keys(): assigned[key] = 0
for i in range(len(peaks)):
    peak = peaks[i]
    if len(peak) > 3:
        chrom, start, stop, name = peaks[i][0:4]
        name = name + ';'
    else:
        chrom, start, stop = peaks[i]
        name = ''
    ost = '{}{}:{}-{}'.format(name, chrom, start, stop)
    row = []
    for state in data.keys(): # AtML1 SPCH MUTE FAMA GC1
        assigned[state] = []
        for filenames in data[state]:
            counts = 0
            for filename in filenames:
                counts += countdata[filename][i]
            row.append(counts)
            ost += '\t{}'.format(counts)
            assigned[state].append(counts)
    ostr.write(ost + '\n')
if outputs_total:
    ost = 'others'
    for state in data.keys():
        for i, filenames in enumerate(data[state]):
            total = 0
            for filename in filenames:
                total += totalreads[filename]
            ost += '\t{}'.format(total - assigned[state][i])
    ostr.write(ost + '\n')

if args.o is not None: ostr.close()

#print(data)
exit()

# for chromosome, position, enr in peaks:

for bam in os.listdir('bamfiles'):
    if bam.endswith('.bam'):
        fn = os.path.join('bamfiles', bam)
        if os.path.getsize(fn) > 100000:
            counts = count_reads(fn, peaks)

p = load_summits('peakcall/AtML1_rep1_summits.bed')
peaks = []
tolerance = 250
for p in load_summits('peakcall/AtML1_rep1_summits.bed'):
    peaks.append((p[0], p[1] - tolerance, p[1] + tolerance))

counts = count_reads('bamfiles/C2_S29.bam', peaks[0:250])
exit()

import collections
peaks = collections.OrderedDict()
combined = []
for bed in sys.argv[1:]:
    name = os.path.basename(bed)
    peaks[name] = load_summits()
    combined += peaks[name]
