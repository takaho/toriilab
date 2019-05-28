import os, sys, re, collections
import sqlite3, subprocess
import argparse, gzip

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
        cur.execute('SELECT chromosome, start, stop, mapped FROM ' + table)
        saved = {}
        n = 0
        for r in cur.fetchall():
            if r[0] not in saved:
                saved[r[0]] = {}
            saved[r[0]][r[1]] = r
            n += 1

        # reuse counted if saved
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
                                touched.add(items[0])
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

def __load_tracks(bamfile, verbose=False):
    dbdir = os.path.expanduser('~/Research/Torii/tracks/')
    filename = os.path.join(dbdir, os.path.basename(bamfile) + '.gz')
    if os.path.exists(dbdir) is False:
        os.makedirs(dbdir)
    data = {}
    if not os.path.exists(filename) or os.path.getsize(filename) < 10000:
        sys.stderr.write('loading from bam file\n')
        if bamfile.endswith('.sam'):
            cmd = ['cat', bamfile]
        else:
            cmd = ['samtools', 'view', bamfile]
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)
        prev = [''] * 12
        n = m = 0
        for line in proc.stdout:
            line = line.decode('ascii')
            if line.startswith('@'): continue
            n += 1
            items = line.split('\t', 10)
            chrom = items[2]
            if chrom == '*': continue
            if chrom not in data:
                sys.stderr.write('reading from bam : {}        \r'.format(chrom))
                data[chrom] = [0] * 1000
            pair = items[6]
            if pair != '=': continue
            start = int(items[3])
            dist = int(items[8])
            if prev[0] == chrom and prev[1] == start and prev[2] == dist and prev[3] == items[9]:
                continue
            if dist > 0:
                pos = start + (len(items[9]) + dist) // 2
                track = data[chrom]
                rev = False
                while pos >= len(track):
                    rev = True
                    track += [0] * len(track)
                track[pos] += 1
                if rev:
                    data[chrom] = track
                m += 1
                if m % 1000 == 0:
                    sys.stderr.write('{} / {} : {}:{}      \r'.format(m // 1000, n // 1000, items[0], start))
            prev = chrom, start, dist, items[9]
        proc.stdout.close()
        proc.wait()
        with gzip.open(filename, 'wb') as fo:
            for chrom, vals in data.items():
                start = stop = 0
                for i in range(len(vals)):
                    if vals[i] > 0:
                        start = i
                        break
                for i in range(len(vals)):
                    j = len(vals) - 1 - i
                    if vals[j] > 0:
                        stop = j + 1
                        break
                sys.stderr.write('{} : {} {} bp \n'.format(bamfile, chrom, stop))
                fo.write('{}\t{}\t{}'.format(chrom, start, stop).encode('ascii'))
                for v in vals[start:stop]:
                    fo.write(' {}'.format(v).encode('ascii'))
                fo.write(b'\n')
    else:
        with gzip.open(filename) as fi:
            for line in fi:
                i = 0
                line = line.decode('ascii')
                l = len(line)
                start = col = 0
                chrom = '?'
                vals = [0]
                colstart = 0
                while i < l:
                    c = line[i]
                    if c <= ' ':
                        if col == 0:
                            chrom = line[colstart:i]
                            if verbose: sys.stderr.write('reading {} '.format(chrom))
                        elif col == 1:
                            start = int(line[colstart:i])
                        elif col == 2:
                            stop = int(line[colstart:i])
                            vals = [0] * (stop + 1)
                            data[chrom] = vals
                        else:
                            vals[col + start - 2] = int(line[colstart:i])
                            # if col > 100000: break
                        colstart = i + 1
                        col += 1
                    i += 1
                if verbose: sys.stderr.write(' {} reads              \n'.format(sum(vals)))
    return data

def count_reads_2(bamfile, regions, verbose=True):
    data = __load_tracks(bamfile, verbose)
    counts = [0] * len(regions)
    index = 0
    for region in regions:
        chrom, start, stop = region[0:3]
        vals = data.get(chrom, ())
        counts[index] = sum(vals[max(start,0):min(len(vals),stop + 1)])
        index += 1
    return counts
    #import hashlib
    # md5 = hashlib.md5()
    # os.path.abspath(bamfile)
    # filename_data = 

parser = argparse.ArgumentParser()
parser.add_argument('-i', nargs='+', default=None)
parser.add_argument('--merge', action='store_true')
parser.add_argument('-d', type=int, default=100, metavar='number', help='tolerance')
parser.add_argument('-o', default=None, metavar='filename', help='output filename')
parser.add_argument('--insert-others', action='store_true')
parser.add_argument('--verbose', action='store_true')
parser.add_argument('--upstream', metavar='number,number', default="500")
parser.add_argument('--downstream', metavar='number,number', default="0")
parser.add_argument('--gtf', default=None)
parser.add_argument('--exclude', default=None, metavar='bed file', help='removing regions such as centromere ')
#default='/mnt/smb/tae/tair10/Arabidopsis_thaliana.TAIR10.40.chrm.gtf')
args = parser.parse_args()
filename_output = args.o
tolerance = args.d
verbose = args.verbose
outputs_total = args.insert_others
merge_request = args.merge#False
peak_sets = collections.OrderedDict()
if args.gtf is not None:
    upstreams = [int(x_) for x_ in args.upstream.split(',')]
    downstreams = [int(x_) for x_ in args.downstream.split(',')]
    import promoterbed
    genes = promoterbed.load_gene_location(args.gtf)
    for upstream in upstreams:
        for downstream in downstreams:
            peaks = []
            for gene_id, info in genes.items():
                name, chrom, start, stop, ori = info
                if ori == '+':
                    p = start - upstream
                    q = start + downstream
                else:
                    p = stop - downstream
                    q = stop + upstream
                if gene_id.startswith('AT2G21280'):
                    name = 'SULA'
                if re.match('AT\\dG\\d+', name):
                    # name = '{};{}:{}-{}'
                    peak_name = gene_id
                else:
                    peak_name = '{}({})'.format(gene_id, name)
                    # name = '{}({});{}:{}-{}'.format(gene_id, name, chrom, p, q)
                peaks.append((chrom, p, q, peak_name))
            peak_sets['u{}d{}'.format(upstream, downstream)] = peaks
else:
    raise Exception('error')
    names = []
    for b in args.i:
        peaks = []
        if b.endswith('_summits.bed'):
            for p in load_summits(b):
                peaks.append((p[0], p[1] - tolerance, p[1] + tolerance))
        elif b.endswith('.bed') or b.endswith('.xls'):
            for p in load_bed(b):
                if p[0][3:].isdigit():
                    peaks.append(p)
        peak_sets[os.path.basename(b)] = peaks
    if merge_request:
        merged = []
        for p in peak_sts.values():
            merged += p
        # peaks = merge_peaks(peaks)
        peak_sets = dict(merged=merge_peaks(merged))
        # peak_sets['+'.join(args.i)] = peaks
#    peak_sets.append(peaks)

if args.exclude is not None:
    if verbose:
        sys.stderr.write('applying exclusion list\n')
    exclusion = []
    with open(args.exclude) as fi:
        for line in fi:
            items = re.split('\\s+', line)#fi.readline())#.split('\t')
            # print(items)
            if len(items) >= 3 and items[1].isdigit() and items[2].isdigit():
                chrom = re.sub('^C', 'c', items[0])
                start = int(items[1])
                stop = int(items[2])
                exclusion.append((chrom, start, stop))
                # print(chrom, start, stop)
    for identifier, peaks in peak_sets.items():
        available = []
        for chrom, p, q, peak_name in peaks:
            # if not re.match('chr\\d+$', chrom): 
            excluded = False
            for c, s, t in exclusion:
                if chrom == c and p <= t and q >= s:
                    excluded = True
                    # print(c, s, t)
                    break
            # if chrom[3:].isdigit() is False:
            #     print(excluded, exclusion)
            # if peak_name.find('SPCH') >= 0:
            #     print(peak_name, chrom, p, q, excluded)
            if not excluded: available.append((chrom, p, q, peak_name))
        if verbose:
            sys.stderr.write('reduce peaks : {} => {} using {}\n'.format(len(peaks), len(available), args.exclude))
    # exit()
        peak_sets[identifier] = sorted(available, key=lambda x_:int(x_[0][3:]) * 1000000000 | x_[1])


# print(len(peaks))
# print(len(peaks))

data = collections.OrderedDict()
bamdir = 'bamfiles'
countdata = {} # dict (key:filename) of dict (key:chrom, val:counts)
samples = []
celltypes = {}
totalreads = {}
with open('description.txt') as fi:
    for line in fi:
        if line.startswith('#'): continue
        items = line.strip().split('\t')
        if len(items) > 1:
            state = items[0]
            names = [x_.strip() for x_ in items[1].split(',')]
            bamfiles = []
            for bam in names:
                fn = os.path.join(bamdir, bam + '.bam')
                if os.path.exists(fn) and os.path.getsize(fn) > 100000:
                    bamfiles.append(fn)
                    sys.stderr.write('loading counts from {}:{}\n'.format(state, bam))
                    countdata[fn] = __load_tracks(fn, verbose)# count_reads_2(fn, peaks)
                    if outputs_total: 
                        totalreads[fn] = get_total_count(fn)
                    # break
            if len(bamfiles) > 0:
                    # celltypes[state] = celltypes.get(state, 0) + 1
                if state in data:
                    data[state].append(bamfiles)
                else:
                    data[state] = [bamfiles, ]
                samples.append('{}_rep{}'.format(state, len(data[state])))
                # break

print(data.keys())
for state in data.keys():
    bamfileset = data[state]
    for bams in bamfileset:
        for bam in bams:
            print(bam)
            print(bam in countdata)
        # print(countdata[fn].__class__)
        # print(countdata[fn])
    # print(val)
# print(samples)
# exit()

ostr = sys.stdout
filename_output = args.o
outdir = args.o
for identifier, peaks in peak_sets.items():
# for peak_index, peaks in enumerate(peak_sets):
    if outdir is not None:
        if os.path.isdir(outdir):
            filename = os.path.join(outdir, identifier + '.tsv')
            # ostr = open(filename_output, 'w')
        else:
            filename_output = outdir
            bn = os.path.basename(filename_output)
            if bn.count('.') > 0:
                filename = os.path.join(os.path.dirname(filename_output), bn[0:bn.rfind('.')], '.{}'.format(peak_index), '.tsv')
            else:
                filename = os.path.join(os.path.dirname(filename_output), bn, '.{}'.format(peak_index), '.tsv')
        ostr = open(filename, 'w')
        if verbose:
            sys.stderr.write('{} : write to {}\n'.format(identifier, filename))

    ostr.write('{}\t'.format(re.sub('\\W', '_', identifier)) + '\t'.join(samples) + '\n')
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
                for bamfile in filenames:
                    try:
                        counts += sum(countdata[bamfile][chrom][start:stop])
                    except Exception as e:
                        sys.stderr.write('{} {} \r'.format(bamfile, bamfile in countdata))
                        raise e
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
# exit()

# # for chromosome, position, enr in peaks:

# for bam in os.listdir('bamfiles'):
#     if bam.endswith('.bam'):
#         fn = os.path.join('bamfiles', bam)
#         if os.path.getsize(fn) > 100000:
#             counts = count_reads(fn, peaks)

# p = load_summits('peakcall/AtML1_rep1_summits.bed')
# peaks = []
# tolerance = 250
# for p in load_summits('peakcall/AtML1_rep1_summits.bed'):
#     peaks.append((p[0], p[1] - tolerance, p[1] + tolerance))

# counts = count_reads('bamfiles/C2_S29.bam', peaks[0:250])
# exit()

# import collections
# peaks = collections.OrderedDict()
# combined = []
# for bed in sys.argv[1:]:
#     name = os.path.basename(bed)
#     peaks[name] = load_summits()
#     combined += peaks[name]
