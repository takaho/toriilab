import os, sys, re, argparse
import argparse
import pandas as pd
import numpy as np
import sklearn.cluster

# def quantile_normalization_(matrix):
#     ranks = matrix.stack().groupby(matrix.rank(method='first').stack().astype(int)).mean()
#     return matrix.rank(method='min').stack().astype(int).map(ranks).unstack()

def quantile_normalization(df):

    # columns from the original dataframe not specified in cols
    non_numeric = df.filter(items=list(filter(lambda col: col not in df.columns, list(df))))

    rank_mean = df.stack().groupby(df.rank(method='first').stack().astype(int)).mean()  

    norm = df.rank(method='min').stack().astype(int).map(rank_mean).unstack()


    result = pd.concat([norm, non_numeric], axis=1)
    return result


def load_normalized_expression(filename, use_zscore=False, aggregation=True, normalization='quantile', rna=None):
    import gzip, pickle, hashlib
    md5 = hashlib.md5()
    md5.update('{}_{}_{}'.format(os.path.basename(filename), int(use_zscore), int(aggregation)).encode('utf-8'))
    if rna is not None:
        md5.update(rna.encode('utf-8'))
    if normalization is not None:
        md5.update(('__' + normalization).encode('utf-8'))
    cache = os.path.join(os.path.dirname(filename), '.' + md5.hexdigest()[0:12] + '.norm.cache')
    # print(cache)
    if os.path.exists(cache):
        with gzip.open(cache) as fi:
            df = pickle.load(fi)
        return df
    if rna is not None:
        name = None
        size = 0
        lengths = {}
        with open(rna) as fi:
            for line in fi:
                if line.startswith('>'):
                    if name is not None and size > 0:
                        lengths[name] = size
                    name = line[1:].strip()
                    size = 0
                elif name is not None:
                    size += len(line.strip())
        if name is not None and size > 0:
            lengths[name] = size
    else:
        lengths = None

    exprdata = pd.read_csv(filename, index_col=0, sep='\t')
    if filename.endswith('.count'):
        tpm = exprdata.copy()
        for col in exprdata.columns:
            sys.stderr.write('Total reads : {}\t{}\n'.format(col, np.sum(exprdata[col])))
            tpm[col] *= 1e6 / np.sum(exprdata[col])
        expradta = tpm
    if filename.endswith('.rpkm') is False and lengths is not None:
        for r in range(exprdata.shape[0]):
            exprdata.iloc[r] *= 1e3 / lengths[exprdata.index[r]]

    genes = {}
    for gene in exprdata.index:
        agi = gene.split('.')[0]
        if agi not in genes:
            genes[agi] = [gene, ]
        else:
            genes[agi].append(gene)
    agigenes = sorted(genes.keys())
    values = []
    for agi in agigenes:
        alt = genes[agi]
        row = exprdata.loc[alt[0]].copy()
        for a in alt[1:]:
            row += exprdata.loc[a]
        values.append(row)
    df = pd.DataFrame(values, columns=exprdata.columns, index=agigenes)

    samples = []
    sample2column = {}
    for i, col in enumerate(df.columns):
        name = col.split('_', 1)[0]
        if name not in samples:
            samples.append(name)
            sample2column[name] = []
        sample2column[name].append(i)

    # data = {}
    # with open(filename) as fi:
    #     columns = fi.readline().strip().split('\t')[1:]
    #     samples = []
    #     sample2column = {}
    #     for i, col in enumerate(columns):
    #         name = col.split('_', 1)[0]
    #         if len(samples) == 0 or samples[-1] != name:
    #             samples.append(name)
    #             sample2column[name] = []
    #         sample2column[name].append(i)
    #     for line in fi:
    #         items = line.strip().split('\t')
    #         gene = items[0].split('.')[0]
    #         if gene not in data: 
    #             data[gene] = np.array([float(x_) for x_ in items[1:]])
    #         else:
    #             data[gene] += np.array([float(x_) for x_ in items[1:]])
    # genes = sorted(data.keys())
    # df = pd.DataFrame([data[g] for g in genes], index=genes, columns=columns)

    # df[df<0.01] = 0.01
    # df = np.log2(df)
    normalized = df
    for col in normalized.columns:
        values = normalized[col]
        n = len([x_ for x_ in values.values if 0.1 < x_ < 1.0])
        print('{}\t{}'.format(col, n))

    if normalization == 'quantile':
        sys.stderr.write('quantile normalization\n')
        normalized = quantile_normalization(df)
    else:
        normalized = df

    for col in normalized.columns:
        values = normalized[col]
        n = len([x_ for x_ in values.values if 0.1 < x_ < 1.0])
        print('{}\t{}'.format(col, n))


    if 0:
        if use_zscore:
            mat = normalized.copy()
            nrow, ncol = mat.shape
            eavg = np.mean(mat, axis=1).values.reshape(nrow, 1)
            sd = np.std(mat, axis=1).values.reshape(nrow, 1)
            mat = np.divide(np.subtract(mat, eavg), sd)#np.std(mat, axis=1).values.reshape(nrow, 1))
            normalized = mat
        else:
            mat[mat<0.01] = 0.01
            mat = np.log(mat)
            normalized = mat
    # aggregation

    if aggregation:
        aggregated = []
        colsets = []
        for sample in samples:
            cols = []
            for col in sample2column[sample]:
                cols.append(col)
            colsets.append(cols)
        for gene in genes:
            row = normalized.loc[gene].values
            vals = []
            for cols in colsets:
                vals.append(np.mean([row[col] for col in cols]))
            aggregated.append(vals)
        expression = pd.DataFrame(aggregated, columns=samples, index=genes)
    else:
        expression = normalized
    with gzip.open(cache, 'wb') as fo:
        pickle.dump(expression, fo)
    return expression

def load_expression(filename):
    data = {}
    with open(filename) as fi:
        columns = fi.readline().strip().split('\t')[1:]
        samples = []
        sample2column = {}
        for i, col in enumerate(columns):
            name = col.split('_', 1)[0]
            if len(samples) == 0 or samples[-1] != name:
                samples.append(name)
                sample2column[name] = []
            sample2column[name].append(i + 1)
        for line in fi:
            items = line.strip().split('\t')
            gene = items[0].split('.')[0]
            if gene not in data: data[gene] = [0] * len(samples)
            row = data[gene]
            for i, name in enumerate(samples):
                cols = sample2column[name]
                val = 0
                for col in cols:
                    val += float(items[col])
                row[i] += val / len(cols)
    return data

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', default='THS/500bp.tsv')
    #parser.add_argument('-o', default=None, help='clustered results')
    parser.add_argument('-e', default='SRP043607.tpm', help='expression filename')
    parser.add_argument('-n', type=int, default=10, metavar='number', help='number of clusters')
    # parser.add_argument('--fdr', default=0.05, type=float)
    # parser.add_argument('--path', default='ebseq.out')
    parser.add_argument('--quantile', action='store_true')
    parser.add_argument('--by-sample', action='store_true')
    parser.add_argument('--index', default='raw', metavar='name', help='zscore/logfc/raw')
    parser.add_argument('-o', default='atac_out')#utdir', default='out')
    parser.add_argument('--verbose', action='store_true')
    #parser.add_argument('--verbose', action='store_true')

    # parser.add_argument('--relative', action='store_true')
    # parser.add_argument('--zscore', action='store_true')
    args = parser.parse_args()

    index_mode = args.index
    outdir = args.o#utdir

    verbose = args.verbose
    if os.path.exists(outdir) is False:
        os.makedirs(outdir)
    filename_input = args.i#'THS/500bp.tsv'
    fdr = 0.05
    n_clusters = args.n
    displayall_mode = args.by_sample
    filename_graph = os.path.join(outdir, 'heatmap.html')
    filename_output = os.path.join(outdir, 'summary.tsv')
    filename_dump = os.path.join(outdir, 'results.tsv')
    filename_expression = args.e

    normalization = 'quantile' if args.quantile else 'none'
    if displayall_mode:
        display_mode = 'all data'
    else:
        display_mode = 'mean of conditoin'

    pthr = 0.05
    preferences = """Input : {input}
Output : {output}
Graph : {graph}
Clusters : {clusters}
Normalization : {normalization}
Index : {index}
Values : {values}
Expression : {expression}
""".format(input=filename_input, graph=filename_graph, output=filename_output, clusters=n_clusters, normalization=normalization, index=index_mode, values=display_mode, expression=filename_expression)
    with open(os.path.join(outdir, 'run_info.txt'), 'w') as fo:
        fo.write(preferences)
    if verbose:
        sys.stderr.write(preferences)

    paths = {}

    conditions = []
    levels = []
    data = {}
    nums = []

    with open(filename_input) as fi:
        header = fi.readline().split('\t')
        for item in header:
            # if item in skipping: continue
            elems = item.split('_', 1)
            # if item == 'MUTE_rep1' or item == 'MUTE_rep2': 
            #     print('skip', item)
            #     continue
            if len(elems) >= 2:
                cnd = elems[0]
                conditions.append(cnd)
                if len(levels) == 0 or levels[-1] != cnd:
                    levels.append(cnd)
                    nums.append(1)
                else:
                    nums[-1] += 1
        num_points = len(nums)
    atac = pd.read_csv(filename_input, sep='\t', index_col=0)
    columns = atac.columns#[c for c in atac if c not in skipping]
    # print(columns)
    atac = atac[columns]
    genic = [x_ for x_ in atac.index if x_.find('rRNA') < 0]
    if args.quantile:
        matrix = quantile_normalization(atac.loc[genic])
    else:
        matrix = atac.loc[genic]

    # hmat = matrix.copy()
    # # hmat[hmat < 0.01] = 0.01
    # # hmat = np.log(hmat)
    # clu = sklearn.cluster.KMeans(5)
    # clu.fit_predict(hmat.T)
    # for i, l in enumerate(clu.labels_):
    #     print('{}\t{}'.format(clu.labels_[i], matrix.columns[i]))
    # exit()

    values = matrix.values
    import scipy.stats
    # import warnings

    significant = []
    valmat = []
    pvalues = []
    pthr = 0.01
    for i, elem in enumerate(matrix.index):
        row = values[i]
        col = 0
        r_ = []
        valset = []
        for j in range(num_points):
            r_.append(np.mean(row[col:col+nums[j]]))
            valset.append(row[col:col+nums[j]])
            col += nums[j]
        vmin = min(r_)
        vmax = max(r_)
        if vmin == vmax:
            pval = 1.0
        else:
            pval = scipy.stats.f_oneway(*valset).pvalue
            # with warnings.catch_warnings():
            #     warnings.filterwarnings('error')
            #     try:
            #         pval = scipy.stats.f_oneway(*valset).pvalue
            #     except RuntimeWarning as e:
            #         print(valset)
            #         print(elem)
            #         pval = 1.0
            #         exit()
        # if pval < 0.01:
        #     print(elem, vmin, vmax, pval, valset)
        g_ = matrix.index[i]
        # if g_.find('(GC1)') >= 0 or g_.find('ATML1') >= 0 or g_.find('SPCH') >= 0 or g_.find('MUTE') >= 0 or g_.find('FAMA') >= 0:
        #     print(g_.split(';')[0], vmin, vmax, pval, valset)
        pvalues.append(pval)
        if pval < pthr:#0.01:# and vmax > vmin * 2.0:
        # if pval < 0.01:# and vmax > vmin * 2.0:
            ostr = matrix.index[i]
            for v in r_:
                ostr += '\t{:.1f}'.format(v)
            # print(ostr)
            significant.append(elem)#matrix.index[i])
        else:
            g_ = matrix.index[i]
            if g_.find('(SPCH)') >= 0 or g_.find('(MUTE)') >= 0 or g_.find('(ATML1)') >= 0 or g_.find('(FAMA)') >= 0:
                significant.append(g_)            
        data[elem] = r_

    preferences = """Genes : {num_genes}
P-value : {pvalue}
Significant genes : {num_significant}
   """.format(num_genes = len(genic), num_significant=len(significant), pvalue=pthr) 
    with open(os.path.join(outdir, 'run_info.txt'), 'a') as fo:
        fo.write(preferences)
    if verbose:
        sys.stderr.write(preferences)

    # if verbose:
    #     sys.stderr.write('number of significant elements : {}\n'.format(len(significant)))

    # import statsmodels.stats.multitest
    # results = statsmodels.stats.multitest.multipletests(pvalues, alpha=0.1)
    # print(len(significant))
    # significant = []
    # for i, rejected in enumerate(results[0]):
    #     if rejected:
    #         significant.append(matrix.index[i])
    # print(len(significant))
    # exit()

    selected_matrix = pd.DataFrame([data[x_] for x_ in significant], index=significant, columns=levels)
    #selected_matrix = matrix.loc[significant]#pd.DataFrame([data[x_] for x_ in significant], index=significant, columns=levels)
    n_clusters = args.n
    clu = sklearn.cluster.AgglomerativeClustering(n_clusters)#, affinity='cosine', linkage='complete')
    #clu = sklearn.cluster.KMeans(n_clusters)#, affinity='cosine', linkage='complete')

    relative = []
    # selected_matrix = matrix
    # selected_matrix[selected_matrix<1] = 1
    # lm = np.log2(selected_matrix)

    if index_mode == 'logfc':
        heatmap_matrix = selected_matrix.copy()
        heatmap_matrix[heatmap_matrix<1] = 1
        heatmap_matrix = np.log2(heatmap_matrix)
        heatmap_matrix = np.subtract(heatmap_matrix, np.median(heatmap_matrix, axis=1).reshape(heatmap_matrix.shape[0], 1))
        vmin, vmax = -3, 3
    elif index_mode == 'zscore':
        nrows = selected_matrix.shape[0]
        mean = np.mean(selected_matrix, axis=1).values.reshape(nrows, 1)
        sd = np.sqrt(np.std(selected_matrix, axis=1)).values.reshape(nrows, 1)
        heatmap_matrix = np.divide(np.subtract(selected_matrix, mean), sd)
        if displayall_mode:
            vmin, vmax = -4, 4
        else:
            vmin, vmax = -6, 6
    elif index_mode == 'log':
        heatmap_matrix = np.log2(selected_matrix)
        vmin, vmax = 0, 8
    else:
        heatmap_matrix = selected_matrix
        vmin, vmax = 0, 1000
    clu.fit_predict(heatmap_matrix)

    #clone = matrix.copy()
    #clone[clone<1] = 1
    if displayall_mode:
        s_ = matrix.loc[significant].copy()
        n_ = s_.shape[0]
        if index_mode == 'logfc':
            s_[s_<1] = 1
            s_ = np.log2(s_)
            display_data = np.subtract(s_, np.median(s, axis=1).reshape(n_, 1))
        elif index_mode == 'zscore':
            display_data = pd.DataFrame(np.divide(np.subtract(s_, np.mean(s_, axis=1).values.reshape(n_, 1)), np.std(s_, axis=1).values.reshape(n_, 1)), index=significant, columns=s_.columns)
            vmin, vmax = -2, 2
        else:
            display_data = np.log2(s_)
        # display_data = np.log2(s_)
        # if args.relative:
        #     display_data = pd.DataFrame(np.subtract(display_data, np.median(display_data, axis=1).reshape,(display_data.shape[0], 1)), index=significant, columns=matrix.columns)
    else:
        display_data = heatmap_matrix
    # display_data = lm
    # display_data = matrix.copy()
    # display_data[display_data<1] = 1
    # display_data =np.log2(display_data)
    # display_data = np.subtract(display_data, np.median(display_data, axis=1).reshape(display_data.shape[0], 1))


    use_zscore_expr = True
    expression = load_normalized_expression(filename_expression, normalization='quantile', rna='/mnt/smb/tae/tair10/cdna.fa')
    #print(expression)
    #exit()

    # display_data = lm
    z = []
    x = display_data.columns
    y = []

    # for c in range(n_clusters):
    #     for i, s in enumerate(significant):
    #         if clu.labels_[i] == c:

    expressedgenes = set(expression.index)

    import scipy.stats
    svalues = heatmap_matrix.values
    clusters = [[] for i in range(n_clusters)]
    n_points = len(svalues[0])
    stack = np.zeros((n_clusters, n_points))
    # print(min(clu.labels_), max(clu.labels_))
    for i, l in enumerate(clu.labels_):
        clusters[l].append(i)
        stack[l] += np.array(svalues[i])
    # print([len(x_) for x_ in clusters])
    slope = np.zeros(n_clusters)
    for i in range(n_clusters):
        # print(stack[i])
        slope[i] = np.argmax(stack[i]) * 10 + scipy.stats.linregress(np.arange(n_points), stack[i] / max(1, len(clusters[i])))[0]
    # print(slope)
    ordered_clusters = sorted(np.arange(0, n_clusters), key=lambda x_:slope[x_], reverse=False)
    # print(ordered_clusters)
    # print([len(clusters[x_]) for x_ in ordered_clusters])
    markers = '(MUTE)', '(SPCH)', '(FAMA)', '(ATML1)', '(GASA9)' # , '(GC1)'
    with open(filename_output, 'w') as fo, open(filename_dump, 'w') as fo_dump:
        ecolumns = ['exp:{}'.format(x_) for x_ in expression.columns]
        fo_dump.write('Cluster\tGene\tLocation\t' + '\t'.join(display_data.columns) + '\t' + '\t'.join(ecolumns))
        fo.write('Cluter\tnum_ATAC\t' + '\t'.join(display_data.columns) + '\tnum_expr\t' + '\t'.join(ecolumns))
        fo_dump.write('\n')
        fo.write('\n')
        for cn in ordered_clusters:
            n_expr = 0    
            cluster_size =0
            expr_genes = []
            locations = []
            found = []

            for index in clusters[cn]:
                location = display_data.index[index]
                gene, loc = location.split(';', 1)
                fo_dump.write('{}\t{}\t{}'.format(cn, gene, loc))
                display_row = display_data.loc[location].values
                rawvalues = matrix.loc[location].values
                for val_ in rawvalues:#display_row:
                    fo_dump.write('\t{:.3f}'.format(val_))
                for m in markers:
                    if location.find(m) >= 0:
                        found.append(m.strip('()'))
                locations.append(location)
                cluster_size += 1
                ma = re.match("AT\\dG\\d+", location)
                if ma is not None:
                    agi = ma.group(0)
                    if agi in expressedgenes:
                        expr_genes.append(agi)
                        n_expr += 1
                        for val_ in expression.loc[agi].values:
                            fo_dump.write('\t{:.3f}'.format(val_))
                    else:
                        for i_ in range(expression.shape[1]):
                            fo_dump.write('\t.')
                elif expression is not None:
                    for i_ in range(expression.shape[1]):
                        fo_dump.write('\t.')
                fo_dump.write('\n')

                y.append(gene)
                # row_ = []
                # # z.append([max(vmin, min(vmax, x_)) for x_ in display_row])
                # for v in display_row:
                #     row_.append(max(vmin, min(vmax, v)))
                #     ostr += '\t{:.1f}'.format(v)
                z.append([max(vmin, min(vmax, v_)) for v_ in display_row])

            ostr = '{}\t{}'.format(cn, len(clusters[cn]))
            for v in np.mean(display_data.loc[locations], axis=0).values:
                ostr += '\t{:.3f}'.format(v)
            ostr += '\t{}'.format(n_expr)
            if n_expr > 0:
                emat = expression.loc[expr_genes]
                #emat[emat<0.1] = 0.1
                #emat = np.log2(emat)
                quvals = np.quantile(emat, [.25, .75], axis=0)
                # quvals = np.quantile(emat, [.25, .75], axis=0)
                # emeans = np.median(emat, axis=0)
                for x in range(quvals.shape[1]):#emeans:
                    ostr += '\t{:.3f},{:.3f}'.format(quvals[0][x], quvals[1][x])#x)
            else:
                ostr += ('\t.') * expression.shape[1]
            if len(found) > 0:
                ostr += '\t{}'.format(','.join(found))        
            else:
                ostr += '\t.'
            fo.write(ostr + '\n')
            if verbose:
                sys.stderr.write(ostr + '\n')
                # if verbose:
                #     fo.write()
                #     estr = '{}\t{}'.format(cn, n_expr)
                #     # for v_ in emat
                #     sys.stderr.write('cluster_size:{}, num_expr:{}, expr:{:.3f}'.format(cluster_size, n_expr, [int(x*100) for x in np.mean(emat, axis=0)]))
    import plotly.graph_objs as go
    import plotly, plotly.offline
    layout = go.Layout(yaxis=dict(autorange='reversed'))
    fig = go.Figure([go.Heatmap(x=display_data.columns, y=y, z=z)], layout)
    plotly.offline.plot(fig, filename=filename_graph)

if __name__ == '__main__':
    main()