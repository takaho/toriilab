import pandas as pd
import numpy as np
import sys, re

for fn in sys.argv[1:]:
    if fn.endswith('.tsv'):
        t = pd.read_csv(fn, sep='\t', index_col=0)
        s = t.copy()
        for c in s.columns:
            s[c] *= 1e6 / np.sum(t[c])
        fn_out = re.sub('\\.tsv$', '.tpm', fn)
        if fn == fn_out:
            continue
        with open(fn_out, 'w') as fo:
            fo.write('Location\t' + '\t'.join(s.columns) + '\n')
            locations = s.index
            values = s.values
            for i, l in enumerate(locations):
                fo.write(l)
                for v in values[i]:
                    fo.write('\t{:.3f}'.format(v))
                fo.write('\n')
            