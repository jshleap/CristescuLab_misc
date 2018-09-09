import pandas as pd
from io import StringIO
from subprocess import Popen, PIPE, check_output
from difflib import SequenceMatcher
from glob import glob
from joblib import Parallel, delayed
from tqdm import tqdm
import sys

COLS = 'prefix type num_seqs min_len avg_len max_len'.split()

def get_prefix(lis):
    first_el = lis.pop(0)
    for i in lis:
        match = SequenceMatcher(None, first_el, i)
        match = match.find_longest_match(0,len(first_el), 0, len(i))
        first_el = first_el[match.a: match.a + match.size]
    return first_el


def loop(d):
    line = 'seqkit stats %s*.fast*' % d
    #seqkit = Popen(line, shell=True, stdout=PIPE)
    #o, e = seqkit.communicate()
    files = glob('%s/*.fast*' % d)
    o = check_output(['seqkit', 'stats'] + files)
    df = pd.read_table(StringIO(o.decode()), delim_whitespace=True, thousands=',')
    df['prefix'] = get_prefix(df.file)
    df['type'] = df.file.apply(lambda x: '.'.join(x.split('.')[-2:]))
    df = df.reindex(columns=COLS)
    df = pd.pivot_table(df, values=COLS[2:], index=['prefix'], columns=['type'])
    return df


def writedf(df, col):
    df.loc[:,col].to_csv('%s_%s.tsv' % (sys.argv[1], col), sep='\t')


dirs = glob('%s*/' % sys.argv[1])
cpus = int(sys.argv[2])
alldfs = Parallel(n_jobs=cpus)(delayed(loop)(d) for d in tqdm(dirs, desc='Getting stats'))

df = pd.concat(alldfs)
cols = 'min_len avg_len max_len num_seqs'.split()
Parallel(n_jobs=cpus)(delayed(writedf)(df, col) for col in tqdm(cols, desc='Writing to file'))

