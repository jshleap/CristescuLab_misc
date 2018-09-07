import pandas as pd
from io import StringIO
from subprocess import Popen, PIPE
from difflib import SequenceMatcher

def get_prefix(lis):
    first_el = lis.pop(0)
    for i in lis:
        match = SequenceMatcher(None, first_el, i)
        match = match.find_longest_match(0,len(first_el), 0, len(i))
        first_el = first_el[match.a: match.a + match.size]
    return first_el


line = 'seqkit stats *.fast*'
seqkit = Popen(line, shell=True, stdout=PIPE, stderr=PIPE)
o, e = seqkit.communicate()
df = pd.read_table(StringIO(o.decode()), delim_whitespace=True)
df['prefix'] = get_prefix(df.file)
df['type'] = df.file.apply(lambda x: '.'.join(x.split('.')[-2:]))
df = df.reindex(columns = 'prefix type num_seqs min_len avg_len max_len'.split())
df = pd.pivot_table(df, values='num_seqs  min_len  avg_len  max_len'.split(),
                    index='prefix', columns=['type'])
for i in 'min_len avg_len max_len'.split():
    sub = df.loc[:, i]
    pref = sub.index.tolist()[0]
    sub.to_csv('%s_%s.tsv' % (pref, i), sep='\t')
