"""
**blast_processing.py
** Copyright (C) 2018  Jose Sergio Hleap

Parse a Blast output on format 6 with the following columns:
'qseqid sseqid pident evalue qcovs qlen length staxid stitle'. It assumes that
the stitle contains the lineage or the first two fields are species. If
lineage is in stitle, it will assume 7 taxon levels:
'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species'.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

E-mail: jose.hleaplozano@mcgill.ca

Python modules:
1. pandas
2. matplotlib
"""
import optparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.style.use('ggplot')
import pandas as pd
import numpy as np
from subprocess import Popen, PIPE

SIX = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']


def get_sps(line, coi=False):
    if coi:
        return line.split('|')[1]
    elif '|' in line:
        return ' '.join(line.split('|')[2].split()[:2])
    elif 'Fish' in line:
       return line.split()[0]
    else:
        return ' '.join(line.split()[1:3])


def parse_blast(fn, filters={}, top_n_hits=None, output_filtered=False,
                coi=False):
    """
    Parse a blast file, and filter it if required

    :param filters: Doctionary with the filters to include
    :param fn: Filename of the blast to parse
    :return: Dataframe with the information
    """

    df = pd.read_table(fn, sep='\t', header=None)
    if df.shape[1] > 6:
        names = 'qseqid sseqid pident evalue qcovs qlen length staxid stitle'
        names = names.split()
        by = 'evalue pident qcovs qlen length'.split()
        asc = [True, False, False, False, False]
    else:
        names = 'qseqid sseqid pident evalue qcovs stitle'.split()
        by = 'evalue pident qcovs'.split()
        asc = [True, False, False]
    df.rename(columns=dict(zip(range(len(names)), names)), inplace=True)
    if filters:
        query = ' & '.join(
            ['(%s > %d)' % (k, v) if k != 'evalue' else '(%s < %e)' % (k, v)
             for k, v in filters.items()])
        df = df.query(query)
    print(df.head())
    if top_n_hits is not None:
        args = dict(by=by, ascending=asc)
        df = df.sort_values(**args).groupby('qseqid').head(top_n_hits)
    if df.stitle.str.count(';').mean() > 2:
        # Lineage present, incorporate it
        ndf = df.stitle.str.split(' ', expand=True)
        ndf = ndf.loc[:1].str.split(';', expand=True)
        # Assume 7 level taxonomy
        ndf.rename(columns=dict(zip(range(1,8), SIX)), inplace=True)
        # Join the dataframes
        df = pd.concat([df, ndf], axis=1)
    else:
        # Assume that species is in the first two fields of stitle
        #def get_sps(x): return ' '.join(x.strip().split()[:2])
        df['species'] = df.stitle.apply(get_sps, coi=coi)
    if output_filtered:
        outfn = fn[:fn.rfind('.')]
        df.to_csv(outfn, sep='\t', index=False, header=False)
    return df


def get_reads_per_group(df, taxlevel='species', min_reads=10):
    """
    Get the number of reads per taxonomic level and the number of unique taxa
    per read

    :param min_reads: Minumum number of reads to retain group
    :param df: filtered blast dataframe
    :param taxlevel: taxonomic
    :return:
    """
    # Get number of reads per taxonomic group
    cou = df.groupby(taxlevel)['qseqid'].count()
    cou.to_csv('Numer_of_reads_in_%s.tsv' % taxlevel, sep='\t')
    # Get number of unique species per read
    re = df.groupby('qseqid')[taxlevel].nunique()
    re.to_csv('Numer_of_reads_in_%s.tsv' % taxlevel, sep='\t')
    # List number of unique species above the min_reads
    sps = cou[cou > min_reads].index.unique().to_series()
    sps.to_csv('List_unique_%s' % taxlevel, header=False, index=False)


def plot_tax(df, n, taxlevel='species', tax_for_pattern=None, pattern=None,
             suffix=None, min_reads=10):
    """
    Create a diversity bar chart
    :param df: Filtered blast datafrae
    :param n: number of records to retain
    :param taxlevel: Taxonomic level to display
    :param tax_for_pattern: Taxonomic level to restrict by pattern
    :param pattern: pattern to restrict the plotting
    :param suffix: Output suffix (before extension)
    """
    prefix = 'top%dhits_%s' % (n, taxlevel)
    if pattern is not None:
        df = df[df.loc[:, tax_for_pattern].str.contains(pattern)]
        prefix += '_%s' % pattern
    if suffix is not None:
        prefix += '_%s' % suffix
    cols = [taxlevel, 'qseqid']
    toplot = df.reindex(columns=cols)
    ntaxa = toplot.loc[:, taxlevel].nunique()
    color = matplotlib.cm.inferno_r(np.linspace(.0, 1., ntaxa))
    toplot = toplot.groupby([taxlevel]).count()
    toplot = toplot[toplot.qseqid > min_reads].unstack()
    fig, ax = plt.subplots()
    toplot.plot(kind='bar', stacked=False, color=color, ax=ax)
    ax.set_xlabel("Sites")
    ax.set_ylabel("Number of reads")
    plt.tight_layout()
    plt.savefig('%s.pdf' % prefix)
    plt.close()


def main(blast_file, pident=None, eval=None, qcov=None, qlen=None, length=None,
         output_filtered=False, taxlevel='species', min_reads=0, plot=False,
         tax_for_pattern=None, pattern=None, suffix_for_plot=None, ntop=None):
    filters = dict(zip('pident evalue qcovs qlen length'.split(),
                       [pident,   eval,   qcov,   qlen,   length]))
    filters = {k: v for k, v in filters.items() if v is not None}
    kwargs = dict(filters=filters, output_filtered=output_filtered,
                  top_n_hits=ntop)
    df = parse_blast(blast_file, **kwargs)
    get_reads_per_group(df, taxlevel=taxlevel, min_reads=min_reads)
    if plot:
        plot_tax(df, ntop, taxlevel=taxlevel, tax_for_pattern=tax_for_pattern,
                 pattern=pattern, suffix=suffix_for_plot, min_reads=min_reads)


if __name__ == '__main__':
    opts = optparse.OptionParser(usage='%prog [options] blast_file')
    opts.add_option('--output_filtered', '-o', action='store_true',
                    default=False, help=('Output a TSV with the filtered table'
                                         ' [default: %default]'))
    opts.add_option('--pident', '-p', action='store', type=float, default=None,
                    help='Minimum percent identity [default: %default]')
    opts.add_option('--eval', '-e', action='store', type=float, default=None,
                    help='Maximum evalue [default: %default]')
    opts.add_option('--qcov', '-q', action='store', type=float, default=None,
                    help='Minimum query coverage [default: %default]')
    opts.add_option('--qlen', '-Q', action='store', type=int, default=None,
                    help='Minimum query length [default: %default]')
    opts.add_option('--length', '-l', action='store', type=int, default=None,
                    help='Minimum alignment length [default: %default]')
    opts.add_option('--taxlevel', '-t', action='store', default='species',
                    help='Taxonomic level to display [default: %default]')
    opts.add_option('--min_reads', '-r', action='store', type=int, default=0,
                    help=('Minimum number of reads to retain group '
                          '[default: %default]'))
    opts.add_option('--plot', '-P', action='store_true', default=False,
                    help='Make a barchart with your group [default: %default]')
    opts.add_option('--tax_for_pattern', '-a', default=None,
                    help=('Parental taxonomic level to subset based on pattern'
                          ' [default: %default]'))
    opts.add_option('--pattern', '-b', default=None,
                    help=('Pattern to subset the tax_for_pattern with '
                          '[default: %default]'))
    opts.add_option('--suffix_for_plot', '-s', default=None,
                    help=('Suffix for plot (before extension) '
                          '[default: %default]'))
    opts.add_option('--ntop', '-n', action='store', type=int, default=None,
                    help='Number of hits per query [default: %default]')
    opts.add_option('--use_coi', '-c', action='store_true', default=False,
                    help=('If no special formating in teh database and using'
                          ' COI [default: %default]'))

    opt, arg = opts.parse_args()
    main(arg[0], opt.pident, opt.eval, opt.qcov, opt.qlen, opt.length,
         opt.output_filtered, opt.taxlevel, opt.min_reads, opt.plot,
         opt.tax_for_pattern, opt.pattern, opt.suffix_for_plot, opt.ntop)