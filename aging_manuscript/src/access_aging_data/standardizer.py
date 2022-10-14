import re
import sys

import pandas as pd

from access_biology_data import annotation


sys.path.append('./../src/')
from aging_tools import inout


def filter_reference_genes(
        ref_genes,
        taxon_ncbi, gene_filter='inflammation_et_al'):
    """
    Filters reference genes against genes that occ

    ref_genes: list of reference genes
    taxon_ncbi: ncbi taxonomy ID
    gene_filter: name of filter to be applied

    """

    if gene_filter == 'inflammation_et_al':

        go = annotation.go(taxon_ncbi)
        go['GO_term'] = go['GO_term'].str.lower()
        go = go[go['gene_ncbi'].isin(ref_genes)]
        go = go[['gene_ncbi', 'GO_term']].drop_duplicates()

        g = go[go['GO_term'].str.contains(
            'immune|stress|inflamm|infect'

        )]
        forbidden = g['gene_ncbi'].unique()
        filtered_genes = [x for x in ref_genes if x not in forbidden]

    elif gene_filter == 'neurons_et_al':

        go = annotation.go(taxon_ncbi)
        go['GO_term'] = go['GO_term'].str.lower()
        go = go[go['gene_ncbi'].isin(ref_genes)]
        go = go[['gene_ncbi', 'GO_term']].drop_duplicates()

        g = go[go['GO_term'].str.contains(
            'brain|neuro|nerv|cerebral|cortex|memory'

        )]
        forbidden = g['gene_ncbi'].unique()
        filtered_genes = [x for x in ref_genes if x not in forbidden]


    else:
        raise ValueError(
            'Currently chosen gene_filter, {}, not supported'.format(
                gene_filter))

    return filtered_genes


def reference_result_folder(version):
    """
    Will return the results subolder that shall be used for individual
    projects;

    Input:
        version         str, referes to a given definiotn, e.g.: 171020

    Output:
        reference_results      dict with run name as key, and results
                                subfolder as string
    """

    if version == '171020':
        reference_results = {
            '160728_NB501488_0018_AHFJJVBGXY': 'Result_20170824_171524',
            '160730_NB501488_0019_AHCK5NBGXY': 'Result_20170828_135335',
            '160801_NB501488_0020_AHW7KJBGXX': 'Result_20170829_155706',
            '160803_NB501488_0022_AHFWHLBGXY': 'Result_20170829_160233',
            '160804_NB501488_0023_AHFY22BGXY': 'Result_20170829_210405',
            '160805_NB501488_0024_AHFY5VBGXY': 'Result_20170830_084853',
            '160811_NB501488_0025_AHFY2GBGXY': 'Result_20170830_085741',
            '160817_NB501488_0026_AHHL7FBGXY': 'Result_20170830_103937',
            '160826_NB501488_0029_AHGH72BGXY': 'Result_20170830_112245',
            '160902_NB501488_0030_AHG7JMBGXY': 'Result_20170830_130329',
            '160903_NB501488_0031_AHHCYFBGXY': 'Result_20170830_162812',
            '160909_NB501488_0033_AHG723BGXY': 'Result_20170830_163103',
            '160922_NB501488_0038_AHLMWCBGXY': 'Result_20170830_183847',
            '161025_NB501488_0047_AHVLGJBGXY': 'Result_20170831_111430',
            '170323_NB501488_0080_AHVV5GBGXY': 'Result_20170828_110107',
            '170329_NB501488_0082_AH2VV5BGX2': 'Result_20170829_155538',
            '170330_NB501488_0083_AHFJY7BGX2': 'Result_20170829_160056',
            '170331_NB501488_0084_AHM2LJBGX2': 'Result_20170829_205929',
            '170405_NB501488_0086_AHMG5YBGX2': 'Result_20170829_210910',
            '170412_NB501488_0092_AHLYWVBGX2': 'Result_20170830_085414',
            '170413_NB501488_0093_AH5KGFBGX2': 'Result_20170830_090101',
            '170414_NB501488_0094_AH5FYKBGX2': 'Result_20170830_111707',
            '170420_NB501488_0097_AHFFHNBGX2': 'Result_20170830_112457',
            '170421_NB501488_0098_AHFG5MBGX2': 'Result_20170830_130530',
            '170424_NB501488_0099_AHCLK7BGX2': 'Result_20170830_162931',
            '170531_NB501488_0114_AHT2FVBGX2': 'Result_20170830_163156',
            '170717_NB501488_0130_AH5T2CBGX3': 'Result_20170830_184719',
            '170731_NB501488_0136_AHFGKTBGX3': 'Result_20171006_093924',
            '170921_NB501488_0147_AHN7CFBGX3': 'Result_20170928_144745'
        }
    else:
        raise ValueError(
            'Did not find the version definition {}'.format(
                version
            ))

    return reference_results


def reference_experimental_sets(keyword):
    """
    Looks up run names for selected set of experiments
    """

    if keyword == 'aging_map':
        run_names = [
            # '160728_NB501488_0018_AHFJJVBGXY', # initial AM run
            '160801_NB501488_0020_AHW7KJBGXX',
            # '160804_NB501488_0023_AHFY22BGXY',  # initial AT run
            '160805_NB501488_0024_AHFY5VBGXY',
            # '160730_NB501488_0019_AHCK5NBGXY',  # initial Lung run
            '160803_NB501488_0022_AHFWHLBGXY',
            '160811_NB501488_0025_AHFY2GBGXY',
            '160817_NB501488_0026_AHHL7FBGXY',
            # '160826_NB501488_0029_AHGH72BGXY',  # initial Brain run
            '160902_NB501488_0030_AHG7JMBGXY',
            '160903_NB501488_0031_AHHCYFBGXY',
            '160909_NB501488_0033_AHG723BGXY',
            '160922_NB501488_0038_AHLMWCBGXY',
            '161025_NB501488_0047_AHVLGJBGXY',
            '170323_NB501488_0080_AHVV5GBGXY',
            '170329_NB501488_0082_AH2VV5BGX2',
            '170330_NB501488_0083_AHFJY7BGX2',
            '170331_NB501488_0084_AHM2LJBGX2',
            # '170405_NB501488_0086_AHMG5YBGX2',  # Verificiation set
            '170412_NB501488_0092_AHLYWVBGX2',
            # '170413_NB501488_0093_AH5KGFBGX2',  # initial LI
            '170414_NB501488_0094_AH5FYKBGX2',
            '170420_NB501488_0097_AHFFHNBGX2',
            # '170421_NB501488_0098_AHFG5MBGX2',  # Verificiation set
            '170424_NB501488_0099_AHCLK7BGX2',
            '170531_NB501488_0114_AHT2FVBGX2',
            '170717_NB501488_0130_AH5T2CBGX3',


        ]

    return run_names


def sample_name_2_meta(sample_name, pattern_name):
    """
    Extracts meta-information about samples from their name.
    If applicable, numeric values are converted to a numeric data
    type.


    Not:
    The returned meta-information might contain different
    kinds of meta-information for different experiments.
    (for aging_map_default it will be mouse_id, tissue, months,
    pfu (plaque forming units of influenza), and replicate_id)

    Input:
        sample_name     str, name of sample
        pattern_name    str, name of scheme that will be used
                            to extract metainformation. Presently
                            supports: 'aging_map_default'
        meta            dict, meta-information as extracted from
                            sample_name.
    """

    m = dict()

    if pattern_name == 'aging_map_default':
        pattern = '^M([0-9]+)_([a-zA-Z0-9]*)_([0-9]{2})([MD])_F([0-9]+)_([0-9]+)$'
        ma = re.search(pattern, sample_name)

        if ma is None:
            raise ValueError(
                '''The sample_name {} can not
                be parsed by the scheme {}.'''.format(
                    sample_name, pattern))
        else:
            m['mouse_id'], m['tissue'], m[
                'age'], m['age_unit'], m['pfu'], m['replicate_id'] = ma.groups(1)

            for k in ['mouse_id', 'age', 'pfu', 'replicate_id']:
                m[k] = int(m[k])

    elif pattern_name == 'bowdish_tnf_condition':
        pattern = '^([A-Z]+)_([A-Z0-9]+)_([YO])_([0-9]+)$'
        ma = re.search(pattern, sample_name)

        if ma is None:
            raise ValueError(
                '''The sample_name {} can not
                be parsed by the scheme {}.'''.format(
                    sample_name, pattern))
        else:
            m['growth_condition'], m['tissue'], m[
                'age'], m['mouse_id'] = ma.groups(1)

            for k in ['mouse_id']:
                m[k] = int(m[k])

    elif pattern_name == 'aging_influenza_brain':

        pattern = '^(Brain)_([YO])_([NF])_([0-9]+)$'
        ma = re.search(pattern, sample_name)

        if ma is None:
            raise ValueError(
                '''The sample_name {} can not
                be parsed by the scheme {}.'''.format(
                    sample_name, pattern))
        else:
            m['tissue'], m['age'], m['influenza'], m['replicate_id'] = ma.groups(
                1)

            for k in ['replicate_id']:
                m[k] = int(m[k])

    else:
        raise ValueError(
            'Could not find any pattern_name called {}'.format(pattern_name))

    return m


def mouse_information():
    """
    Will return the information on individual mice
    """

    p = inout.get_internal_path(
        'datasets/general/sequencing/sync_from_google_sheets_171020/Aging_Map_samples.xlsx')

    df = pd.read_excel(p, sheetname='All mice')

    df['Mouse ID'] = df['Mouse ID'].copy().str.extract(
        'M(.*)', expand=False).apply(lambda x: float(x))
    df['Condition'] = df['Condition'].copy().str.extract(
        'F(.*)', expand=False).apply(lambda x: float(x))

    df = df.rename(columns={
        'Mouse ID': 'mouse_id',
        'Age (mo)': 'age_months',
        'Harvest date': 'harvest_date',
        'Batch': 'experimental_batch',
        'Condition': 'pfu',
        'Comments': 'comments'
    }).drop('comments', axis=1)

    return df
