
import os
import sys

import numpy as np
import pandas as pd

from access_biology_data import meta, relations
from access_science_shared import inout as inout_scisci
from access_science_shared import mapper, utils

sys.path.append('./../src/')
from aging_tools import inout


def anage():
    """
    Database containing lifespan related information for several organisns
    """
    p = inout_scisci.get_path('anage')
    p = os.path.join(p, 'anage_data.txt')
    df = pd.read_table(p)

    return df


def bodymap():
    """
    Load bodymap, a large-scale survey of Rat
    """

    p = inout.get_internal_path(
        'datasets/general/resources/bodymap/Ratbodymap_AgePattern_v2.txt')

    df_map = pd.read_table(p)

    p = inout.get_internal_path(
        (
            'datasets/general/resources/bodymap/'
            'SEQC_Ratbodymap_GeneInfo_20130813updated.txt'))
    df_gene_info = pd.read_table(p, usecols=['GeneSymbol', 'EntrezID'])
    df_body_map = pd.merge(df_gene_info, df_map, how='inner')

    df_body_map = df_body_map[['EntrezID', 'Tissue', 'Age', 'Pattern']].rename(
        columns={
            'EntrezID': 'gene_ncbi',
            'Tissue': 'tissue',
            'Age': 'age_comparison',
            'Pattern': 'age_pattern'
        })

    df_body_map['age_pattern'] = df_body_map['age_pattern'].replace(
        {
            'M': 'maintain',
            'U': 'up',
            'D': 'down'})

    df_body_map['tissue'] = df_body_map['tissue'].replace(
        {
            'Adrenal Gland': 'Adrenal'
        })

    return df_body_map


def fleischer_2018():
    """
    Loads FPKM data from Fleischer et al. 2018, and
    mapps to genes

    """

    p = inout.get_internal_path(
        (
            'datasets/general/resources/publications/'
            'fleischer_2018/GSE113957_fpkm.txt.gz'
        )
    )

    df = pd.read_table(p)

    df = df.drop('Copies', 1)
    df = df.drop('chr', 1)
    df = df.drop('start', 1)
    df = df.drop('end', 1)
    df = df.drop('strand', 1)
    df = df.drop('Annotation/Divergence', 1)
    df = df.drop('Length', 1)

    df = df.rename(columns={'Transcript ID': 'rna_ncbi'})
    df = mapper.rna_ncbi_2_gene_ncbi(df, 'sum')

    return df


def deMagalhaes_2009(direction):
    """
    Loads a meta-analysis by de Magalhaes et al. 2009, who
    found a signature of some genes that
    appear to be stable makers of aging at transcripitonal level
    in many tissues and experiments.

    input:
        direction   str, overexpression or underexpression

    """

    p = inout.get_internal_path(
        (
            'datasets/general/resources/publications/'
            'deMagalhaes_2009/supplementary_tables.xls'
        )
    )

    if direction == 'overexpression':

        df = pd.read_excel(
            p,
            sheet_name='Genes_overexpressed',
            header=14).rename(columns={'EntrezGeneID': 'gene_ncbi'})

    elif direction == 'underexpression':

        df = pd.read_excel(
            p,
            sheet_name='Genes_underexpressed',
            header=9).rename(columns={'EntrezGeneID': 'gene_ncbi'})

    else:
        raise ValueError(
            'direction must either be overexpression or underexpression')

    df = df.rename(columns={
        'EntrezGeneID': 'gene_ncbi'
    })

    return df


def hagr_mapped_summary(target_taxon):
    """
    Loads HAGR (which includes the AnAge database) and mappes
    documents associations with age to the target_taxon.

    Unites two databases:
        genage_models:  model organisms (only preserved in this
            function, if documented influence)
        longevity: genetic associations on human (only preserved
            in this function, if significant influence)
    Processing:
        - mapping through homologene
        - if applicable: add records on aging on genes, which
            are not part of homologene (e.g.: not conserved)

    Input:
        - target_taxon:     ncbi taxonomy identified of target species

    Output:
        - hagr:     df: 'gene_ncbi', 'influence' (pro-longevity, anti-
            longevity, human_association), 'taxon_ncbi_source' (source
            of annotation)

    """

    p_hagr = inout_scisci.get_path('hagr')

    p = os.path.join(p_hagr, 'models_genes/genage_models.csv')
    df = pd.read_csv(p)
    df = df[['entrez gene id', 'longevity influence']]
    df = df.rename(columns={'entrez gene id': 'gene_ncbi',
                            'longevity influence': 'influence'})
    present_taxa = [
        6239, 10090, 559292, 7227, 10036, 5145, 284812, 9606,
        7955, 515849, 6238]

    agg = []
    for t in present_taxa:
        gi = meta.gene_info(taxon_id=t, usecols=['gene_ncbi', 'taxon_ncbi'])
        agg.append(gi)
    gi = pd.concat(agg)
    models = pd.merge(df, gi)
    models['influence'] = models['influence'].replace(
        {'Unclear': 'Unannotated'})  # corresponds to ~30 records
    models['influence'] = models['influence'].replace(
        {'Pro-Longevity': 'pro_longevity', 'Anti-Longevity': 'anti_longevity'})
    models = models[models['influence'].isin(
        ['pro_longevity', 'anti_longevity'])]

    p = os.path.join(p_hagr, 'longevity_genes/longevity.csv')
    df = pd.read_csv(p)
    df = df[['Association', 'Gene(s)']]
    df = df.astype(str)
    df = df.rename(columns={'Gene(s)': 'symbol_ncbi'})
    df = utils.split_text_to_multiple_rows(df, 'symbol_ncbi', '\,')
    df = df.drop_duplicates()

    gi = meta.gene_info(9606, usecols=['gene_ncbi', 'symbol_ncbi'])
    h = pd.merge(df, gi)
    h = h[h['Association'] == 'significant']
    h.loc[:, 'taxon_ncbi'] = 9606
    h = h.rename(columns={'Association': 'influence'})
    h['influence'] = h['influence'].copy().replace(
        {'significant': 'human_association'})

    joint = pd.concat([
        models,
        h[['gene_ncbi', 'influence', 'taxon_ncbi']]
    ]
    ).rename(columns={
        'gene_ncbi': 'gene_ncbi_source',
        'taxon_ncbi': 'taxon_ncbi_source'})

    hg = relations.homologene()
    master = pd.merge(joint, hg[['homologene_group', 'gene_ncbi']].rename(
        columns={'gene_ncbi': 'gene_ncbi_source'}))

    f = hg['taxon_ncbi'] == target_taxon
    if any(f):
        master = pd.merge(master, hg.loc[f, ['homologene_group', 'gene_ncbi']])
        patch = joint[joint['taxon_ncbi_source'] == target_taxon].copy()
        patch['gene_ncbi'] = patch['gene_ncbi_source']
        master = master[['gene_ncbi_source', 'influence',
                         'taxon_ncbi_source', 'gene_ncbi']]
        master = pd.concat([master, patch]).drop_duplicates()

    master = master[['gene_ncbi', 'influence', 'taxon_ncbi_source']]
    master = master.drop_duplicates().sort_values([
        'gene_ncbi', 'influence', 'taxon_ncbi_source']).reset_index(drop=True)

    return master


def angelidis_2018(dataset):
    """
    Loads reference datasets from Anglidis, 2018, paper

    dataset:
        simple_protein
        simple_simulated_bulk_rna
    """

    if dataset == 'simple_protein':
        df = _load_angelidis_2018_simple_protein()
    elif dataset == 'simulated_bulk_rna':
        df = _load_angelidis_2018_simple_simulated_bulk_rna()

    return df


def ori_2015(dataset, map_to_mouse=True):
    """
    Loads differential expression of Ori et al. 2015, rat brain and liver.
    Young: 6 months
    Old: 24 months
    """

    p_base = os.path.join(
        inout_scisci.get_path('publications'),
        'ori2015')

    if dataset == 'transcripts':

        p = os.path.join(p_base, 'mmc2.xlsx')

        agg = []
        for tissue in ['liver', 'brain']:

            df = pd.read_excel(
                p, sheet_name=tissue, header=1)

            df = df.rename(columns={   # bypass typo
                'Is change in transcription significant?': 'significant',
                'Is change in transcription significant? ': 'significant',
            })

            df = df[[
                'EnsTranscID',
                'logFC(transcription)',
                'adj.pval(transcription)',
                'significant'
            ]].rename(columns={
                'EnsTranscID': 'rna_ensembl',
                'logFC(transcription)': 'logFC_OoverY',
                'adj.pval(transcription)': 'padj',
            })

            df = pd.concat([
                mapper.rna_ensembl_2_gene_ncbi(df[
                    ['rna_ensembl', 'logFC_OoverY', 'padj']], 'median'),
                mapper.rna_ensembl_2_gene_ncbi(df[
                    ['rna_ensembl', 'significant']], 'any')
            ],
                axis=1)

            df = df.reset_index()
            df.loc[:, 'tissue'] = tissue.capitalize()
            agg.append(df)

        df = pd.concat(agg)
        df = df.sort_values(['tissue', 'gene_ncbi'])

    elif dataset == 'proteins':

        p = os.path.join(p_base, 'mmc3.xlsx')

        agg = []
        for tissue in ['liver', 'brain']:

            df = pd.read_excel(
                p,
                sheet_name='{}_proteins_quantified'.format(tissue),
                header=1)

            df = df.rename(columns={   # bypass typo
                'Is change significant?': 'significant',
                'Is change significant? ': 'significant',
            })

            df = df[[
                'EnsGeneID',
                'logFC',
                'pval',
                'significant'
            ]].rename(columns={
                'EnsGeneID': 'gene_ensembl',
                'logFC': 'logFC_OoverY'
            })

            df = pd.concat([
                mapper.gene_ensembl_2_gene_ncbi_unambiguously(
                    df[['gene_ensembl', 'logFC_OoverY', 'pval']],
                    taxon_id=10116
                ),
                mapper.gene_ensembl_2_gene_ncbi_unambiguously(
                    df[['gene_ensembl', 'significant']],
                    taxon_id=10116
                )
            ],
                axis=1)

            df = df.reset_index()
            df.loc[:, 'tissue'] = tissue.capitalize()
            agg.append(df)

        df = pd.concat(agg)
        df = df.sort_values(['tissue', 'gene_ncbi'])

    else:
        raise ValueError(
            'Only supports transcripts or proteins as dataset.')

    if map_to_mouse:

        hg = relations.homologene()

        hg_mapper = pd.merge(
            hg.loc[hg['taxon_ncbi'] == 10116, [
                'homologene_group', 'gene_ncbi']],
            hg.loc[hg['taxon_ncbi'] == 10090, [
                'homologene_group', 'gene_ncbi']],
            left_on='homologene_group',
            right_on='homologene_group',
            suffixes=('', '_m')
        )[['gene_ncbi', 'gene_ncbi_m']].drop_duplicates()

        df = pd.merge(df, hg_mapper).drop('gene_ncbi', 1).rename(
            columns={'gene_ncbi_m': 'gene_ncbi'}).groupby(
            ['tissue', 'gene_ncbi']).agg(np.median).reset_index()

    return df


def zahn_2007():
    """
    Loads data from Zahn et al. (AGEMAP), which used
    microarrays to identify age-related genes across several
    time points for several tissues of C57BL/6 mice.

    Only gene symbols or gene synonyms which are unambiguous
    within NCBI (2017) are mapped. In case that several
    records map to the same NCBI (Entrez) gene, the
    median value is considered.


    outputs:
        df_correlation    a correlation metric for age-relatedness
                            introduced by authors
        df_pvalues        note that they used 0.001 as threshold
                            for their paper
    """

    p = inout_scisci.get_path(
        'publications',
        'zahn2007/journal.pgen.0030201.st001.XLS'
    )
    df = pd.read_excel(p, header=5)

    f = df == 'Age Coef.a'
    df[f] = 'Age Coef.'

    f = df == 'pb <'
    df[f] = 'p <'

    df = df.rename(columns={
        'Unnamed: 0': 'unigene',
        'Unnamed: 1': 'gene_symbol',
        'Unnamed: 2': 'GO',
        'Adrenals.1': 'Adrenals',
        'Bone Marrow.1': 'Bone Marrow',
        'Cerebellum.1': 'Cerebellum',
        'Cerebrum.1': 'Cerebrum',
        'Eye.1': 'Eye',
        'Gonads.1': 'Gonads',
        'Heart.1': 'Heart',
        'Hippocampus.1': 'Hippocampus',
        'Kidney.1': 'Kidney',
        'Liver.1': 'Liver',
        'Lung.1': 'Lung',
        'Muscle.1': 'Muscle',
        'Spinal Cord.1': 'Spinal Cord',
        'Spleen.1': 'Spleen',
        'Striatum.1': 'Striatum',
        'Thymus.1': 'Thymus'})
    df = df.drop('unigene', 1)
    df = df.drop('GO', 1)

    from access_biology_data import meta

    from access_science_shared import utils

    gi = meta.gene_info(10090)[['gene_ncbi', 'symbol_ncbi', 'Synonyms']]

    gi = utils.split_text_to_multiple_rows(gi, 'Synonyms', '\|')

    gi = pd.concat([
        gi[['gene_ncbi', 'symbol_ncbi']].rename(
            columns={'symbol_ncbi': 'symbol'}),
        gi[['gene_ncbi', 'Synonyms']].rename(columns={'Synonyms': 'symbol'})
    ])
    gi = gi[gi['symbol'] != '-'].dropna()
    gi = gi.drop_duplicates()  # all symbol gene combinations
    gi = gi.drop_duplicates(subset='symbol', keep=False)  # discard ambiguous

    # get data relating to p-values (and meta)
    is_not_coef = df.loc[0, :] != 'Age Coef.'
    dp = pd.merge(df.loc[:, is_not_coef], gi,
                  left_on='gene_symbol', right_on='symbol', how='inner')
    dp = dp.drop('symbol', 1)
    dp = dp.drop('gene_symbol', 1)
    dp = dp.drop_duplicates()
    dp = dp.astype(float)
    dp = dp.groupby('gene_ncbi').agg(np.median)

    # get data relating to correlations (and meta)
    is_not_pvalue = df.loc[0, :] != 'p <'
    dc = pd.merge(df.loc[:, is_not_pvalue], gi,
                  left_on='gene_symbol', right_on='symbol', how='inner')
    dc = dc.drop('symbol', 1)
    dc = dc.drop('gene_symbol', 1)
    dc = dc.drop_duplicates()
    dc = dc.astype(float)
    dc = dc.groupby('gene_ncbi').agg(np.median)

    # format outptu
    dp = dp.sort_index()
    dc = dc.sort_index()

    # sanity check
    if not all(dp.index == dc.index):
        raise AssertionError(
            'Something got lost during processing. Need to check code')

    return dc, dp


def _load_angelidis_2018_simple_protein():

    gene_info = meta.gene_info(10090)
    p = inout.get_internal_path(
        (
            'datasets/tstoeger/180718_angelidis_pre_publications/'
            '172963_0_supp_617204_phm79s.xlsx'))

    df = pd.read_excel(p)

    df = df[[
        "Gene names",
        "Student's T-test Difference_old/young [log2]",
        "-Log Student's T-test p-value 24m",
        "Student's T-test Significant 24m"
    ]]

    df = df.rename(columns={
        "Gene names": 'symbol_ncbi',
        "Student's T-test Difference_old/young [log2]": 'fold_change',
        "-Log Student's T-test p-value 24m": "p_value",
        "Student's T-test Significant 24m": 'significant'
    })

    df['significant'] = df['significant'] == '+'

    a = pd.merge(df, gene_info[['gene_ncbi', 'symbol_ncbi']])
    a = a.groupby(['gene_ncbi', 'symbol_ncbi']).agg(np.median)
    a = a.reset_index()

    if a['gene_ncbi'].value_counts().max() > 1:
        raise ValueError('something wrong')
    else:
        df = a.copy().drop('symbol_ncbi', axis=1).set_index('gene_ncbi')

    df = df.sort_values('fold_change', ascending=False)
    return df


def _load_angelidis_2018_simple_simulated_bulk_rna():

    gene_info = meta.gene_info(10090)
    p = inout.get_internal_path(
        (
            'datasets/tstoeger/180718_angelidis_pre_publications/'
            '172963_0_supp_617205_p1m791.xlsx'))
    df = pd.read_excel(p)

    df = df[[
        "Gene name",
        "log2FoldChange",
        "pvalue",
        "padj"
    ]]

    df = df.rename(columns={
        "Gene name": 'symbol_ncbi',
        "log2FoldChange": 'fold_change',
        "pvalue": "p_value",
        "padj": 'padj'
    })

    a = pd.merge(df, gene_info[['gene_ncbi', 'symbol_ncbi']])
    a = a.groupby(['gene_ncbi', 'symbol_ncbi']).agg(np.median)
    a = a.reset_index()
    if a['gene_ncbi'].value_counts().max() > 1:
        raise ValueError('something wrong')
    else:
        df = a.copy().drop('symbol_ncbi', axis=1).set_index('gene_ncbi')

    df = df.sort_values('fold_change', ascending=False)

    df = df.dropna(how='all')

    return df
