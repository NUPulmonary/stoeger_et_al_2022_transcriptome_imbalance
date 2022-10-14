import multiprocessing
import os
import sys

import feather

import numpy as np
import pandas as pd

from scipy.stats import mannwhitneyu

from statsmodels.sandbox.stats.multicomp import multipletests

from access_biology_data import meta
from access_science_shared import standardizer

sys.path.append('./../src/')
from aging_tools import inout
from access_aging_data import companions

from narrative import nar181026_heat_confidence_genes


def load_181106_chached(ref_code):
    """
    Load cached data on p-values and fold changes of lung
    subpopulations

    Input:
        ref_code      reference gene code, e.g: 'orp' for
                        official gene symbol, reported in
                        literature, protein-coding

    Output:
        df_out_r      fold change of mean
        df_out_p      Mann-Whitney U
        df_out_padj   Mann-Whitney U (Bonferroni corrected)
    """

    p_folder = inout.get_internal_path(
        'dynamic/tstoeger/181106_cache_de_for_all_lung_subpopulations/')

    agg = []
    for j in np.arange(7):
        p = os.path.join(p_folder, 'df_out_r_{}.csv.gz').format(int(j))
        df = pd.read_csv(p)
        df = df.set_index('gene_ncbi')
        agg.append(df)

    df_out_r = pd.concat(agg, sort=True)

    agg = []
    for j in np.arange(7):
        p = os.path.join(p_folder, 'df_out_p_{}.csv.gz').format(int(j))
        df = pd.read_csv(p)
        df = df.set_index('gene_ncbi')
        agg.append(df)

    df_out_p = pd.concat(agg, sort=True)

    if ref_code is not None:
        ref_genes = standardizer.reference_genes(10090, ref_code)

        df_out_r = df_out_r[df_out_r.index.isin(ref_genes)]
        df_out_p = df_out_p[df_out_p.index.isin(ref_genes)]

    df_out_padj = pd.DataFrame(
        index=df_out_p.index,
        columns=df_out_p.columns)

    for c in df_out_p.columns:
        d = df_out_p.loc[:, c].dropna()
        padj = multipletests(d, method='bonferroni')[1]
        df_out_padj.loc[d.index, c] = padj

    df_out_padj = df_out_padj.astype(float)

    return df_out_r, df_out_p, df_out_padj


def get_differential_expression_18_over_4(sc_data, sc_meta):
    """
    Will perform differential gene expression using
    Mann-Whitney U

    """

    cell_types = sc_meta['cluster_label'].unique()
    genes_of_interest = sc_data.index.unique()

    sc_meta = sc_meta.reset_index()
    df_out_p = pd.DataFrame(index=genes_of_interest, columns=cell_types)
    df_out_r = pd.DataFrame(index=genes_of_interest, columns=cell_types)

    for cluster in cell_types:
        for gene in genes_of_interest:

            y_cells = sc_meta[(sc_meta['cluster_label'] == cluster) & (
                sc_meta['months'] == 4)]['cell_index']
            o_cells = sc_meta[(sc_meta['cluster_label'] == cluster) & (
                sc_meta['months'] == 18)]['cell_index']

            f = sc_data.index == gene
            if any(f):
                y = sc_data.loc[f, sc_data.columns.isin(y_cells)].values[0]
                o = sc_data.loc[f, sc_data.columns.isin(o_cells)].values[0]

            m_y = np.mean(y)
            m_o = np.mean(o)

            if m_y == 0:
                if m_o > 0:
                    f_c = np.inf
                else:
                    f_c = np.nan
            elif m_o == 0:
                f_c = -np.inf
            else:
                f_c = np.log2(m_o / m_y)

            if ((o.max() > 0) | (y.max() > 0)):
                _, pvalue = mannwhitneyu(y, o, alternative='two-sided')
                df_out_p.loc[gene, cluster] = pvalue

            df_out_r.loc[gene, cluster] = f_c

    return df_out_r, df_out_p


def prime_genes(all_de, tissues, pfus):
    """
    Returns double-replicate confirmed genes found in a given tissue under
    under a certain pfu with conditions reaching 100% precision

    Input:
        all_de      stacked differential expression table
        tissues     list of tissues
        pfus        list of pfus

    """
    d = cofirmed_sub_prime_genes(all_de, tissues, pfus, required_fraction=1)

    return d


def cofirmed_sub_prime_genes(all_de, tissues, pfus, required_fraction):
    """
    Returns double-replicate confirmed genes found in a given tissue under
    under a certain pfu with conditions reaching less than 100% precision

    Input:
        all_de      stacked differential expression table
        tissues     list of tissues
        pfus        list of pfus
        required_fraction  fraction required

    """

    focus = all_de[
        (all_de['tissue'].isin(tissues)) & (all_de['pfu'].isin(pfus))
    ].copy()

    settings = {
        'cis_series': ['both'],
        'trans_series': ['first', 'second'],
        'trans_logic': 'all',
        'significance_metric': 'log_padj',
        'trans_significance': -np.log10(0.05),
        'require_ubiquitous_detection': True,
        'cis_range_of_significance': np.arange(0, 20 + 0.01, 0.1)
    }

    fractions_confirmed, genes_confirmed = \
        nar181026_heat_confidence_genes.obtain_information_for_grouping(
            focus.groupby('condition'), settings)

    target = pd.Series(index=fractions_confirmed.columns)
    for t in list(fractions_confirmed.index)[::-1]:
        f = fractions_confirmed.loc[t, :] >= required_fraction
        target[f] = t

    target = target.dropna()

    d = pd.merge(        # above dropna, and  inner join filter for presence
        genes_confirmed,
        target.to_frame('necessary_threshold').rename_axis(
            'condition').reset_index()
    )

    d = d[d['log_padj'] > d['necessary_threshold']]

    return d


def putative_conversion_of_soupx_normalized_to_seurat_normalized(df_data):
    """
    According to hte best understanding of the reading of the code of
    SoupX and Seurat, they seem to be doing the normalization differently,
    despite the help suggesting that they did the same. Most
    critically, SoupX does not seem to divide by the total number of
    transcripts detected in a cell.

    Input:
        df_data   dataframe, rows: genes, columns: cellIDS
                    normalized by SoupX

    Output
        df_data   dataframe, rows: genes, columns: cellIDS,
                    mimicking Seurat normalization


    """

    df_data = df_data.apply(lambda x: 10**x - 1)
    total_per_cell = df_data.sum()
    df_data = df_data.divide(total_per_cell)
    df_data = df_data * 10000
    df_data = df_data.apply(lambda x: np.log10(x + 1))

    return df_data


def load_sc181101(dataset, enforce_ncbi_gene_ids):

    def _load_prior_cluster_assignment():

        p = inout.get_internal_path(
            'dynamic/tstoeger/180918_free_aginglungsoupx2_from_r/'
        )

        df_meta = _read_feather(os.path.join(p, 'meta.feather'))
        df_meta = df_meta.rename(columns={
            'feather_index': 'cell_index',
            'orig.ident': 'sample_id',
            'res.3': 'cluster_id',
            'percent.mito': 'percent_mito'
        })

        p = inout.get_box_path(
            'Aging_single_cell_analysis/02_Analysis_Naive_V1/s02_Subgroup_Analysis/cluster_cell_soupXV2.txt'
        )

        cluster_ids = pd.read_table(p, names=['imported'])
        cluster_ids['imported'] = cluster_ids['imported'].astype(str)
        cluster_ids.loc[:, 'cluster_id'] = cluster_ids['imported'].str.extract(
            '#C([0-9]*)', expand=False)
        cluster_ids.loc[:, 'cluster_label'] = cluster_ids[
            'imported'].str.extract(
            '#C[0-9]* (.*)', expand=False)
        cluster_ids = cluster_ids.drop('imported', axis=1)
        cluster_ids['cluster_id'] = cluster_ids['cluster_id'].astype(float)

        # fix inconsistent naming in cluster labels,
        # Confirmed by A. Misharin on Nov 2nd on Slack that
        # they refer to same entity
        f = cluster_ids['cluster_label'] == 'AT2 cells'
        cluster_ids.loc[f, 'cluster_label'] = 'AT2'

        df_meta['cluster_id'] = df_meta['cluster_id'].astype(float)
        df_meta = pd.merge(
            df_meta,
            cluster_ids,
            how='left',
            left_on='cluster_id',
            right_on='cluster_id',
        )

        df_original_cluster_ids = df_meta[[
            'cell_index', 'cluster_id', 'cluster_label']]

        return df_original_cluster_ids

    def _load_meta():
        p = inout.get_internal_path(
            'dynamic/tstoeger/181031_naive_after_soupx/'
        )

        df_meta = _read_feather(os.path.join(p, 'meta.feather'))
        df_meta = df_meta.rename(columns={
            'feather_index': 'cell_index',
            'orig.ident': 'sample_id',
            'res.3': 'cluster_id',
            'percent.mito': 'percent_mito'

        })

        # instead use original curation
        df_meta = df_meta.drop('cluster_id', axis=1)

        df_meta['nUMI'] = df_meta['nUMI'].astype(float)
        df_meta['nGene'] = df_meta['nGene'].astype(float)
        df_meta['percent_mito'] = df_meta['percent_mito'].astype(float)

        p = inout.get_internal_path(
            'dynamic/tstoeger/180918_free_aginglungsoupx2_from_r/'
        )

        age_labels = pd.read_csv(
            os.path.join(
                p,
                'labelling_of_age_extracted_from_comments_in_code_named_mastercode_r.csv'))

        df_meta = pd.merge(df_meta, age_labels, how='left')

        return df_meta

    def _load_data(dataset, enforce_ncbi_gene_ids):

        if dataset == 'raw':
            file = 'raw_data.feather'
        elif dataset == 'normalized':
            file = 'norm_data.feather'
        elif dataset == 'scaled':
            file = 'scaled_data.feather'
        else:
            raise EnvironmentError(
                'Dataset must be named raw, normalized, or scaled')

        p = os.path.join(
            inout.get_internal_path(
                'dynamic/tstoeger/181031_naive_after_soupx'), file)

        df_data = _read_feather(p).rename(columns={
            'feather_index': 'symbol_ncbi'
        })

        if enforce_ncbi_gene_ids is True:
            df_data = pd.merge(
                df_data,
                meta.gene_info(10090, ['gene_ncbi', 'symbol_ncbi'])).drop(
                'symbol_ncbi', axis=1).set_index('gene_ncbi')
        elif enforce_ncbi_gene_ids is False:
            pass
        else:
            raise ValueError(
                'enforce_ncbi_gene_ids must be either True or False')

        return df_data

    df_meta = _load_meta()
    df_original_cluster_ids = _load_prior_cluster_assignment()

    df_meta = pd.merge(df_meta, df_original_cluster_ids, how='left')
    df_meta = df_meta.set_index('cell_index')

    df_data = _load_data(dataset, enforce_ncbi_gene_ids)

    detected_cells = set(df_data.columns).intersection(set(df_meta.index))

    df_data = df_data.loc[:, df_data.columns.isin(detected_cells)]
    df_meta = df_meta[df_meta.index.isin(detected_cells)]

    return df_data, df_meta


def _read_feather(p):
        """
        Make dedicated loading routine for feather since synatax will
        change with future releases of feather
        """

        # cpus = multiprocessing.cpu_count()
        # if cpus > 1:
        #     cpus = cpus - 1

        # df = pd.read_feather(p, nthreads=1)

        df = feather.read_dataframe(p)

        return df
