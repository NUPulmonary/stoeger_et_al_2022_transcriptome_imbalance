import sys

import numpy as np
import pandas as pd

from natsort import natsorted


sys.path.append('./../src/')
from aging_tools import inout
from access_aging_data import companions


def load_rgrant_181123_detection(detection_column='pvalue'):
    """
    Wrapper. Loads datasets used for creating the figure, using dataset of
    rgrant 181123, but mimicking output of load_detection

        Input:
            detection_column  defaut: pvalue: column used to determine
                                detection of gene

        Output:
            all_de      differential expression resluts
    """

    all_de = companions.rgrant_181123_continuous_de()

    all_de = all_de.loc[:, [
        'gene_ncbi', 'de_unit', 'o_over_y', 'padj', 'pfu', 'pvalue', 'tissue'
    ]].copy()

    h = all_de['padj'].dropna().copy()
    f = all_de['padj'] == 0
    m = h[h > 0].min()
    all_de.loc[f, 'padj'] = m

    h = all_de['pvalue'].dropna().copy()
    f = all_de['pvalue'] == 0
    m = h[h > 0].min()
    all_de.loc[f, 'pvalue'] = m

    f = all_de['padj'].notnull()
    all_de.loc[f, 'log_padj'] = -np.log10(all_de.loc[f, 'padj'])
    all_de.loc[f, 'log_pvalue'] = -np.log10(all_de.loc[f, 'pvalue'])
    all_de.loc[:, 'is_detected'] = all_de.loc[:, detection_column].notnull()

    if all_de.groupby([
        'gene_ncbi', 'pfu', 'tissue', 'de_unit']).size(
    ).max() > 1:
        raise EnvironmentError('Some wrong redundancy.')

    all_de['condition'] = all_de[
        'tissue'] + '_' + all_de[
        'pfu'].astype(int).astype(str)

    return all_de


def load_detection(detection_column='pvalue'):
    """
    Loads several datasets used for creating the figure

    Input:
        detection_column  defaut: pvalue: column used to determine detection
                            of gene

    Output:
        all_de      differential expression resluts
        detection   presence of genes within individual comparisons
        mice_in_comparisons   numer of mice included in a given comparison
        triplicate_series     triplicate series

    """

    allowed_ages = [4, 9, 12, 18, 24]

    rohan = companions.tstoeger_181024_age_differential_expression(
        'inclusive_any_detection')

    thomas = companions.tstoeger_181024_age_differential_expression_by_batch(
        'inclusive_any_detection')

    # define exprimental series
    for j in [1, 3, 5]:
        thomas.loc[thomas['batch'] == j, 'series'] = 'first'
        thomas.loc[thomas['batch'] == (j + 1), 'series'] = 'second'

    interest = ['gene_ncbi', 'tissue', 'pfu',
                'younger', 'older', 'padj', 'series', 'pvalue', 'o_over_y']
    rohan.loc[:, 'series'] = 'both'

    all_de = pd.concat(
        [rohan[interest], thomas[interest]], sort=True).rename(
        columns={'series': 'de_unit'})

    f = (all_de['younger'].isin(allowed_ages)) & (
        all_de['older'].isin(allowed_ages))

    all_de = all_de[f]

    h = all_de['padj'].dropna().copy()
    f = all_de['padj'] == 0
    m = h[h > 0].min()
    all_de.loc[f, 'padj'] = m

    h = all_de['pvalue'].dropna().copy()
    f = all_de['pvalue'] == 0
    m = h[h > 0].min()
    all_de.loc[f, 'pvalue'] = m

    f = all_de['padj'].notnull()
    all_de.loc[f, 'log_padj'] = -np.log10(all_de.loc[f, 'padj'])
    all_de.loc[f, 'log_pvalue'] = -np.log10(all_de.loc[f, 'pvalue'])
    all_de.loc[:, 'is_detected'] = all_de.loc[:, detection_column].notnull()

    if all_de.groupby([
        'gene_ncbi', 'pfu', 'tissue', 'de_unit', 'younger', 'older']).size(
    ).max() > 1:
        raise EnvironmentError('Some wrong redundancy.')

    h = all_de.copy()
    detection = h[
        ['gene_ncbi',
         'tissue',
         'pfu',
         'is_detected',
         'younger',
         'older',
         'de_unit']].pivot_table(
        values='is_detected',
        index=['gene_ncbi', 'tissue', 'pfu', 'younger', 'older'],
        columns='de_unit'
    ).fillna(False).reset_index().rename_axis('', axis=1)

    p = inout.get_internal_path(
        'dynamic/tstoeger/181022_inclusive_any_detected_per_batch/sample_meta.csv')
    df_meta = pd.read_csv(p)

    df_meta = df_meta[df_meta['age'].isin(allowed_ages)]

    h = df_meta[[
        'age',
        'mouse_id',
        'pfu',
        'tissue',
        'experimental_batch']].drop_duplicates(
    ).groupby(
        ['age', 'pfu', 'tissue', 'experimental_batch']
    ).size().reset_index().rename(columns={0: 'mice'})

    # define exprimental series
    for j in [1, 3, 5]:
        h.loc[h['experimental_batch'] == j, 'series'] = 'first'
        h.loc[h['experimental_batch'] == (j + 1), 'series'] = 'second'

    if h.groupby(
            ['age', 'pfu', 'tissue', 'experimental_batch']).size().max() == 1:
        h = h.drop('experimental_batch', 1)
    else:
        raise ValueError(
            'There must only be a single entry per gene in one condition and batch.')

    h_b = pd.merge(
        h,
        h,
        left_on=['pfu', 'tissue', 'series'],
        right_on=['pfu', 'tissue', 'series'],
        suffixes=('_x', '_y')
    )
    h_b = h_b[(h_b['age_x'] < h_b['age_y'])]
    h_b = h_b.rename(columns={
        'age_x': 'younger',
        'age_y': 'older',
        'mice_x': 'mice_younger',
        'mice_y': 'mice_older'
    })
    mice_in_comparisons = h_b.copy()[[
        'tissue',
        'pfu',
        'younger',
        'older',
        'series',
        'mice_younger',
        'mice_older'
    ]]
    mice_in_comparisons['minimal_mice'] = mice_in_comparisons[[
        'mice_younger', 'mice_older'
    ]].min(1)

    he = mice_in_comparisons[[
        'tissue', 'pfu', 'younger', 'older', 'minimal_mice'
    ]].groupby(['tissue', 'pfu', 'younger', 'older']).agg(
        lambda x: sum(x == 3)
    ).rename(
        columns={'minimal_mice': 'triplicate_series'}).reset_index()
    he.loc[:, 'label'] = he['tissue'] + \
        '_' + he['pfu'].astype(int).astype(str)
    he.loc[:, 'has_two_triplicates'] = he['triplicate_series'] == 2
    triplicate_series = he

    # add condition for better readable exploratory analysis
    helper = all_de[['tissue', 'pfu', 'younger', 'older']].drop_duplicates()
    helper['condition'] = helper[
        'tissue'] + '_' + helper[
        'pfu'].astype(int).astype(str) + '_' + helper[
        'younger'].astype(int).astype(str) + '_' + helper[
        'older'].astype(int).astype(str)

    all_de = pd.merge(all_de, helper, how='left')

    return all_de, detection, mice_in_comparisons, triplicate_series


def obtain_information_for_grouping(grouping, settings):
    """
    Will determine the fraction of genes that also confirmed within
    another differential expression unit (used in the differential
    expression analysis, such as "first" or "second" batch or
    "both" in case that mice of both batches were pooled). Anticipates
    an inverse of p-value (higher numeric values corresponding to lower
    p-values)

    Input:
    - grouping            dataframe, containing 'de_unit', 'gene_ncbi',
                            and the defined significance_metric,
                            and 'is_detected', which would be grouped
    - settings dictionary containing:
        cis_series     list containing the series very p-value will increment
        trans_series   list containing the series very static p-value will
                                be considered
        trans_logic    either "all" or "any" (all or any sample in trans
                            needs to be significant)
        significance_metric  str, e.g.: "log_padj" column of signficance
                                    value that should be considered
        trans_significance  float, lower bound of significance; Note: not
                                included; ATTENTION:
                                The function expects the inverse of p-values,
                                so that  "higher" values of "significance
                                metric" would correspond to "higher"
                                significance (and initially smaller p-values)
        require_ubiquitous_detection  boolean, if TRUE the genes must be
                                            above detection threshold
                                            in cis_series and trans_series
        cis_range_of_significance  numpy range array, containing lower
                                        bounds to test

    Output:
    - fractions_confirmed      df, confirmed genes
    - genes_confirmed          genes that were confirmed in trans
    """

    agg_fraction_confirmed = {}
    agg_number_of_genes = {}
    agg_genes = {}

    for name, group in grouping:
        confirmed, amount, genes = obtain_confirmation_in_trans(
            group, settings)
        agg_fraction_confirmed[name] = confirmed
        agg_number_of_genes[name] = amount
        agg_genes[name] = genes

    collected = dict()
    has_gene = False
    for name in agg_genes.keys():
        genes = agg_genes[name]
        if genes.shape[0] > 0:
            genes.loc[:, 'condition'] = name
            collected[name] = genes
            has_gene = True

    if has_gene:
        genes_confirmed = pd.concat(
            collected, sort=False).reset_index(drop=True)
    else:
        genes_confirmed = []

    fractions_confirmed = pd.concat(agg_fraction_confirmed, 1)
    fractions_confirmed = fractions_confirmed.loc[:, natsorted(
        fractions_confirmed.columns)]

    return fractions_confirmed, genes_confirmed


def obtain_confirmation_in_trans(group, settings):
    """
    Will determine the fraction of genes that also confirmed within
    another differential expression unit (used in the differential
    expression analysis, such as "first" or "second" batch or
    "both" in case that mice of both batches were pooled). Anticipates
    an inverse of p-value (higher numeric values corresponding to lower
    p-values)

    Input:
    - group            dataframe, containing 'de_unit', 'gene_ncbi',
                            and the defined significance_metric,
                            and 'is_detected'
    - settings dictionary containing:
        cis_series     list containing the series very p-value will increment
        trans_series   list containing the series very static p-value will
                                be considered
        trans_logic    either "all" or "any" (all or any sample in trans
                            needs to be significant)
        significance_metric  str, e.g.: "log_padj" column of signficance
                                    value that should be considered
        trans_significance  float, lower bound of significance; Note: not
                                included; ATTENTION:
                                The function expects the inverse of p-values,
                                so that  "higher" values of "significance
                                metric" would correspond to "higher"
                                significance (and initially smaller p-values)
        require_ubiquitous_detection  boolean, if TRUE the genes must be
                                            above detection threshold
                                            in cis_series and trans_series
        cis_range_of_significance  numpy range array, containing lower
                                        bounds to test

    Output:
    - fraction_of_genes_confirmed_in_trans      series, containing the
                                                fraction of genes confirmed
                                                in trans at
                                                "cis_range_of_significance"
    - number_of_genes_beyond_bound              series, containing the number
                                                of cis genes beyond
                                                "cis_range_of_significance"
    """

    cis_series = settings['cis_series']
    trans_series = settings['trans_series']
    trans_logic = settings['trans_logic']
    significance_metric = settings['significance_metric']
    trans_significance = settings['trans_significance']
    require_ubiquitous_detection = settings['require_ubiquitous_detection']
    cis_range_of_significance = settings['cis_range_of_significance']

    if group[['gene_ncbi', 'de_unit']].groupby(
            ['gene_ncbi', 'de_unit']).size().max() > 1:
        raise EnvironmentError(
            'Genes must be unique within a given series and condition!')

    if len(cis_series) != 1:
        raise EnvironmentError(
            'The cis series must be a list with length of one.')

    amount_of_trans_series = len(trans_series)
    if amount_of_trans_series > 2:
        raise EnvironmentError(
            'The present code is only ready to process up to two trans series')

    if require_ubiquitous_detection:

        conditions_to_query = cis_series + trans_series
        h = group.loc[
            group['de_unit'].isin(conditions_to_query),
            ['gene_ncbi', 'de_unit', 'is_detected']].groupby(
            ['gene_ncbi']).agg(sum)
        always_detected_genes = set(
            h[h['is_detected'] == len(conditions_to_query)].index)
        group = group[group['gene_ncbi'].isin(always_detected_genes)]
    elif require_ubiquitous_detection == False:
        pass
    else:
        raise EnvironmentError(
            'require_ubiquitous_detection needs to be TRUE or FALSE.')

    cis = group[group['de_unit'].isin(cis_series)].copy()
    trans = group[group['de_unit'].isin(trans_series)].copy()

    trans.loc[:, 'is_significant'] = trans[
        significance_metric] > trans_significance

    significant_occurences_per_gene_in_trans = trans[
        trans['is_significant']]['gene_ncbi'].value_counts(
    )

    if trans_logic == 'all':
        if amount_of_trans_series == 2:
            genes_confirmed_in_trans = list(
                significant_occurences_per_gene_in_trans[
                    significant_occurences_per_gene_in_trans == 2
                ].index)
        elif amount_of_trans_series == 1:
            genes_confirmed_in_trans = list(
                significant_occurences_per_gene_in_trans[
                    significant_occurences_per_gene_in_trans == 1
                ].index)
    elif trans_logic == 'any':
        genes_confirmed_in_trans = list(
            significant_occurences_per_gene_in_trans[
                significant_occurences_per_gene_in_trans >= 1
            ].index)
    else:
        raise ValueError('trans_logic must either be all or any')

    cis.loc[:, 'confirmed_in_trans'] = cis['gene_ncbi'].isin(
        genes_confirmed_in_trans)

    fraction_of_genes_confirmed_in_trans = pd.Series(
        index=cis_range_of_significance)
    number_of_genes_beyond_bound = pd.Series(index=cis_range_of_significance)

    for bound_cis_significance in cis_range_of_significance:
        beyond_bound = cis.loc[
            cis[significance_metric] > bound_cis_significance, [
                'gene_ncbi', 'confirmed_in_trans']]
        fraction_of_genes_confirmed_in_trans[
            bound_cis_significance] = beyond_bound['confirmed_in_trans'].mean()
        number_of_genes_beyond_bound[
            bound_cis_significance] = beyond_bound['confirmed_in_trans'].shape[0]

    confirmed_genes = cis[cis['confirmed_in_trans']
                          ][['gene_ncbi', significance_metric]]

    return fraction_of_genes_confirmed_in_trans, number_of_genes_beyond_bound, confirmed_genes
