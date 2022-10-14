import glob
import os
import sys

import pandas as pd

from access_science_shared import mapper, utils


sys.path.append('./../src/')
from aging_tools import inout
from access_aging_data import sequencing


def tstoeger_190427_gtex():
    """
    Aggregate differential gene expression for GTEx
    created in pulmonary style; here: all tissues (except
    cells) for male and female, if at least 2 samples
    per condition (those were processed for DE)
    """

    p_male = inout.get_internal_path(
        'dynamic/tstoeger/190427_gtex_m/DE/Flu_imac/')
    df_male = _load_gtex_subset(p_male)
    df_male.loc[:, 'gender'] = 'male'

    p_female = inout.get_internal_path(
        'dynamic/tstoeger/190427_gtex_f/DE/Flu_imac/')
    df_female = _load_gtex_subset(p_female)
    df_female.loc[:, 'gender'] = 'female'

    df = pd.concat([df_male, df_female]).reset_index(drop=True)
    return df


def _load_gtex_subset(p):
    """
    Imports a GTEXx dataset as created on 190427 (mimicking
    pulmonary's style of DE)
    """

    df = pd.DataFrame(
        data={'fn': os.listdir(p)})

    df['tissue'] = df['fn'].str.extract('(.*)_pfu', expand=False)

    df['older'] = df['fn'].str.extract(
        '.*_ages_([0-9])', expand=False).astype(float)
    df['younger'] = df['fn'].str.extract(
        '.*_ages_[0-9] ([0-9])_DE', expand=False).astype(float)

    f = df['tissue'].str.startswith('Cells')   # do not load cells
    df = df[~f]

    agg = []
    for j in df.index:
        pa = os.path.join(p, df.loc[j, 'fn'])
        tissue = df.loc[j, 'tissue']
        older = df.loc[j, 'older']
        younger = df.loc[j, 'younger']

        d = pd.read_csv(
            pa,
            usecols=['Symbol', 'log2FoldChange', 'pvalue', 'padj']).rename(
            columns={
                'Symbol': 'gene_ensembl',
                'log2FoldChange': 'o_over_y'
            })

        d.loc[:, 'tissue'] = tissue
        d.loc[:, 'younger'] = younger
        d.loc[:, 'older'] = older

        agg.append(d)

    df = pd.concat(agg)

    df = mapper.gene_ensembl_2_gene_ncbi_unambiguously(
        df, 9606).reset_index().sort_values(
        ['tissue', 'younger', 'older']
    )
    return df


def rgrant_181123_continuous_de():
    '''
    Loads DE expression analysis performed by Rogan Grant
    on 181123, which modeled age as a continous variable.


    Output:
        all_de          stacked dataframe, which with de_unit to
                            indicate if "first", or "second"
                            replicate batch, or "both"
    '''

    pa = glob.glob(
        inout.get_internal_path(
            'datasets/rgrant/181019_meta_analysis/181123_amap_reanalysis_split_by_tissue/*.csv'))

    m = pd.Series(pa).to_frame('full_path')

    def get_fn(x):
        _, x = os.path.split(x)
        return x

    m.loc[:, 'fn'] = m['full_path'].apply(lambda x: get_fn(x))

    def get_short_fn(x):
        x = x.replace('182311AgingMap_reanalysis_', '')
        x = x.replace('182311181115_AgingMap_reanalysis_', '')
        x = x.replace('182411181115_AgingMap_reanalysis_', '')
        x = x.replace('_age_dge', '')
        x = x.replace('.csv', '')
        x = x.replace('182411AgingMap_reanalysis_', '')
        return x

    m['short'] = m['fn'].apply(lambda x: get_short_fn(x))

    f = m['short'].str.startswith('A_')
    m.loc[f, 'de_unit'] = 'first'

    f = m['short'].str.startswith('B_')
    m.loc[f, 'de_unit'] = 'second'

    f = m['short'].str.startswith('all_')
    m.loc[f, 'de_unit'] = 'both'

    m['de_unit'].value_counts()

    f = m['short'].str.contains('_150pfu_')
    m.loc[f, 'pfu'] = 150

    f = m['short'].str.contains('_10pfu_')
    m.loc[f, 'pfu'] = 10

    # Slack: Rogan: And correct, anything without a
    # “XXXpfu” flag is the the output from yesterday (no flu) (edited)
    f = m.loc[:, 'pfu'].isnull()
    m.loc[f, 'pfu'] = 0

    m['tissue'] = m['short'].str.extract('^.*_(.*)$')

    if m.groupby(['de_unit', 'pfu', 'tissue']).size().max() > 1:
        raise ValueError('Meta is not unambiguous')

    m = m[~m['tissue'].isin(['age'])]

    agg = []
    for _, r in m.iterrows():
        rr = r

        df = pd.read_csv(
            rr['full_path'],
            usecols=[   # only load ensemble genes to avoid ambiguity
                'ensembl_gene_id',
                'baseMean',
                'log2FoldChange',
                'lfcSE',
                'stat',
                'pvalue',
                'padj'])

        df.loc[:, 'tissue'] = rr['tissue']
        df.loc[:, 'pfu'] = rr['pfu']
        df.loc[:, 'de_unit'] = rr['de_unit']
        agg.append(df)

    all_de = pd.concat(agg, sort=True)

    all_de = all_de[all_de['pvalue'].notnull()]

    all_de = all_de.drop_duplicates()

    all_de = mapper.gene_ensembl_2_gene_ncbi_unambiguously(
        all_de.rename(columns={'ensembl_gene_id': 'gene_ensembl'}),
        taxon_id=10090).reset_index()

    return all_de


def tstoeger_190415_pfu_differential_expression_to_4_months(kind):
    """
    Loads differential gene expression analysis
    for distinct pfu doses. All relative to 0 PFU
    in 4 months old animals.

    Will map to genes to NCBI gene identifiers.

    """

    if 'inclusive_any_detection':

        df = _load_pfu_de(
            'datasets/tstoeger/190415_flu_relative_to_4_month_no_infection/age_groups.csv.gz'
        )

    else:
        raise ValueError('kind inclusive_any_detection')

    return df


def tstoeger_190327_pfu_differential_expression(kind):
    """
    Loads differential gene expression analysis
    for distinct pfu doses

    Will map to genes to NCBI gene identifiers.

    """

    if 'inclusive_any_detection':

        df = _load_pfu_de(
            'datasets/tstoeger/190327_pooled_differential_expression_flu_any_detected/age_groups.csv.gz'
        )

    else:
        raise ValueError('kind inclusive_any_detection')

    return df


def tstoeger_181024_age_differential_expression(kind):
    """
    Loads differential gene expression analysis
    for distinct ages

    Will map to genes to NCBI gene identifiers.

    """

    if kind == 'inclusive':

        df = _load_age_de(
            'datasets/tstoeger/181024_pooled_differential_expression_inclusive_de/age_groups.csv.gz'
        )

    elif kind == 'inclusive_any_detection':

        df = _load_age_de(
            'datasets/tstoeger/181024_pooled_differential_expression_inclusive_de_any_detected/age_groups.csv.gz'
        )

    else:
        raise ValueError('kind must be inclusive or inclusive_any_detection')

    return df


def _load_pfu_de(p):
    """
    Loads differential gene expression

    Will map to genes to NCBI gene identifiers.

    """

    columns_to_use = [
        'gene_ncbi',
        'tissue',
        'age',
        'lower',
        'higher',
        'h_over_l',
        'pvalue',
        'padj']

    df = pd.read_csv(
        inout.get_internal_path(
            p),
        usecols=columns_to_use,
    )

    output_order = [
        'gene_ncbi',
        'tissue',
        'age',
        'lower',
        'higher',
        'h_over_l',
        'pvalue',
        'padj'
    ]
    df = df[output_order]

    df = df.dropna(subset=['gene_ncbi'])

    return df


def _load_age_de(p):
    """
    Loads differential gene expression analysis performed by Rohan
    Verma on 180105 for distinct ages, and packaged by Thomas
    Stoeger on 180924.

    Will map to genes to NCBI gene identifiers.

    """

    columns_to_use = [
        'gene_ncbi',
        'tissue',
        'pfu',
        'youngest',
        'oldest',
        'o_over_y',
        'pvalue',
        'padj']

    df = pd.read_csv(
        inout.get_internal_path(
            p),
        usecols=columns_to_use,
    ).rename(
        columns={'youngest': 'younger', 'oldest': 'older'})

    output_order = [
        'gene_ncbi',
        'tissue',
        'pfu',
        'younger',
        'older',
        'o_over_y',
        'pvalue',
        'padj'
    ]
    df = df[output_order]

    df = df.dropna(subset=['gene_ncbi'])

    return df


def tstoeger_181024_age_differential_expression_by_batch(kind):
    """
    Loads differential gene expression analysis on aging map,
    separately analysed per replicate.

    Will map to genes to NCBI gene identifiers.

    """

    if kind == 'inclusive':

        df = _load_age_de_by_batch(
            'datasets/tstoeger/181024_pool_differential_expression_inclusive_de_by_batch/age_groups.csv.gz'
        )

    elif kind == 'inclusive_any_detection':

        df = _load_age_de_by_batch(
            'datasets/tstoeger/181024_pool_differential_expression_inclusive_any_detected_by_batch/age_groups.csv.gz'
        )
    else:
        raise ValueError('kind must be inclusive or inclusive_any_detection')

    return df


def _load_age_de_by_batch(p):
    """
    Loads differential exrpression by age batch

    Input:
        p           path to differential expression

    """

    columns_to_use = [
        'gene_ncbi',
        'tissue',
        'pfu',
        'youngest',
        'oldest',
        'o_over_y',
        'pvalue',
        'padj',
        'batch']

    df = pd.read_csv(
        inout.get_internal_path(
            p),
        usecols=columns_to_use,
    ).rename(
        columns={'youngest': 'younger', 'oldest': 'older'})

    output_order = [
        'gene_ncbi',
        'tissue',
        'pfu',
        'younger',
        'older',
        'o_over_y',
        'pvalue',
        'padj',
        'batch'
    ]
    df = df[output_order]

    df = df.dropna(subset=['gene_ncbi'])
    return df


def sal_protein_181008():
    """
    Loads age-dependent protein expression of AM and AT2 cells
    shared by Salvatore Loguercio on October 3rd.

    Will use median in case that multiple genes map to same protein.

    """

    p = inout.get_internal_path('proteomics/PPG_AM_AT_0818_MSSTATS.xlsx')

    sheets_of_interest = (
        'PPG_AM_Old_AM_Yng',
        'PPG_AT_Old_AT_Yng',
    )

    sheet_meta = {
        'PPG_AM_Old_AM_Yng': {
            'younger': 4, 'older': 18, 'tissue': 'AM', 'pfu': 0},
        'PPG_AT_Old_AT_Yng': {
            'younger': 4, 'older': 18, 'tissue': 'AT2', 'pfu': 0},
    }

    agg = []

    for sheet in sheets_of_interest:
        df = pd.read_excel(p, sheet_name=sheet)

        df = utils.split_text_to_multiple_rows(df[[
            'pvalue',
            'adj.pvalue',
            'Protein',
            'log2FC'
        ]], 'Protein', ';').rename(
            columns={
                'Protein': 'protein_uniprot',
                'adj.pvalue': 'padj',
                'log2FC': 'o_over_y'
            }).dropna()

        df = mapper.uniprot_protein_2_gene_ncbi(df, 'median')

        df = df.reset_index()
        df.loc[:, 'tissue'] = sheet_meta[sheet]['tissue']
        df.loc[:, 'pfu'] = sheet_meta[sheet]['pfu']
        df.loc[:, 'younger'] = sheet_meta[sheet]['younger']
        df.loc[:, 'older'] = sheet_meta[sheet]['older']
        agg.append(df)

    df = pd.concat(agg)

    return df


def flu_brain_181002_differential_expression():
    """
    Loads differential gene expression analysis of "flu brain" data,
    which uses counts as computed by pipeline of agingmap. Will
    perform DE according to template used by Rohan Verma for
    agingmousemap.

    Will map to genes to NCBI gene identifiers.

    """

    columns_to_use = [
        'gene_ncbi',
        'tissue',
        'pfu',
        'youngest',
        'oldest',
        'o_over_y',
        'pvalue',
        'padj']

    df = pd.read_csv(
        inout.get_internal_path(
            'datasets/tstoeger/181002_pooled_differential_expression_flu_brain_as_in_agingmap/age_groups.csv.gz'),
        usecols=columns_to_use,
    ).rename(
        columns={'youngest': 'younger', 'oldest': 'older'})

    output_order = [
        'gene_ncbi',
        'tissue',
        'pfu',
        'younger',
        'older',
        'o_over_y',
        'pvalue',
        'padj'
    ]
    df = df[output_order]

    df = df.dropna(subset=['gene_ncbi'])

    return df


def bowdish_181002_differential_expression():
    """
    Loads differential gene expression analysis of Bowdish data,
    which uses counts as computed by pipeline of agingmap. Will
    perform DE according to template used by Rohan Verma for
    agingmousemap.

    Will map to genes to NCBI gene identifiers.

    """

    columns_to_use = [
        'gene_ncbi',
        'condition',
        'tissue',
        'pfu',
        'youngest',
        'oldest',
        'o_over_y',
        'pvalue',
        'padj']

    df = pd.read_csv(
        inout.get_internal_path(
            'datasets/tstoeger/181002_pooled_differential_expression_bowdish_as_in_agingmap/age_groups.csv.gz'),
        usecols=columns_to_use,
    ).rename(
        columns={'youngest': 'younger', 'oldest': 'older'})

    output_order = [
        'gene_ncbi',
        'tissue',
        'pfu',
        'condition',
        'younger',
        'older',
        'o_over_y',
        'pvalue',
        'padj'
    ]
    df = df[output_order]

    f = df['tissue'] == 'WL'
    df.loc[f, 'tissue'] = 'Lung'

    df = df.dropna(subset=['gene_ncbi'])

    return df


def tstoeger_181013_age_differential_expression_by_batch():
    """
    Loads differential gene expression analysis on aging map,
    separately analysed per batch.

    Will map to genes to NCBI gene identifiers.

    """

    columns_to_use = [
        'gene_ncbi',
        'tissue',
        'pfu',
        'youngest',
        'oldest',
        'o_over_y',
        'pvalue',
        'padj',
        'batch']

    df = pd.read_csv(
        inout.get_internal_path(
            'datasets/tstoeger/181013_pooled_differential_expression_by_batch/age_groups.csv.gz'),
        usecols=columns_to_use,
    ).rename(
        columns={'youngest': 'younger', 'oldest': 'older'})

    output_order = [
        'gene_ncbi',
        'tissue',
        'pfu',
        'younger',
        'older',
        'o_over_y',
        'pvalue',
        'padj',
        'batch'
    ]
    df = df[output_order]

    df = df.dropna(subset=['gene_ncbi'])

    return df


def tstoeger_181012_age_differential_expression_by_batch():
    """
    Loads differential gene expression analysis on aging map,
    separately analysed per replicate.

    Will map to genes to NCBI gene identifiers.

    """

    columns_to_use = [
        'gene_ncbi',
        'tissue',
        'pfu',
        'youngest',
        'oldest',
        'o_over_y',
        'pvalue',
        'padj',
        'batch']

    df = pd.read_csv(
        inout.get_internal_path(
            'datasets/tstoeger/181012_pooled_differential_expression_by_replicate/age_groups.csv.gz'),
        usecols=columns_to_use,
    ).rename(
        columns={'youngest': 'younger', 'oldest': 'older'})

    output_order = [
        'gene_ncbi',
        'tissue',
        'pfu',
        'younger',
        'older',
        'o_over_y',
        'pvalue',
        'padj',
        'batch'
    ]
    df = df[output_order]

    df = df.dropna(subset=['gene_ncbi'])

    return df


def rverma_180924_age_differential_expression():
    """
    Loads differential gene expression analysis performed by Rohan
    Verma on 180105 for distinct ages, and packaged by Thomas
    Stoeger on 180924.

    Will map to genes to NCBI gene identifiers.

    """

    columns_to_use = [
        'gene_ncbi',
        'tissue',
        'pfu',
        'youngest',
        'oldest',
        'o_over_y',
        'pvalue',
        'padj']

    df = pd.read_csv(
        inout.get_internal_path(
            'datasets/tstoeger/180924_pool_differential_expression/age_groups.csv.gz'),
        usecols=columns_to_use,
    ).rename(
        columns={'youngest': 'younger', 'oldest': 'older'})

    output_order = [
        'gene_ncbi',
        'tissue',
        'pfu',
        'younger',
        'older',
        'o_over_y',
        'pvalue',
        'padj'
    ]
    df = df[output_order]

    df = df.dropna(subset=['gene_ncbi'])

    return df


def rverma_180924_flu_differential_expression():
    """
    Loads differential gene expression analysis performed by Rohan
    Verma on 180105 for distinct flu doses, and packaged by Thomas
    Stoeger on 180924.

    Will map to genes to NCBI gene identifiers.

    """

    columns_to_use = [
        'gene_ncbi',
        'tissue',
        'age',
        'lowest_pfu',
        'highest_pfu',
        'h_over_l',
        'pvalue',
        'padj']

    df = pd.read_csv(
        inout.get_internal_path(
            'datasets/tstoeger/180924_pool_differential_expression/flu_groups.csv.gz'),
        usecols=columns_to_use,
    ).rename(
        columns={'lowest_pfu': 'lower', 'highest_pfu': 'higher'})

    output_order = [
        'gene_ncbi',
        'tissue',
        'age',
        'lower',
        'higher',
        'h_over_l',
        'pvalue',
        'padj'
    ]
    df = df[output_order]

    df = df.dropna(subset=['gene_ncbi'])

    return df


def bowdish_180718_differential_expression():
    """
    Loads differential gene expression of Bowdish data set.
    Shared by A. Misharin on 180718 to T. Stoeger through
    Email. Note that this dataset only contains signficiant genes
    and that preprocessing has been done with Ceto rather
    than with the bioinformatics pipline used for aging
    mouse map at that given point.
    """

    p_base = inout.get_internal_path(
        'datasets/general/complementary/180718_lung_mail_from_amisharin/Bowdish_pairwise_result_ZR_0623')

    agg = []
    potential_comparisons = ['SPF', 'GF', 'TNF', 'O', 'Y']
    for tissue in ['AM', 'AT2', 'WL']:
        for condition in potential_comparisons:
            for dividend in potential_comparisons:
                for divisor in potential_comparisons:

                    p = os.path.join(
                        p_base,
                        '{}_result'.format(tissue),
                        'Annotated_{}_{}_{}vs{}_count.csv'.format(
                            tissue,
                            condition,
                            dividend,
                            divisor,
                        )
                    )

                    if os.path.isfile(p):
                        df = pd.read_csv(p)

                        df = df[[
                            'Row.names',
                            'logFC',
                            'PValue',
                            'FDR',
                        ]].rename(columns={'Row.names': 'gene_ensembl'})

                        df.loc[:, 'tissue'] = tissue
                        df.loc[:, 'condition'] = condition
                        df.loc[:, 'dividend'] = dividend
                        df.loc[:, 'divisor'] = divisor

                        agg.append(df)

    df = pd.concat(agg)
    df = mapper.gene_ensembl_2_gene_ncbi_unambiguously(
        df, taxon_id=10090).reset_index()
    df = df.sort_values(['tissue', 'condition', 'dividend', 'divisor', 'FDR'])
    df = df.reset_index(drop=True)

    return df


def bowdish_180921_differential_expression():
    """
    Loads differential gene expression of Bowdish data set.
    Shared by Z. Ren on 180921 to T. Stoeger through
    Email. Note that preprocessing has been done with Ceto rather
    than with the bioinformatics pipline used for aging
    mouse map at that given point. It further only contains
    one condition, "AM".
    """

    p_base = inout.get_internal_path(
        'datasets/general/complementary/180921_bowdish_mail_from_zren')

    agg = []
    potential_comparisons = ['SPF', 'GF', 'TNF', 'O', 'Y']
    for tissue in ['AM']:  # ['AM', 'AT2', 'WL']:
        for condition in potential_comparisons:
            for dividend in potential_comparisons:
                for divisor in potential_comparisons:

                    p = os.path.join(
                        p_base,
                        '{}_{}{}vs{}_wholegene_list.csv'.format(
                            tissue,
                            condition,
                            dividend,
                            divisor,
                        )
                    )

                    if os.path.isfile(p):
                        df = pd.read_csv(p).rename(
                            columns={'Unnamed: 0': 'gene_ensembl'})

                        df = df[[
                            'gene_ensembl',
                            'logFC',
                            'logCPM',
                            'LR',
                            'PValue',
                            'FDR',
                        ]]

                        df.loc[:, 'tissue'] = tissue
                        df.loc[:, 'condition'] = condition
                        df.loc[:, 'dividend'] = dividend
                        df.loc[:, 'divisor'] = divisor
                        agg.append(df)

    df = pd.concat(agg)
    df = mapper.gene_ensembl_2_gene_ncbi_unambiguously(
        df, taxon_id=10090).reset_index()
    df = df.sort_values(['tissue', 'condition', 'dividend', 'divisor', 'FDR'])
    df = df.reset_index(drop=True)

    return df


def bowdish_180921f_differential_expression():
    """
    Loads differential gene expression of Bowdish data set.
    Almelolar macrophages shared by Z. Ren on 180921 to
    T. Stoeger through Email. Alveolar Type II, and whole lung
    shared by Z. Ren on 180928 to T. Stoeger trough Slack.
    Note that preprocessing has been done with Ceto rather
    than with the bioinformatics pipline used for aging
    mouse map at that given point.
    """

    p_base_am = inout.get_internal_path(
        'datasets/general/complementary/180921_bowdish_mail_from_zren')

    p_base_at2_and_wl = inout.get_internal_path(
        'datasets/general/complementary/180928_bowdish_slack_from_zren/result')

    agg = []
    potential_comparisons = ['SPF', 'GF', 'TNF', 'O', 'Y']
    for tissue in ['AM', 'AT2', 'WL']:
        for condition in potential_comparisons:
            for dividend in potential_comparisons:
                for divisor in potential_comparisons:

                    if tissue == 'AM':
                        p = os.path.join(
                            p_base_am,
                            '{}_{}{}vs{}_wholegene_list.csv'.format(
                                tissue,
                                condition,
                                dividend,
                                divisor,
                            )
                        )
                    else:
                        p = os.path.join(
                            p_base_at2_and_wl,
                            '{}_{}{}vs{}_wholegene_list.csv'.format(
                                tissue,
                                condition,
                                dividend,
                                divisor,
                            )
                        )

                    if os.path.isfile(p):
                        df = pd.read_csv(p).rename(
                            columns={'Unnamed: 0': 'gene_ensembl'})

                        df = df[[
                            'gene_ensembl',
                            'logFC',
                            'logCPM',
                            'LR',
                            'PValue',
                            'FDR',
                        ]]

                        df.loc[:, 'tissue'] = tissue
                        df.loc[:, 'condition'] = condition
                        df.loc[:, 'dividend'] = dividend
                        df.loc[:, 'divisor'] = divisor
                        agg.append(df)

    df = pd.concat(agg)
    df = mapper.gene_ensembl_2_gene_ncbi_unambiguously(
        df, taxon_id=10090).reset_index()
    df = df.sort_values(['tissue', 'condition', 'dividend', 'divisor', 'FDR'])
    df = df.reset_index(drop=True)

    # Enforce nomenclature consistent with aging map
    f = df['tissue'] == 'WL'
    df.loc[f, 'tissue'] = 'Lung'

    return df


### START OF DEPRECATED FUNCTIONS ####

def rverma_180105_age_differential_expression():
    """
    Loads differential gene expression analysis performed by Rohan
    Verma on 180105 for distinct ages.

    Will map to genes to NCBI gene identifiers.

    """

    raise ValueError(
        'This function has been deprecated. The direction of the fold change '
        'o_over_y is wrong for some of the samples. There is no reason '
        'to doubt the p-values.'
    )

    _, _, df_genes = sequencing.load_cached_aging_map(
        dataset_name='aging_map_tmm_180105',
        unambiguous_to_entrez=True,
        as_entrez=True
    )

    df_diff_expression = pd.read_csv(
        inout.get_internal_path(
            'datasets/tstoeger/180117_pooled_differential_expression/age_groups.csv.gz'),
    )

    df_diff_expression = pd.merge(
        df_diff_expression,
        df_genes[['gene_ensembl', 'gene_ncbi']]).drop(
        'gene_ensembl', 1).set_index('gene_ncbi').reset_index()

    return df_diff_expression

### END OF DEPRECATED FUNCTIONS ####
