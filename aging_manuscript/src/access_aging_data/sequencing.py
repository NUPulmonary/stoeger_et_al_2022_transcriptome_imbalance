import os
import sys

import numpy as np
import pandas as pd

sys.path.append('./../src/')
from aging_tools import inout

from access_aging_data import standardizer as aging_standardizer


def load_cached_aging_map(dataset_name, unambiguous_to_entrez, as_entrez):
    """
    Will load a cached version of aging map:

    Input:
        - data_set_name         str, name of dataset
        - unambiguous_to_entrez default: False; If True: only genes with
                                    1:1 of ensembl and Entrez Gene ID are shown
        - as_entrez             default: False; If True: will return Entrez
                                    as Index


    recommendations:
        aging_map_qc_171031     raw counts
        aging_map_tmm_180105    TMM normalized per condition, after QC of AM/SB
                                (171219f_overview_of_metrics_v4_180101.xlsx)

    """

    # Input Check:
    if (unambiguous_to_entrez == False) and (as_entrez == True):
        raise AssertionError(
            'unambigous_to_entrez must be True if as_entrez is True.'
            'This ensures that there is no misunderstaning in returned'
            'Index values.'
        )

    # Sample metadata
    if dataset_name == '181028_inclusive_tmm':
        p = inout.get_internal_path(
            'datasets/general/sequencing/aging_map_qc_171031/df_meta.csv')
        df_meta = pd.read_csv(p).rename(
            columns={'Unnamed: 0': 'sample_name'}).set_index(
            'sample_name')
        df_meta = df_meta.drop('passed_qc', axis=1)   # ignore QC by Sophia Liu

        # Gene meta info
        p = inout.get_internal_path(
            'datasets/general/sequencing/aging_map_qc_171031/df_genes.csv')
        df_genes = pd.read_csv(p).set_index('gene_ensembl')

        p = inout.get_internal_path(
            'datasets/general/complementary/transcript_info_171103.csv')
        df_transcript_info = pd.read_csv(p)
        df_transcript_info.loc[:, 'present_in_genbank_info'] = True

        df_genes_extended = pd.merge(
            df_genes.reset_index(),
            df_transcript_info,
            left_on='gene_ensembl',
            right_on='gene_ensembl',
            how='left')
        df_genes_extended['present_in_genbank_info'] = df_genes_extended[
            'present_in_genbank_info'].fillna(False)
        df_genes = df_genes_extended

        # Count data (TMM normalized per condition)
        p_datasets = inout.get_internal_path(
            'dynamic/tstoeger/181028_inclusive_tmm/normalized_out')

        # # get conditions (note: MuscSat 150 PFU is absent)
        conditions_to_consider = df_meta[['pfu', 'tissue']].dropna(
        ).drop_duplicates().reset_index(drop=True)

        agg = []
        for _, content in conditions_to_consider.iterrows():
            pfu = content['pfu']
            tissue = content['tissue']

            p = os.path.join(
                p_datasets,
                '{}_f{}_counts.csv'.format(tissue, int(pfu))
            )

            df = pd.read_csv(p).rename(
                columns={
                    'Unnamed: 0': 'gene_ensembl'
                }).set_index('gene_ensembl')

            agg.append(df)
        df_counts = pd.concat(agg, axis=1, join='outer', verify_integrity=True)

        f = df_counts.index.isin(df_genes['gene_ensembl'])
        condition_cis = (f.sum() / len(f)) == 1
        f = df_genes['gene_ensembl'].isin(df_counts.index)
        condition_trans = (f.sum() / len(f)) == 1
        if not (condition_cis & condition_trans):
            raise ValueError(
                'Different genes listed in counts and gene meta-info.')

        # Arrange output
        df_meta = df_meta.loc[df_counts.columns, :]
        df_meta = df_meta.sort_values(['tissue', 'pfu', 'mouse_id'])
        df_counts = df_counts.loc[:, df_meta.index]

    elif dataset_name == 'aging_map_tmm_180105':

        p = inout.get_internal_path(
            'datasets/general/sequencing/aging_map_qc_171031/df_meta.csv')
        df_meta = pd.read_csv(p).rename(
            columns={'Unnamed: 0': 'sample_name'}).set_index(
            'sample_name')
        df_meta = df_meta.drop('passed_qc', axis=1)   # ignore QC by Sophia Liu

        # Gene meta info
        p = inout.get_internal_path(
            'datasets/general/sequencing/aging_map_qc_171031/df_genes.csv')
        df_genes = pd.read_csv(p).set_index('gene_ensembl')

        p = inout.get_internal_path(
            'datasets/general/complementary/transcript_info_171103.csv')
        df_transcript_info = pd.read_csv(p)
        df_transcript_info.loc[:, 'present_in_genbank_info'] = True

        df_genes_extended = pd.merge(
            df_genes.reset_index(),
            df_transcript_info,
            left_on='gene_ensembl',
            right_on='gene_ensembl',
            how='left')
        df_genes_extended['present_in_genbank_info'] = df_genes_extended[
            'present_in_genbank_info'].fillna(False)
        df_genes = df_genes_extended

        # Count data (TMM normalized per condition)
        p_datasets = inout.get_internal_path(
            'datasets/general/sequencing/180105_TMM_normalized')

        # # get conditions (note: MuscSat 150 PFU is absent)
        conditions_to_consider = df_meta[['pfu', 'tissue']].dropna(
        ).drop_duplicates().reset_index(drop=True)
        # f = (
        #     conditions_to_consider['pfu'] == 150) & (
        #     conditions_to_consider['tissue'] == 'MuscSat')
        # conditions_to_consider = conditions_to_consider.loc[~f, :]

        # read in count data
        agg = []
        for _, content in conditions_to_consider.iterrows():
            pfu = content['pfu']
            tissue = content['tissue']

            p = os.path.join(
                p_datasets,
                '{}_f{}_counts.csv'.format(tissue, int(pfu))
            )

            df = pd.read_csv(p).rename(
                columns={
                    'Unnamed: 0': 'gene_ensembl'
                }).set_index('gene_ensembl')

            agg.append(df)
        df_counts = pd.concat(agg, axis=1, join='outer', verify_integrity=True)

        # # tidy up counts: the genes list entries that aren't genes
        # f = df_counts.index.isin([
        #     '__no_feature',
        #     '__ambiguous',
        #     '__too_low_aQual',
        #     '__not_aligned',
        #     '__alignment_not_unique'
        # ])
        # df_counts = df_counts.loc[~f, :]

        f = df_counts.index.isin(df_genes['gene_ensembl'])
        condition_cis = (f.sum() / len(f)) == 1
        f = df_genes['gene_ensembl'].isin(df_counts.index)
        condition_trans = (f.sum() / len(f)) == 1
        if not (condition_cis & condition_trans):
            raise ValueError(
                'Different genes listed in coutns and gene meta-info.')

        # Arrange output
        df_meta = df_meta.loc[df_counts.columns, :]
        df_meta = df_meta.sort_values(['tissue', 'pfu', 'mouse_id'])
        df_counts = df_counts.loc[:, df_meta.index]

    elif dataset_name == 'aging_map_qc_171031':
        """
        This corresponds to the raw counts

        """

        p = inout.get_internal_path(
            'datasets/general/sequencing/aging_map_qc_171031/df_counts.csv')
        df_counts = pd.read_csv(p).set_index('gene_ensembl')
        df_counts.columns.name = 'sample_name'

        p = inout.get_internal_path(
            'datasets/general/sequencing/aging_map_qc_171031/df_meta.csv')
        df_meta = pd.read_csv(p).rename(
            columns={'Unnamed: 0': 'sample_name'}).set_index(
            'sample_name')

        p = inout.get_internal_path(
            'datasets/general/sequencing/aging_map_qc_171031/df_genes.csv')
        df_genes = pd.read_csv(p).set_index('gene_ensembl')

        # f = np.array([x in df_counts.columns for x in df_meta.index])
        # if f.sum() == df_counts.shape[1]:
        df_meta = df_meta.loc[df_counts.columns, :]

        p = inout.get_internal_path(
            'datasets/general/complementary/transcript_info_171103.csv')
        df_transcript_info = pd.read_csv(p)
        df_transcript_info.loc[:, 'present_in_genbank_info'] = True

        df_genes_extended = pd.merge(
            df_genes.reset_index(),
            df_transcript_info,
            left_on='gene_ensembl',
            right_on='gene_ensembl',
            how='left')
        df_genes_extended['present_in_genbank_info'] = df_genes_extended[
            'present_in_genbank_info'].fillna(False)
        df_genes = df_genes_extended

    elif dataset_name == 'edgeR_TMM_normalized_171031':

        #     edgeR_TMM_normalized_171031:
        # x only contains sampels without influenza
        # - df_counts: TMM normalized
        # - df_meta: Metainfo of samples
        # - df_genes: Gene info as obtained by HTSeq, and further
        #     summary statistics on genes from Genbank / NCI
        #     note that there can be multiple records / rows for
        #     one gene_ensembl

        p = inout.get_internal_path(
            'datasets/general/sequencing/edgeR_TMM_normalized_171031/TMM_normalized.csv')
        df_counts = pd.read_csv(p).rename(
            columns={'Unnamed: 0': 'gene_ensembl'}).set_index(
            'gene_ensembl')
        df_counts.columns.name = 'sample_name'

        p = inout.get_internal_path(
            'datasets/general/sequencing/aging_map_qc_171031/df_meta.csv')
        df_meta = pd.read_csv(p).rename(
            columns={'Unnamed: 0': 'sample_name'}).set_index(
            'sample_name')

        p = inout.get_internal_path(
            'datasets/general/sequencing/aging_map_qc_171031/df_genes.csv')
        df_genes = pd.read_csv(p).set_index('gene_ensembl')

        # f = np.array([x in df_counts.columns for x in df_meta.index])
        # if f.sum() == df_counts.shape[1]:
        df_meta = df_meta.loc[df_counts.columns, :]

        p = inout.get_internal_path(
            'datasets/general/complementary/transcript_info_171103.csv')
        df_transcript_info = pd.read_csv(p)
        df_transcript_info.loc[:, 'present_in_genbank_info'] = True

        df_genes_extended = pd.merge(
            df_genes.reset_index(),
            df_transcript_info,
            left_on='gene_ensembl',
            right_on='gene_ensembl',
            how='left')
        df_genes_extended['present_in_genbank_info'] = df_genes_extended[
            'present_in_genbank_info'].fillna(False)
        df_genes = df_genes_extended

    elif dataset_name == 'batchcorrection_TMM_combat_171031':

        agg = []
        for t in [
            'AM',
            'AT2',
            'Blood',
            'Brain',
            'Cerebellum',
            'GutEP',
            'Heart',
            'Kidney',
            'Lung',
            'MoDC'
        ]:
            p = inout.get_internal_path(
                'datasets/general/sequencing/batchcorrection_TMM_combat_171031/{}_TMM.csv'.format(t))
            agg.append(pd.read_csv(p).rename(
                columns={'Unnamed: 0': 'gene_ensembl'}).set_index('gene_ensembl'))

        df_counts = pd.concat(agg, axis=1)
        df_counts.index.name = 'gene_ensembl'

        p = inout.get_internal_path(
            'datasets/general/sequencing/aging_map_qc_171031/df_meta.csv')
        df_meta = pd.read_csv(p).rename(
            columns={'Unnamed: 0': 'sample_name'}).set_index(
            'sample_name')

        p = inout.get_internal_path(
            'datasets/general/sequencing/aging_map_qc_171031/df_genes.csv')
        df_genes = pd.read_csv(p).set_index('gene_ensembl')

        # f = np.array([x in df_counts.columns for x in df_meta.index])
        # if f.sum() == df_counts.shape[1]:
        df_meta = df_meta.loc[df_counts.columns, :]

        p = inout.get_internal_path(
            'datasets/general/complementary/transcript_info_171103.csv')
        df_transcript_info = pd.read_csv(p)
        df_transcript_info.loc[:, 'present_in_genbank_info'] = True

        df_genes_extended = pd.merge(
            df_genes.reset_index(),
            df_transcript_info,
            left_on='gene_ensembl',
            right_on='gene_ensembl',
            how='left')
        df_genes_extended['present_in_genbank_info'] = df_genes_extended[
            'present_in_genbank_info'].fillna(False)
        df_genes = df_genes_extended

    elif dataset_name == 'TMM_normalized_bytissue_171113_combatch':

        agg = []
        for t in [
            'Adrenal',
            'AM',
            'AT2',
            'BAT',
            'Blood',
            'Brain',
            'Cerebellum',
            'Esophagus',
            'GutEP',
            'Heart',
            'Kidney',
            'LI',
            'Liver',
            'Lung',
            'MoDC',
            'MuscSat',
            'SI',
            'Skin',
            'Stomach',
            'WAT'
        ]:
            p = inout.get_internal_path(
                'datasets/general/sequencing/edgeR_TMM_normalized_bytissue_171113/combatch_corrected/{}.csv'.format(t))
            agg.append(pd.read_csv(p).rename(
                columns={'Unnamed: 0': 'gene_ensembl'}).set_index('gene_ensembl'))

        df_counts = pd.concat(agg, axis=1)
        df_counts.index.name = 'gene_ensembl'

        p = inout.get_internal_path(
            'datasets/general/sequencing/aging_map_qc_171031/df_meta.csv')
        df_meta = pd.read_csv(p).rename(
            columns={'Unnamed: 0': 'sample_name'}).set_index(
            'sample_name')

        p = inout.get_internal_path(
            'datasets/general/sequencing/aging_map_qc_171031/df_genes.csv')
        df_genes = pd.read_csv(p).set_index('gene_ensembl')

        # f = np.array([x in df_counts.columns for x in df_meta.index])
        # if f.sum() == df_counts.shape[1]:
        df_meta = df_meta.loc[df_counts.columns, :]

        p = inout.get_internal_path(
            'datasets/general/complementary/transcript_info_171103.csv')
        df_transcript_info = pd.read_csv(p)
        df_transcript_info.loc[:, 'present_in_genbank_info'] = True

        df_genes_extended = pd.merge(
            df_genes.reset_index(),
            df_transcript_info,
            left_on='gene_ensembl',
            right_on='gene_ensembl',
            how='left')
        df_genes_extended['present_in_genbank_info'] = df_genes_extended[
            'present_in_genbank_info'].fillna(False)
        df_genes = df_genes_extended

    elif dataset_name == 'TMM_normalized_bytissue_171113_171120_combatch':

        agg = []
        for t in [
            'Adrenal',
            'AM',
            'AT2',
            'BAT',
            'Blood',
            'Brain',
            'Cerebellum',
            'Esophagus',
            'GutEP',
            'Heart',
            'Kidney',
            'LI',
            'Liver',
            'Lung',
            'MoDC',
            'MuscSat',
            'SI',
            'Skin',
            'Stomach',
            'WAT'
        ]:
            p = inout.get_internal_path(
                'datasets/general/sequencing/edgeR_TMM_normalized_bytissue_171113/combatch_corrected/{}.csv'.format(t))
            agg.append(pd.read_csv(p).rename(
                columns={'Unnamed: 0': 'gene_ensembl'}).set_index('gene_ensembl'))

        for pfu_type in ['10_pfu', '150_pfu']:
            for cell_type in [
                'AM',
                'AT2',
                'Blood',
                'Lung',
                'MoDC',
            ]:
                p = inout.get_internal_path(
                    'datasets/general/sequencing/edgeR_TMM_normalized_bytissue_Flu_171120/combatch_corrected/{}/{}.csv'.format(
                        pfu_type, cell_type)
                )
                agg.append(pd.read_csv(p).rename(
                    columns={'Unnamed: 0': 'gene_ensembl'}).set_index('gene_ensembl'))

        df_counts = pd.concat(agg, axis=1)
        df_counts.index.name = 'gene_ensembl'

        p = inout.get_internal_path(
            'datasets/general/sequencing/aging_map_qc_171031/df_meta.csv')
        df_meta = pd.read_csv(p).rename(
            columns={'Unnamed: 0': 'sample_name'}).set_index(
            'sample_name')

        p = inout.get_internal_path(
            'datasets/general/sequencing/aging_map_qc_171031/df_genes.csv')
        df_genes = pd.read_csv(p).set_index('gene_ensembl')

        # f = np.array([x in df_counts.columns for x in df_meta.index])
        # if f.sum() == df_counts.shape[1]:
        df_meta = df_meta.loc[df_counts.columns, :]

        p = inout.get_internal_path(
            'datasets/general/complementary/transcript_info_171103.csv')
        df_transcript_info = pd.read_csv(p)
        df_transcript_info.loc[:, 'present_in_genbank_info'] = True

        df_genes_extended = pd.merge(
            df_genes.reset_index(),
            df_transcript_info,
            left_on='gene_ensembl',
            right_on='gene_ensembl',
            how='left')
        df_genes_extended['present_in_genbank_info'] = df_genes_extended[
            'present_in_genbank_info'].fillna(False)
        df_genes = df_genes_extended

    else:
        raise ValueError(
            'Did not specify input')

    if unambiguous_to_entrez:
        df_genes = df_genes[df_genes['genes_with_same_ensembl'] == 1]
        df_counts = df_counts[df_counts.index.isin(df_genes['gene_ensembl'])]

    if as_entrez:
        df_counts = pd.merge(
            df_counts.reset_index(),
            df_genes[['gene_ensembl', 'gene_ncbi']]).set_index(
            'gene_ncbi', verify_integrity=True).drop('gene_ensembl', axis=1)
        df_counts = df_counts.sort_index()

    return df_counts, df_meta, df_genes


def load_aging_map():
    """
    Loads aging map counts and meta-info

    Output:
        df_data     count table
        df_meta     meta-info about specimen
        df_gene     information about the used gene coordinates

    """

    # main definitons of meta-info
    run_main_folder = inout.get_internal_path(
        'datasets/general/sequencing/sync_from_quest_171020/mock_questmon2')
    ref_runs = aging_standardizer.reference_experimental_sets('aging_map')
    ref_results = aging_standardizer.reference_result_folder('171020')
    run_def = run_definitions()
    mouse_information = aging_standardizer.mouse_information()

    # define speciment that don't correspond to main biological samples
    non_samples = ['UMRNAERCC1', 'UMRNAERCC2', 'undefined', 'Control_skin']

    agg = []
    agg_meta = []
    agg_gene = []
    for run in ref_runs:
        p = os.path.join(
            run_main_folder,
            run,
            ref_results[run],
            '05_quantification',
            'htseq.all.counts.txt')
        df, df_gene = _raw_count_table_from_file(p)
        agg.append(df)

        mapping_scheme = run_def.loc[run, 'meta_pattern_name']

        anticipated_columns = [
            'age',
            'age_unit',
            'mouse_id',
            'pfu',
            'replicate_id',
            'tissue',
            'is_specimen']
        df_meta = pd.DataFrame(index=df.columns, columns=anticipated_columns)

        for s in df_meta.index:
            if s not in non_samples:
                df_meta.loc[s, 'is_specimen'] = True
                for k, v in aging_standardizer.sample_name_2_meta(s, mapping_scheme).items():
                    df_meta.loc[s, k] = v
            else:
                df_meta.loc[s, 'is_specimen'] = False

        df_meta.loc[:, 'run_name'] = run
        agg_meta.append(df_meta)
        agg_gene.append(df_gene)

    # Double check that all gene definitions are the same
    df_gene = agg_gene[0]
    for g in agg_gene:
        same = all(df_gene == g)
        if not same:
            raise ValueError('At least one gene detail differs.')

    df_data = pd.concat(agg, axis=1)
    df_meta = pd.concat(agg_meta, axis=0)

    df_meta['mouse_id'] = df_meta['mouse_id'].astype(float)
    mouse_information['mouse_id'] = mouse_information['mouse_id'].astype(float)

    extended_meta = pd.merge(
        df_meta.reset_index(),
        mouse_information,
        left_on='mouse_id',
        right_on='mouse_id',
        how='left')

    # check for potential mistakes
    f = extended_meta['pfu_x'] == extended_meta['pfu_y']
    bb = extended_meta[~f]
    if any(~bb['sample'].isin(non_samples)):
        raise ValueError('meta_doesnot fit')

    f = (extended_meta['age_unit'] == 'M') & (
        extended_meta['age_months'] != extended_meta['age'])
    if extended_meta[f].shape[0] > 0:
        raise ValueError('Age mapped wrongly')

    f = (extended_meta['age_unit'] == 'D') & (extended_meta['age'] == 14)
    if not all(extended_meta.loc[f, 'age_months'] == 0.5):
        raise ValueError('age wrong')

    # define meta-information to use
    mouse_info_to_use = [
        'sample',
        'age',
        'age_unit',
        'mouse_id',
        'pfu_x',
        'replicate_id',
        'tissue',
        'is_specimen',
        'run_name',
        'harvest_date',
        'experimental_batch',
        'clotted',
        'died_during_intubation',
        'tumor'
    ]

    extended_meta = extended_meta.loc[:, mouse_info_to_use].rename(
        columns={'pfu_x': 'pfu'})

    df_meta = extended_meta.set_index('sample')

    # Safety check, that counts and meta have use same positions
    if not all(df_data.columns == df_meta.index):
        raise ValueError('some sorting mistake')

    # Cross check with Sophia's data for conistency, and for finding samples passing quality control
    p = inout.get_internal_path(
        'datasets/other/direct_communication_from_sophia/171020_expression/all_expression_raw.csv')
    df_sophia_counts = pd.read_csv(p).set_index('Unnamed: 0')

    p = inout.get_internal_path(
        'datasets/other/direct_communication_from_sophia/171020_expression/all_metadata.csv')
    df_sophia_meta = pd.read_csv(p)

    df_sophia_meta = df_sophia_meta[df_sophia_meta['Exp'] == 'S']
    u = df_sophia_meta['Unnamed: 0']

    df_sophia_counts = df_sophia_counts.loc[:, u]

    a = df_data.loc[
        df_sophia_counts.index,
        u
    ].fillna(0)

    m = (df_sophia_counts != a).sum()
    if any(m[m > 0]):
        raise ValueError("Inconsistent with Sophia")
    else:
        df_meta.loc[:, 'passed_qc'] = df_meta.index.isin(
            df_sophia_meta['Unnamed: 0'])

    # additions (Oct27th):
    # - remove undefined samples
    # - add information on i7 and i5 primer

    f = np.array(df_meta.index == 'undefined')
    df_data = df_data.loc[:, ~f]
    df_meta = df_meta.loc[~f, :]

    agg = []
    for run in df_meta['run_name'].unique():

        p = inout.get_internal_path(
            (
                'datasets/general/sequencing/'
                'sync_from_quest_171020/sample_sheets/'
                'ceto_sample_sheets/'
                'SampleSheet_{}.csv'.format(run)))

        df = pd.read_csv(p, header=19)
        dff = df[['Sample_Name', 'I7_Index_ID', 'I5_Index_ID']].copy()
        dff.loc[:, 'run_name'] = run
        agg.append(dff)
    df = pd.concat(agg, axis=0)

    f = df['Sample_Name'] == 'undefined'
    df = df.loc[~f, :]

    f = df[['Sample_Name', 'run_name']].duplicated()
    if any(f):
        df[f]
        raise EnvironmentError('At least one sample duplicated.')

    df_meta = df_meta.reset_index()
    df_meta = pd.merge(
        df_meta,
        df.rename(columns={
            'Sample_Name': 'sample'}),
        left_on=['sample', 'run_name'],
        right_on=['sample', 'run_name'],
        how='left')
    df_meta = df_meta.set_index('sample')

    # additions (Oct 31st):
    # - replace age format
    # - enforce unambiguous sample IDs (append run name)

    f = (df_meta['age'] == 14) & (df_meta['age_unit'] == 'D')
    df_meta.loc[f, 'age'] = 0.5
    df_meta.loc[f, 'age_unit'] = 'M'

    df_meta.loc[:, 'run_id'] = df_meta.loc[
        :, 'run_name'].str.extract(
        '^[0-9]*_[A-Z0-9]*_([0-9]*)_[A-Z]*',
        expand=False)

    if all(df_meta.index == df_data.columns):
        labels = df_meta.index.values
    else:
        raise ValueError(
            'Labels for metainfo do not match samples')

    df_meta.loc[:, 'orig_sample_name'] = df_meta.index.values
    f = df_meta.index.duplicated(keep=False)
    labels[f] = df_meta.loc[f, 'orig_sample_name'] + \
        '_r' + df_meta.loc[f, 'run_id']
    df_data.columns = labels
    df_meta.index = labels

    return df_data, df_meta, df_gene


def run_definitions():

    p_project_definitions = [
        inout.get_internal_path(
            (
                'datasets/general/sequencing/'
                'sync_from_quest_171020/sample_sheets/170719/'
                'experimental_annotation.xlsx')),
        inout.get_internal_path(
            (
                'datasets/general/sequencing/'
                'sync_from_quest_171020/sample_sheets/'
                'sample_sheets_171004_catchup/'
                'experimental_annotation.xlsx'))
    ]

    agg = []
    for p in p_project_definitions:
        df = pd.read_excel(p, 'runs')
        agg.append(df)

    df = pd.concat(agg, axis=0)
    df = df.set_index('run_name', verify_integrity=True)
    df = df.sort_index()

    return df


def _raw_count_table_from_file(p):
    """
    Loads a raw count table

    Output:
        df      counts
        df_gene non-count information on genes
    """
    df = pd.read_table(p)
    df.index.name = 'gene_ensembl'
    columns_with_gene_annotation = ['Chr', 'Start', 'End', 'Length', 'Strand']

    df_gene = df.loc[:, columns_with_gene_annotation]

    df = df.drop(columns_with_gene_annotation, axis=1)
    df.columns.name = 'sample'
    return df, df_gene
