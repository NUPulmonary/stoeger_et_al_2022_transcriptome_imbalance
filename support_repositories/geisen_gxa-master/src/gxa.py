import glob
import os
import re

import numpy as np
import pandas as pd


import statistics

from shared import get_path, comma_separated_string_to_list, flatten
from shared import get_gene_ncbi_and_gene_ensembl_with_taxon
from shared import ensure_presence_of_directory

from os import listdir
from os.path import isdir


def export_as_table(df, name, split_by_taxon_ncbi=False):

    p = '/Users/tstoeger/Dropbox/Work/resources/workspace/gxa_190622'

    # if split_by_taxon_ncbi:
    #     for taxon in df['taxon_ncbi'].unique():
    #         dff = df[df['taxon_ncbi'] == taxon]
    #         p_out = os.path.join(
    #             p, '{}_{}.csv.gz'.format(
    #                 name,
    #                 int(taxon)
    #             )
    #         )
    #         ensure_presence_of_directory(p_out)
    #         dff.to_csv(p_out)
    # else:
    p_out = os.path.join(
        p, '{}.csv.gz'.format(
            name,
        )
    )
    ensure_presence_of_directory(p_out)
    df.to_csv(p_out, compression='gzip', index=False)


def parse_all_idf_files():
    """
    Parses all IDF files. These files appear to be proprietary
    format of EBI-GXA. They contain meta-information about
    experiments, which is not biological (e.g.: email of authors,
    pubmed id of publication)

    Note that, besides parsing, no information is changed
    (e.g.: original empty values and vocabulary is maintained)

    Only considers records with at least two elements per record
    in the format of IDF e.g.: qualifer /tab something but not
    solely qualifer

    Output:
        stacked table containing:
        experiment       experimental ID
        qualifier    e.g.: "author first name" or "pubmed"
        value        values for a given qualifier
        position_in_record  starting with 1, the position if multiple
                            values for a rcecord within a given experiment
    """

    agg_experiments = []

    mypath = get_path('gtx_atlas')
    paths = [f for f in listdir(mypath) if isdir(os.path.join(mypath, f))]
    paths = [os.path.join(get_path('gtx_atlas'), x) for x in paths]

    for pe in paths:
        source = os.path.split(pe)[1]
        fname = os.path.join(mypath, source, '{}.idf.txt'.format(source))

        if os.path.exists(fname):

            with open(fname, encoding='latin-1') as f:
                content = f.readlines()
            content = [x.strip() for x in content]

            agg = []
            for line in content:

                if line is not '':
                    elements = line.split('\t')
                    amount_of_elements = len(elements)

                    if amount_of_elements > 1:
                        to_process = range(1, amount_of_elements)
                        for j in to_process:
                            agg.append(
                                {
                                    'qualifier': elements[0],
                                    'value': elements[j],
                                    'position_in_record': j
                                })

            df = pd.DataFrame(agg)
            df.loc[:, 'experiment'] = source
            agg_experiments.append(df)

    d = pd.concat(agg_experiments).sort_values(
        ['experiment', 'qualifier', 'position_in_record']).reset_index(
            drop=True)[
        ['experiment', 'qualifier', 'value', 'position_in_record']
    ]

    return d


def get_all_de_comparisons():
    """
    Get results of all differential expressions - on per-gene level,
    per-GO-term, and per-Interpro-term -- all as precomputed by
    EBI-GXA; During processing a custom comparison key is generated
    for every differential expression (multiple keys if experiment
    has multiple sub-experiments, e.g.: conditions)

    Output:
        df_genes      (contains custom comparison_key)
        df_go         (contains custom comparison_key)
        df_interpro   (contains custom comparison_key)
        df_contrasts  links custom comparison key to contrast details

    """

    # --- Settings --- #
    require_canonical_id = True   # for GO and Interpro. Some carry other names

    # --- Processing  --- #

    agg_genes, agg_go, agg_interpro, agg_contrasts = _parse_all_differential_expressions()

    # organize results of differential gene expression analysis
    df_genes = pd.concat(agg_genes, axis=0, ignore_index=True).rename(
        columns={'Gene ID': 'gene_ensembl'})
    df_genes = df_genes.loc[:, [
        'gene_ensembl',
        'log2foldchange',
        'p-value',
        'comparison_key']]

    gi = get_gene_ncbi_and_gene_ensembl_with_taxon(
        taxon_id='all',
        discard_ambiguous=True)
    df_genes = pd.merge(df_genes, gi, how='left')

    df_genes_ignoring_non_defined_taxa = df_genes.dropna(subset=['gene_ncbi'])

    taxa_per_comparison = df_genes_ignoring_non_defined_taxa[
        ['comparison_key', 'taxon_ncbi']].drop_duplicates().groupby(
            'comparison_key').size()

    if taxa_per_comparison.max() > 1:
        raise EnvironmentError("""
        Some comparisons map to multiple taxa. Handing of these conditions
        has not been considered and/or implemented. Neither has it been
        necessary for earlier versions of GXA
        """)
    else:
        # keep dominant taxon for meta (and/or data storage)
        comparison_2_taxon = df_genes_ignoring_non_defined_taxa[
            ['comparison_key', 'taxon_ncbi']
        ].drop_duplicates()

    df_genes = df_genes.drop('taxon_ncbi', axis=1)

    df_genes = df_genes.dropna(subset=['gene_ncbi'])[
        ['comparison_key', 'gene_ncbi', 'log2foldchange', 'p-value'
         ]]  # do not consider gene ensembl

    df_go = pd.concat(agg_go, axis=0, ignore_index=True)
    df_go = df_go.loc[:, [
        'Term',
        'p adj (non-dir.)',
        'effect.size',
        'comparison_key']].rename(columns={
            'Term': 'GO_ID'
        })
    if require_canonical_id:
        f = df_go['GO_ID'].str.contains('^GO:')
        df_go = df_go.loc[f, :]

    df_interpro = pd.concat(agg_interpro, axis=0, ignore_index=True)
    df_interpro = df_interpro.loc[:, [
        'Term',
        'p adj (non-dir.)',
        'effect.size',
        'comparison_key']].rename(columns={
            'Term': 'interpro_id'
        })
    if require_canonical_id:
        f = df_interpro['interpro_id'].str.contains('^IPR[0-9]')
        df_interpro = df_interpro.loc[f, :]

    del agg_interpro

    df_contrasts = pd.DataFrame(
        agg_contrasts)
    df_contrasts = df_contrasts.loc[:, [
        'comparison_key',
        'experiment',
        'analysis',
        'comparison']]

    df_contrasts = pd.merge(
        df_contrasts, comparison_2_taxon, how='left')

    # only consider contrasts with taxon
    df_contrasts = df_contrasts.dropna(subset=['taxon_ncbi'])

    return df_genes, df_go, df_interpro, df_contrasts


def get_all_absolute_expressions(unit):
    """


    """

    all_expressions = _list_all_absolute_expressions(unit)

    agg_expr = []
    agg_conditions = []

    for pe in all_expressions:

        df, exists, labeller = _parse_single_absolute_expression(pe, unit)

        agg_expr.append(df)
        agg_conditions.append(labeller)

    group_to_annotation = pd.concat(agg_conditions)
    df_genes = pd.concat(agg_expr)

    gi = get_gene_ncbi_and_gene_ensembl_with_taxon(
        taxon_id='all',
        discard_ambiguous=True)
    df_genes = pd.merge(df_genes, gi, how='left')

    df_genes_ignoring_non_defined_taxa = df_genes.dropna(subset=['gene_ncbi'])

    group_to_annotation = group_to_annotation.drop_duplicates(
        ['assay_id', 'assay_label'], keep=False
    ).reset_index(drop=True).rename_axis('assay_key', axis=0).reset_index()

    group_to_annotation = group_to_annotation.rename(columns={
        'assay_id': 'condition'

    })

    df_genes_ignoring_non_defined_taxa = pd.merge(
        df_genes_ignoring_non_defined_taxa,
        group_to_annotation,
        on=['experiment', 'condition']
    )

    df_genes_ignoring_non_defined_taxa = df_genes_ignoring_non_defined_taxa[
        [unit, 'gene_ncbi', 'taxon_ncbi', 'assay_key']]

    taxa_per_comparison = df_genes_ignoring_non_defined_taxa[
        ['assay_key', 'taxon_ncbi']].drop_duplicates().groupby(
        'assay_key').size()

    if taxa_per_comparison.max() > 1:
        raise EnvironmentError("""
        Some comparisons map to multiple taxa. Handing of these conditions
        has not been considered and/or implemented. Neither has it been
        necessary for earlier versions of GXA
        """)
    else:
        # keep dominant taxon for meta (and/or data storage)
        comparison_2_taxon = df_genes_ignoring_non_defined_taxa[
            ['assay_key', 'taxon_ncbi']
        ].drop_duplicates()

    df_genes_ignoring_non_defined_taxa = df_genes_ignoring_non_defined_taxa.drop(
        'taxon_ncbi', 1)

    df_assays = pd.merge(
        group_to_annotation, comparison_2_taxon, how='left')

    df_assays = df_assays[df_assays['assay_key'].isin(
        df_genes_ignoring_non_defined_taxa['assay_key'])]

    return df_genes_ignoring_non_defined_taxa, df_assays


def _parse_all_differential_expressions():

    mypath = get_path('gtx_atlas')
    paths = [f for f in listdir(mypath) if isdir(os.path.join(mypath, f))]
    paths = [os.path.join(get_path('gtx_atlas'), x) for x in paths]

    agg_genes = []
    current_comparison_key = 0
    agg_comparison_key = []
    agg_go = []
    agg_interpro = []

    for pe in paths:
        experiment = os.path.split(pe)[1]
        pa = glob.glob(os.path.join(
            pe, '{}[-_]*analytics.tsv'.format(experiment)))
        if len(pa) > 0:
            for paa in pa:   # sometimes there are several analyses

                df_genes = pd.read_table(paa, low_memory=False)

                _, current_analysis = os.path.split(paa)

                helper = pd.DataFrame(
                    data=df_genes.columns,
                    columns=['name'])
                comparisons = helper[helper['name'].str.endswith(
                    'log2foldchange')].loc[
                        :, 'name'].str.extract(
                            '^(.*)\.', expand=False).values

                for co in comparisons:

                    agg_comparison_key.append({
                        'comparison_key': current_comparison_key,
                        'experiment': experiment,
                        'comparison': co,
                        'analysis': current_analysis
                    })

                    usecols = [
                        'Gene ID',
                        '{}.log2foldchange'.format(co),
                        '{}.p-value'.format(co),
                    ]
                    dff = df_genes.loc[:, usecols]
                    dff = dff.rename(columns={
                        '{}.log2foldchange'.format(co): 'log2foldchange',
                        '{}.p-value'.format(co): 'p-value',
                    })

                    dff = dff.sort_values(
                        'p-value', ascending=True, na_position='last')

                    dff.loc[:, 'comparison_key'] = current_comparison_key
                    agg_genes.append(dff)

                    p_go = os.path.join(
                        pe, '{}.{}.go.gsea.tsv'.format(
                            experiment,
                            co))
                    if os.path.exists(p_go):
                        contains_data, df = _read_go_or_interpro(p_go)
                        if contains_data:
                            df.loc[
                                :, 'comparison_key'] = current_comparison_key
                            agg_go.append(df)

                    p_interpro = os.path.join(
                        pe, '{}.{}.interpro.gsea.tsv'.format(
                            experiment,
                            co))
                    if os.path.exists(p_interpro):
                        contains_data, df = _read_go_or_interpro(p_interpro)
                        if contains_data:
                            df.loc[
                                :, 'comparison_key'] = current_comparison_key
                            agg_interpro.append(df)

                    current_comparison_key += 1

    return agg_genes, agg_go, agg_interpro, agg_comparison_key


def _parse_single_absolute_expression(pe, unit):

    # mypath = get_path('gtx_atlas')
    # paths = [f for f in listdir(mypath) if isdir(os.path.join(mypath, f))]
    # paths = [os.path.join(get_path('gtx_atlas'), x) for x in paths]

    # agg_genes = []

    # counter = 0  # only for development: part of mechanism to limit search

    # for pe in paths:

        # if counter < 5:  # only for development: part of mechanism to limit search

    experiment = os.path.split(pe)[1]
    p_abundance = os.path.join(
        pe, '{}-{}.tsv'.format(experiment, unit))

    p_annotation = os.path.join(
        pe, '{}-configuration.xml'.format(experiment))

    exists = False
    df = []

    if (os.path.exists(p_abundance)) & (os.path.exists(p_annotation)):

        exists = True

        # Load group labels
        pattern = '^.*assay_group id=\"(.*)\" label=\"(.*)\">.*$'
        group_to_annotation = dict()
        with open(p_annotation) as configuration:
            for li in configuration:
                if '<assay_group id' in li:
                    catch = re.search(pattern, li)
                    group_key = catch.group(1)
                    label = catch.group(2)
                    group_to_annotation[group_key] = label

        with open(p_abundance) as f:
            first_line = f.readline().strip()

        agg = []
        genes = []

        a_line_had_been_read = False

        fh = open(p_abundance)
        for line in fh:

            if a_line_had_been_read:
                line = line.strip()
                split_line = line.split('\t')
                genes.append(split_line[0])

                elements = []
                for element in split_line[2:]:
                    within_group = [float(x) for x in element.split(',')]
                    elements.append(statistics.median(within_group))

                agg.append(elements)

            a_line_had_been_read = True
        fh.close()

        df = pd.DataFrame(index=genes, data=agg)
        df.index.name = 'gene_ensembl'
        df.columns = first_line.split('\t')[2:]
        # df.columns = [group_to_annotation[x] for x in df.columns]

        df = df.stack().to_frame(unit).reset_index().rename(
            columns={'level_1': 'condition'})

        df.loc[:, 'experiment'] = experiment

        labeller = pd.Series(
            group_to_annotation).to_frame(
            'assay_label').rename_axis('assay_id', 0).reset_index()

        labeller.loc[:, 'experiment'] = experiment

    return df, exists, labeller


def _list_all_absolute_expressions(unit):

    mypath = get_path('gtx_atlas')
    paths = [f for f in listdir(mypath) if isdir(os.path.join(mypath, f))]
    paths = [os.path.join(get_path('gtx_atlas'), x) for x in paths]

    agg = []

    counter = 0  # only for development: part of mechanism to limit search

    for pe in paths:

        experiment = os.path.split(pe)[1]
        p_abundance = os.path.join(
            pe, '{}-{}.tsv'.format(experiment, unit))

        p_annotation = os.path.join(
            pe, '{}-configuration.xml'.format(experiment))

        if (os.path.exists(p_abundance)) & (os.path.exists(p_annotation)):

            counter = counter + 1  # only for development: part of mechanism to limit search

            agg.append(pe)

    return agg


def _read_go_or_interpro(p):
    """
    Loads columns ['Term', 'p adj (non-dir.)', 'effect.size']
    from GO and Interpro enrichment analyses - while also
    capturing the exception - that files contain no data,
    but have single line starting with "V" (note that this
    is observation of Thomas Stoeger, rather than part
    of documentation of EBI GXA)


    Output:
        file_contains_data   True if file contains enrichments
        df                   Results of enrichment analysis
                                as stored by EBI-GXA
    """

    file_contains_data = False
    df = []

    lines = 0   # count whether a) 0; b) 1; c) >=2 lines
    with open(p, encoding='latin-1') as text_file:
        for line in text_file:
            lines = lines + 1
            if lines == 2:
                break

    if lines == 0:
        raise ValueError(
            'Format of {} does not match any expectation.'.format(
                p
            ))
    elif lines == 1:
        expected_for_alternate = line.startswith('V')
        expected_for_common = 'Term' in line

        if expected_for_alternate & ~expected_for_common:
            file_contains_data = False
        else:
            raise ValueError(
                'Format of {} does not match any expectation.'.format(
                    p
                ))
    elif lines > 1:
        df = pd.read_table(p, usecols=[
            'Term', 'p adj (non-dir.)', 'effect.size'],
            low_memory=False)
        file_contains_data = True
    else:
        raise EnvironmentError(
            'Number of lines in {} does not match any expectation.'.format(
                p
            ))

    return file_contains_data, df


def parse_all_configuration_xmls():

    def _defensive_update(
            dictionary, key, is_in_contrast, value_to_be_stored):

        if is_in_contrast is False:
            raise ValueError(
                'Lies outside of contrast.')

        if key in dictionary.keys():
            raise ValueError(
                'Key has already been set.')
        else:
            dictionary[key] = value_to_be_stored

        return dictionary

    agg_experiments = []

    mypath = get_path('gtx_atlas')
    paths = [f for f in listdir(mypath) if isdir(os.path.join(mypath, f))]
    paths = [os.path.join(get_path('gtx_atlas'), x) for x in paths]

    to_monitor_for_unexpected = []

    for pe in paths:
        experiment = os.path.split(pe)[1]

        p = os.path.join(pe, '{0}-configuration.xml'.format(experiment))

        if os.path.exists(p):

            with open(p, 'r') as file:
                config_file = file.readlines()

            config_file = [x.strip() for x in config_file]

            is_in_contrast = False
            is_in_comment = False

            for line in config_file:

                if line.startswith('<!--'):
                    is_in_comment = True
                elif line.startswith('-->'):
                    is_in_comment = False

                if is_in_comment is False:

                    if line.startswith('<contrast id='):
                        is_in_contrast = True
                        contrast = dict()
                        contrast['experiment'] = experiment

                        potential_contrast_id = re.findall(
                            '<contrast id="(.*)">', line)
                        if len(potential_contrast_id) == 0:
                            raise ValueError(
                                'Found no value for contrast.')
                        elif len(potential_contrast_id) > 1:
                            raise ValueError(
                                'Found multiple values for contrast.')
                        else:
                            contrast['contrast_id'] = potential_contrast_id[0]

                    elif line.startswith("<name>"):
                        label = line.replace(
                            "<name>", '').replace(
                            "</name>", '').replace(
                            "'", '').replace(
                            ' vs ', ' VS '
                        )
                        contrast = _defensive_update(
                            dictionary=contrast,
                            key='name',
                            is_in_contrast=is_in_contrast,
                            value_to_be_stored=label)

                    elif line.startswith('<reference_assay_group>'):
                        label = line.replace(
                            '<reference_assay_group>', '').replace(
                            "</reference_assay_group>", '')
                        contrast = _defensive_update(
                            dictionary=contrast,
                            key='reference_assay_group',
                            is_in_contrast=is_in_contrast,
                            value_to_be_stored=label)

                    elif line.startswith('<test_assay_group>'):
                        label = line.replace(
                            '<test_assay_group>', '').replace(
                            "</test_assay_group>", '')
                        contrast = _defensive_update(
                            dictionary=contrast,
                            key='test_assay_group',
                            is_in_contrast=is_in_contrast,
                            value_to_be_stored=label)

                    elif line.startswith('</contrast'):
                        is_in_contrast = False
                        agg_experiments.append(contrast)

                    elif is_in_contrast:
                        to_monitor_for_unexpected.append(line)

    monitor = pd.Series(to_monitor_for_unexpected).to_frame('suspicious')

    permitted_skips = [
        '^<batch.*',
        '^</batch.*',
        '^<assay.*',
        '^</assay.*',
        '^<--',
        '-->$'
    ]

    monitor.loc[:, 'allowed_skip'] = False

    for pattern in permitted_skips:
        f = monitor.loc[:, 'suspicious'].str.contains(pattern)
        monitor.loc[:, 'allowed_skip'] = monitor.loc[:, 'allowed_skip'] | f

    if all(monitor['allowed_skip']) is False:
        raise ValueError(
            'Skipped line with unanticipated pattern.'
        )

    parsed_configuration_xmls = pd.DataFrame(agg_experiments)
    return parsed_configuration_xmls
