
import os

import pandas as pd


from access_biology_data import meta


def get_path(dataset=None, extension=None):
    '''
    Returns subfolder within internal part of geisen.

    Input:
        dataset     str, name of dataset, e.g.: geisen
        extension   str, optional, subfolder
    Output:
        outpath     str, folder within internal part of geisen
    '''

    # is_found = False
    # for l in locations:    # upper paths will be preferred
    #     if not is_found:
    #         p = os.path.join(l, datasets[dataset])
    #         if os.path.exists(p):
    #             ref_folder = p
    #             is_found = True
    # if not is_found:
    #     raise EnvironmentError(
    #         'Could not find location of {} dataset at {}'.format(
    #             dataset, p))

    if dataset == 'gtx_atlas':
        ref_folder = '/Users/tstoeger/outside_time_machine/190318/atlas-latest-data'

    if extension is not None:
        extension = str.replace(extension, '\\', os.path.sep)
        extension = str.replace(extension, '/', os.path.sep)

        outpath = os.path.join(ref_folder, extension)
    else:
        outpath = ref_folder

    return outpath


def flatten(l):
    return [item for sublist in l for item in sublist]


def comma_separated_string_to_list(x):
    unique_counts = list(map(float, x.split(',')))
    return unique_counts


def stack_by_delimiter_in_column(df, column, delimiter):
    """
    Stacks dataframe according to delimiter in column

    Input:
        df          dataframe
        column      column with delimiter
        delimiter   delimiter (note: no regular expression)

    Output:
        stacked_df  stacked dataframe

    """

    df.loc[:, column] = df.loc[:, column].astype(str)

    orig_index_name = df.index.name
    orig_column_order = df.columns

    df.index.name = 'original_index_used_before_splitting'
    df = df.reset_index()
    df.index.name = 'helper_index'

    f = (df[column].str.contains(
        delimiter, regex=False)) | (df[column].isnull())
    df_no_delimiter = df[~f]
    df_with_delimiter = df[f]

    ser_with_delimiter = df.loc[:, column]

    agg_values = []
    agg_indices = []

    for i, v in ser_with_delimiter.iteritems():
        vi = v.split(delimiter)
        indices = [i] * len(vi)

        agg_values.append(vi)
        agg_indices.append(indices)

    agg_values = flatten(agg_values)
    agg_indices = flatten(agg_indices)

    g = pd.DataFrame(data={'helper_index': agg_indices, column: agg_values})

    df_with_delimiter = pd.merge(
        df_with_delimiter.drop(column, 1).reset_index(),
        g)

    joined = pd.concat([
        df_no_delimiter.reset_index(),
        df_with_delimiter],
        sort=True
    )

    joined = joined.sort_values(
        ['original_index_used_before_splitting', column])
    joined = joined.drop('helper_index', 1)
    joined = joined.set_index('original_index_used_before_splitting')
    joined.index.name = orig_index_name
    joined = joined.loc[:, orig_column_order]

    return joined


def get_gene_ncbi_and_gene_ensembl_with_taxon(
        taxon_id='all', discard_ambiguous=True):
    """
    Obtains a dataframe mapping ncbi_gene to ncbi_ensembl

    Input:
        taxon_id                'all' (default) or ncbi taxonomoy id
                                or list of ncbi taxonomy ids: gene
                                information which shall be considered
                                (for gene_ncbi and gene_ensembl)

        discard_ambiguous        optional, default: True
                                removes amgiguous ncbi_gene and/or
                                ncbi_ensembl

    Output:
        gene_ncbi_ensembl       Dataframe with gene_ncbi and gene_ensembl
                                and taxon_ncbi

    """

    gi = meta.gene_info(
        taxon_id=taxon_id,
        usecols=['taxon_ncbi', 'gene_ncbi', 'dbXrefs'])
    gi = stack_by_delimiter_in_column(gi, 'dbXrefs', '|')

    gi = gi[gi['dbXrefs'] != '-'].copy()
    f = gi['dbXrefs'].str.contains(':')
    if f.mean() == 1:
        gi['base'] = gi['dbXrefs'].str.extract('(.*):.*')
        gi['tail'] = gi['dbXrefs'].str.extract('.*:(.*)')
    else:
        raise EnvironmentError('database references by ncbi are ambiguous')

    # for some organisms ensembl inherits
    # database identifiers from others.
    # NCBI still sometimes lists them differently
    allowed_bases = [
        'Ensembl',
        'Araport',
        'TAIR',
        'FLYBASE',
        'WormBase'
    ]

    gi = gi.loc[
        gi['base'].isin(allowed_bases),
        ['taxon_ncbi', 'gene_ncbi', 'tail']].drop_duplicates(
    ).rename(columns={'tail': 'gene_ensembl'})

    counts_in_ncbi = gi['gene_ncbi'].value_counts()
    counts_in_ensembl = gi['gene_ensembl'].value_counts()

    if discard_ambiguous:
        gi = gi[
            (gi['gene_ncbi'].isin(counts_in_ncbi[counts_in_ncbi == 1].index)) &
            (gi['gene_ensembl'].isin(
                counts_in_ensembl[counts_in_ensembl == 1].index))
        ]

    gene_ncbi_ensembl = gi.reset_index(drop=True)  # format output

    return gene_ncbi_ensembl


def ensure_presence_of_directory(directory_path=None, ):
    '''
    Ensure that the directory of exists. Creates dictionary with cognate
    name in case that directory does not exist. Anticipates that files have
    extensions separated by '.' symbol (can not create directory with . in
    its name); If file does not have an extension separated by '.' a folder
    will with its filname will be created, a behavior that can be avoided
    by calling os.path.dirname prior this function.

    Input:
        directory_path      str; Name of a directory or the full path of a file
    '''
    if directory_path is None:
        raise ValueError('No input specfied for ensure_presence_of_directory')

    directory_path_n, ext = os.path.split(directory_path)

    if '.' in ext:
        directory_path = directory_path_n

    if not os.path.exists(directory_path):
        os.makedirs(directory_path)
