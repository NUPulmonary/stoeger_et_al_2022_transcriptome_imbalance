import pandas as pd


from access_biology_data import meta
from access_science_shared import inout


def gtex_all_tpm():
    """
    Will download the full GTEX data set;
    Will only consider uambiguous linkages between NCBI gene and Ensembl gene;

    Note that presently this function directly parses GTEX and is thus
    slow (~3-5 min), which could be made much faster, possibly seconds,
    through caching (e.g.: feather), if the function was frequently used

    Output:
        df_gtex     per gene TPM measurements; rows: genes; columns: samples
        df_meta     sample and donor meta information

    """

    # import GTEX dataset
    p = inout.get_path(
        'gtex', 'GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct')

    df_gtex = pd.read_table(   # This is kind of slow
        p,
        sep='\t',
        skiprows=2    # row 1: GTEX version; row 2: dimensions
    ).rename(
        columns={
            'Name': 'gene_ensembl.version',
            'Description': 'gene_symbol'}).drop(
                'gene_symbol', axis=1)

    df_gtex['gene_ensembl'] = df_gtex['gene_ensembl.version'].str.extract(
        '(.*)\.[0-9]*', expand=False)
    df_gtex = df_gtex.drop('gene_ensembl.version', axis=1)

    # get information on donors
    p = inout.get_path(
        'gtex', 'GTEx_v7_Annotations_SubjectPhenotypesDS.txt')

    df_sujects = pd.read_table(p)
    df_sujects = df_sujects.rename(columns={' SUBJID': 'SUBJID'})   # typo

    df_sujects['SEX'] = df_sujects['SEX'].replace({
        1: 'm',    # manually update according to accompanying GTEX excel file
        2: 'f'
    }
    )

    # get information on specimen
    p = inout.get_path(
        'gtex', 'GTEx_v7_Annotations_SampleAttributes.txt')
    df_sample_attributes = pd.read_table(p)

    # Add donor information to specimen information
    df_sample_attributes.loc[:, 'SUBJID'] = df_sample_attributes[
        'SAMPID'].str.extract(
        '^([^-]*-[^-]*).*',      # Donor appears encoded in part of sample name
            expand=False)

    df_meta = pd.merge(df_sample_attributes, df_sujects,
                       how='left', left_on='SUBJID', right_on='SUBJID')
    df_meta = df_meta.set_index('SAMPID', verify_integrity=True)

    # Update genex data to unambiguous entrez genes
    gene_info = meta.gene_info(9606, usecols=['gene_ncbi', 'dbXrefs'])

    gene_info['gene_ensembl'] = gene_info[
        'dbXrefs'].str.extract('.*Ensembl:(ENSG[0-9]*)', expand=False)

    occurences_of_ensg = gene_info['gene_ensembl'].value_counts().to_frame(
        'occurences_of_ensg').reset_index().rename(
        columns={'index': 'gene_ensembl'})

    gene_info = pd.merge(gene_info, occurences_of_ensg, how='left')
    gene_info = gene_info[['gene_ncbi', 'gene_ensembl', 'occurences_of_ensg']]

    unambiguous_ncbi_2_ensembl = gene_info.loc[
        gene_info['occurences_of_ensg'] == 1, [
            'gene_ncbi', 'gene_ensembl']]

    df_gtex = pd.merge(
        df_gtex,
        unambiguous_ncbi_2_ensembl,
        how='inner')

    df_gtex = df_gtex.drop('gene_ensembl', axis=1)
    df_gtex = df_gtex.set_index('gene_ncbi', verify_integrity=True)
    df_gtex = df_gtex.sort_index()

    return df_gtex, df_meta
