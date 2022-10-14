import pandas as pd

from access_science_shared import inout


def gtr(allowed_type):
    """
    Accesses the genetic testing registry from NIH.
    Will ignore tests agains non-canonical genes, e.g.:
    BCL-ABL fusion gne

        Input:
                allowed_type   str or list with 'Clinical' and / or Research

    """

    if isinstance(allowed_type, str):
        allowed_type = [allowed_type]

    allowed_type = [x.capitalize() for x in allowed_type]

    p = inout.get_path(
        dataset='geisen',
        extension='nih/gtr/test_condition_gene.txt')
    df = pd.read_table(p)

    f = df['test_type'].isin(allowed_type)
    df = df.loc[f, :]

    # Ignore non-canonical genes, e.g.: against BCR-ABL fusion gene
    f = df['gene_or_SNOMED_CT_ID'] < 0
    accession_not_gene_specific = df.loc[f, '#accession_version'].values

    f = df['#accession_version'].isin(accession_not_gene_specific)
    df = df.loc[~f, :]

    return df
