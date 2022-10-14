import pandas as pd

from access_science_shared import inout


def names_and_targets():
    """
    Loads synonyms (including those with speeling mistakes) and
    target genes (NCBI gene for positive nubmers) (or code of
    inapplicability for 0 and lower)

    Output:
        df_names
        df_targets

    """

    p = inout.get_path(
        dataset='antibodies',
        extension='antibody_id_2_names.csv')
    df_names = pd.read_csv(p)

    df_names = df_names.rename(columns={
        'id': 'antibody_id',
        'name': 'antibody_name'})

    p = inout.get_path(
        dataset='antibodies',
        extension='antibody_id_2_target_gene.csv')
    df_targets = pd.read_csv(p)

    return df_names, df_targets
