import glob
import os

import pandas as pd

from access_science_shared import inout


def dais(subset, allowed_dais='all', allowed_wos='all'):
    """
    Loads all disambiguated author data for Web Of Science,
    together with authorship information; Will iterate through
    batched data; Thus filtering might be used to reduce memory
    footprint (e.g.: if wos IDs are arleady known)

    Input:
        subset      wos_dais that should be loaded; options:
                        'gene-linked' and 'all'
        allowed_dais    list of dais that should be loaded
        allowed_wos     list of wos IDs that should be loaded

    Output:
        df_wos_dais

    """

    p = inout.get_path(
        'rbusa',
        'disambiguation/wos_dais')

    if subset == 'gene-linked':
        mask = os.path.join(p, 'wos_dais_gene_mapped_batch_*.csv.gz')
    elif subset == 'all':
        mask = os.path.join(p, 'wos_dais_all_batch*.csv.gz')
    else:
        raise EnvironmentError('subset not specified')

    # Web of Science IDs can come in two formats, one which is
    # leaderless, and another one which fills up to 15 digits
    # Check explicitly to account for different input formats
    target_amount_of_numbes = 15
    if allowed_wos is not 'all':
        m = [len(x) == target_amount_of_numbes for x in allowed_wos]
        if all(m):
            has_leader = False
        else:
            has_leader = True

    agg = []
    for fi in glob.glob(mask):
        df = pd.read_csv(fi)

        if allowed_dais is not 'all':
            f = df.loc[:, 'DAIS'].isin(allowed_dais)
            df = df.loc[f, :]

        if allowed_wos is not 'all':

            if has_leader is False:
                df.loc[:, 'WOS'] = df.loc[:, 'WOS'].astype(str).apply(
                    lambda x: x.zfill(target_amount_of_numbes))
                df.loc[:, 'WOS'] = df.loc[:, 'WOS'].astype(str)

            f = df.loc[:, 'WOS'].isin(allowed_wos)
            df = df.loc[f, :]

        agg.append(df)

    df = pd.concat(agg)
    df = df.rename(columns={'WOS': 'wos_id', 'DAIS': 'dais_id'})
    df.loc[:, 'dais_id'] = df.loc[:, 'dais_id'].astype(int)
    df.loc[:, 'wos_id'] = df.loc[:, 'wos_id'].astype(str)

    df.loc[:, 'wos_id'] = df.loc[:, 'wos_id'].apply(
        lambda x: x.zfill(target_amount_of_numbes))
    df.loc[:, 'wos_id'] = df.loc[:, 'wos_id'].astype(str)

    return df
