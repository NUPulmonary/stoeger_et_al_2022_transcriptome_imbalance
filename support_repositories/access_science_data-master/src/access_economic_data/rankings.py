import os
import pandas as pd

from access_science_shared import inout


def topuniversities():
    """
    Will load ranking of topuniviersities

    Output:
        df_top_universities

    """

    p = inout.get_path(
        'webpages',
        os.path.join(
            'topuniversities',
            'topuniversities_biology.txt'))

    with open(p) as fi:
        top_universities = fi.readlines()

    top_universities = top_universities[3:]

    top_universities = [x.strip() for x in top_universities]

    df_top_universities = pd.DataFrame(columns=['position', 'name'])
    df_top_universities.loc[:, 'position'] = top_universities[::2]
    df_top_universities.loc[:, 'name'] = top_universities[1::2]
    df_top_universities['name'] = df_top_universities['name'].str.replace(
        '\t', ' ')

    f = df_top_universities['name'].str.contains(' Logo')
    df_top_universities.loc[f, 'name'] = df_top_universities.loc[
        f, 'name'].str.extract(
        '(.*) Logo.*', expand=False)

    f = df_top_universities['name'].str.contains('\t')
    df_top_universities.loc[f, 'name'] = df_top_universities.loc[
        f, 'name'].str.extract(
        '(.*)\t', expand=False)

    df_top_universities['name'] = df_top_universities['name'].str.strip()

    return df_top_universities
