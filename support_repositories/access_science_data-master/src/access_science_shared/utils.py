import pandas as pd


def split_text_to_multiple_rows(df, column, separator):
    """
    Separates entries, that have been separted by separator
    within a given
    column of the dataframe df into separate rows

    Input:
        df          DataFrame
        column      String; Name of the column,
                        which contains records separated by separator
        separator   String: Regular Expression that
                        separates the records in column

    Output:
        df_stacked  DataFrame, where records in
                        column have been separated into individual rows

    """

    # Separate rows, that should be processed (to save overhead)
    f = df[column].str.contains(separator)
    dff_s = df.loc[f, :]    # _s  separate
    dff_ns = df.loc[~f, :]  # _ns not separate
    orig_column_order = df.columns  # track original order of columns

    # Separate records and place them in separate rows
    s = dff_s[column].str.split(separator).apply(pd.Series, 1).stack()
    s.index = s.index.droplevel(-1)
    s.name = column
    del dff_s[column]
    dff_s = dff_s.join(s)

    # Recreate DataFrame in original input format
    dff_s = dff_s.reindex(orig_column_order, axis=1)
    df = pd.concat([dff_s, dff_ns])
    df_stacked = df

    return df_stacked


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


def flatten(l):
    return [item for sublist in l for item in sublist]


def comma_separated_string_to_list(x):
    unique_counts = list(map(float, x.split(',')))
    return unique_counts
