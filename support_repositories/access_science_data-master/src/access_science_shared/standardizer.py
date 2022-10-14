from access_biology_data import meta, relations
from access_literature_data import medline
from access_science_shared import utils


def filter_by_paper_kind(df, paper_kind):
    """
    Filters a dataframe for pubmed_id that also
    appear in paper_kind, which allows filtering of medline
    by potentially diverse options; prsently supported:
    'research'  --> research artciles of diverse kinds


    Input:
        df          dataframe, containing 'pubmed_id' in column
        paper_kind  input option for get_pubmed_ids_passing_medline_filter
                        e.g.: research; will filter according to
                        MedLine
    """

    if 'pubmed_id' not in df.columns:
        raise EnvironmentError(
            'Failed to find pubmed_id in provided input df')
    else:
        r = medline.select_medline_records(
            columns_sql='''
                medline.pubmed_id''',
            taxon_id=None,
            kind='research')
        r = r['pubmed_id'].values

        df = df[df.loc[:, 'pubmed_id'].isin(r)]
    return df


def reference_genes(taxon_id, ref_code):
    """
    Obtains a list of reference genes

    Input:
        taxon_id    int
        ref_code    str;  if it contains
                        l   -> at least one medline paper
                        o   -> official nomenclature require
                        p   -> protein-coding only

    Output:
        ref_genes   sorted list of gene identifiers
    """

    df = meta.gene_info(taxon_id)

    if df.shape[0] == 0:
        raise EnvironmentError("""
            Did not find gene info for taxon {}""".format(
            int(taxon_id)))

    if 'e' in ref_code:
        unambiguous = _get_gene_ncbi_and_gene_ensembl_with_taxon(
            taxon_id,
            discard_ambiguous=True)
        f = df.loc[:, 'gene_ncbi'].isin(unambiguous['gene_ncbi'])
        df = df.loc[f, :]

        if df.shape[0] == 0:
            raise EnvironmentError("""
                After filtering for genes unambiguous betewen ensembl
                and ncbi gene, no gene is left.""")

    if 'l' in ref_code:
        genes_in_medline = medline.gene2pubmed(taxon_id, ['gene_ncbi'])
        f = df.loc[:, 'gene_ncbi'].isin(genes_in_medline['gene_ncbi'])
        df = df.loc[f, :]

        if df.shape[0] == 0:
            raise EnvironmentError("""
                After filtering for genes with at least one paper,
                no gene is left.""")

    if 'o' in ref_code:  # official nomeclature
        f = df.loc[:, 'Nomenclature_status'] == 'O'
        df = df.loc[f, :]

        if df.shape[0] == 0:
            raise EnvironmentError("""
                After filtering for genes with official nomeclature,
                no gene is left.""")

    if 'p' in ref_code:  # protein-coding
        f = df.loc[:, 'type_of_gene'] == 'protein-coding'
        df = df.loc[f, :]

        if df.shape[0] == 0:
            raise EnvironmentError("""
                After filtering for protein-coding, no gene is
                left.""")

    if 'r' in ref_code:
        genes_in_medline = medline.gene2pubmed(
            taxon_id, ['pubmed_id', 'gene_ncbi'], paper_kind='research')
        f = df.loc[:, 'gene_ncbi'].isin(genes_in_medline['gene_ncbi'])
        df = df.loc[f, :]

        if df.shape[0] == 0:
            raise EnvironmentError("""
                After filtering for genes with at least one research paper,
                no gene is left.""")

    if 'h' in ref_code:
        hg = relations.homologene()
        f = df.loc[:, 'gene_ncbi'].isin(hg['gene_ncbi'])
        df = df.loc[f, :]

        if df.shape[0] == 0:
            raise EnvironmentError("""
                After filtering for genes in homologene,
                no gene is left.""")

    ref_genes = sorted(df.loc[:, 'gene_ncbi'].values)

    return ref_genes


def _get_gene_ncbi_and_gene_ensembl_with_taxon(
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
    gi = utils.stack_by_delimiter_in_column(gi, 'dbXrefs', '|')

    f = gi['dbXrefs'] != '-'
    gi = gi.loc[f, :].copy()
    if f.sum() == 0:
        raise EnvironmentError(
            'No dbXrefs provided')
    else:
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
