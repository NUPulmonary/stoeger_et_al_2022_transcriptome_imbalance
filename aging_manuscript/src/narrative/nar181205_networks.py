import sys

import numpy as np
import pandas as pd


from access_biology_data import relations

sys.path.append('./../src/')


def complemented_biogrid_bioplex_edges(target_taxon_id):
    """
    edges with all genes in biogrid and bioplex,
    mapped to target_taxon_id

    For interations reported in target_taxon_id only
    individual genes are considered, whereas for
    for other taxa all orthologs are considered
    (mapped to all their orthologs in target_taxon,
    and combined in case of several paralogs in
    other taxa)

    self interactions are removed, and only one
    edge per interaction is kept (undirectional,
    unweighted edges)

    Input:
        target_taxon_id   NCBI taxonomy ID

    Output:
        edges           dataframe containing interactions,
                            columns are called 'source' and
                            'target' for compatibility with
                            cytoscape

    """

    # Get Homologene
    homologene = relations.homologene()

    non_target_to_target = pd.merge(
        homologene.loc[
            homologene['taxon_ncbi'] == target_taxon_id, [
                'homologene_group', 'gene_ncbi']].rename(
            columns={'gene_ncbi': 'gene_ncbi_target_taxon'}),
        homologene.loc[
            homologene['taxon_ncbi'] != target_taxon_id, [
                'homologene_group', 'gene_ncbi']].rename(
            columns={'gene_ncbi': 'gene_ncbi'}))[[
                'gene_ncbi', 'gene_ncbi_target_taxon']]

    biogrid = relations.biogrid('all')

    biogrid = biogrid[
        ['Entrez Gene Interactor A', 'Entrez Gene Interactor B']
    ].rename(columns={
        'Entrez Gene Interactor A': 'source',
        'Entrez Gene Interactor B': 'target'})

    biogrid = pd.concat(
        [
            biogrid,
            biogrid.rename(columns={'source': 'target', 'target': 'source'})
        ],
        sort=True
    )

    biogrid = biogrid.replace('-', np.nan)
    biogrid = biogrid.dropna()
    biogrid = biogrid.astype(float)

    bioplex = relations.bioplex2()

    bioplex = pd.merge(
        bioplex,
        bioplex, left_on='bioplex2_id', right_on='bioplex2_id').rename(
        columns={'gene_ncbi_x': 'source', 'gene_ncbi_y': 'target'})

    bioplex = bioplex[['source', 'target']]

    edges = pd.concat([biogrid, bioplex])

    edges_foreign = pd.merge(
        edges, non_target_to_target, left_on='source', right_on='gene_ncbi')
    edges_foreign = edges_foreign[
        ['gene_ncbi_target_taxon', 'target']].rename(
        columns={'gene_ncbi_target_taxon': 'source'})

    edges_foreign = pd.merge(
        edges_foreign,
        non_target_to_target,
        left_on='target',
        right_on='gene_ncbi')
    edges_foreign = edges_foreign[
        ['gene_ncbi_target_taxon', 'source']].rename(
        columns={'gene_ncbi_target_taxon': 'target'})

    edges = pd.concat([edges, edges_foreign], sort=True)

    edges = pd.concat([
        edges,
        edges.rename(columns={'source': 'target', 'target': 'source'})],
        sort=True
    )

    edges = edges[
        edges['source'] < edges['target']
    ].drop_duplicates()

    return edges
