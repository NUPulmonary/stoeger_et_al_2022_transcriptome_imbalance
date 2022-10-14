import sys

import numpy as np
import pandas as pd

from access_biology_data import properties, relations
from access_science_shared import mapper, standardizer

sys.path.append('./../src/')
from access_aging_data import earlier_studies
from aging_tools import inout


# Aggregate loading functions

def agg_features_190402_with_mirnas(taxon_id):
    """
    Obtains several features on transcript- and gene-
    and protein-architecture, as well as the binding
    of transcription factors and miRNAs

    """

    miRNAs = get_miRNAs(taxon_id)

    biomart_transcripts = get_biomart_transcript_features(taxon_id)

    protein = get_protein(taxon_id)
    ncbi_rna = get_validated_rna(taxon_id)
    # exac = get_exac(taxon_id)
    chromosomes = get_chromosomes(taxon_id)
    bm = get_biomart_transcript_features(taxon_id)
    # hagr = get_hagr(taxon_id)
    tfs = get_transcription_factors(taxon_id)

    g_helper = pd.concat([
        protein,
        ncbi_rna,
        biomart_transcripts,
        #    exac,    # difficult to interpret
        chromosomes,
        bm
    ], axis=1, join='inner').dropna()

    g = pd.concat([
        g_helper,
        #    hagr,     # remove controversial
        tfs,
        miRNAs
    ], axis=1, join='outer').fillna(0)

    g = g[g.index.isin(g_helper.index)]
    g = g.astype(float)

    f = np.array([len(set(g.loc[:, x])) >= 2 for x in g.columns])
    g = g.loc[:, f].copy()

    return g

def agg_features_190402_no_mirnas(taxon_id):
    """
    Obtains several features on transcript- and gene-
    and protein-architecture, as well as the binding
    of transcription factors and miRNAs

    """

    #miRNAs = get_miRNAs(taxon_id)

    biomart_transcripts = get_biomart_transcript_features(taxon_id)

    protein = get_protein(taxon_id)
    ncbi_rna = get_validated_rna(taxon_id)
    # exac = get_exac(taxon_id)
    chromosomes = get_chromosomes(taxon_id)
    bm = get_biomart_transcript_features(taxon_id)
    # hagr = get_hagr(taxon_id)
    tfs = get_transcription_factors(taxon_id)

    g_helper = pd.concat([
        protein,
        ncbi_rna,
        biomart_transcripts,
        #    exac,    # difficult to interpret
        chromosomes,
        bm
    ], axis=1, join='inner').dropna()

    g = pd.concat([
        g_helper,
        #    hagr,     # remove controversial
        tfs
    ], axis=1, join='outer').fillna(0)

    g = g[g.index.isin(g_helper.index)]
    g = g.astype(float)

    f = np.array([len(set(g.loc[:, x])) >= 2 for x in g.columns])
    g = g.loc[:, f].copy()

    return g


def agg_features_190406_with_mirnas(taxon_id):
    """
    Obtains several features on transcript- and gene-
    and protein-architecture, as well as the binding
    of transcription factors and miRNAs

    """

    g_helper = pd.concat([
        get_validated_rna_no_codon_bias(taxon_id),
        get_biomart_transcript_features(taxon_id),
        get_chromosomes(taxon_id),
        get_gene(taxon_id)
    ], axis=1, join='inner').dropna()

    g = pd.concat([
        g_helper,
        get_transcription_factors(taxon_id),
        get_miRNAs(taxon_id)
    ], axis=1, join='outer').fillna(0)

    g = g[g.index.isin(g_helper.index)]
    g = g.astype(float)

    f = np.array([len(set(g.loc[:, x])) >= 2 for x in g.columns])
    g = g.loc[:, f].copy()

    return g


# Detail loading functions for specific datasets

def get_protein(taxon_id):
    r4 = properties.aminoacids_swissprot(taxon_id).set_index('gene_ncbi')
    r4.columns = [x.replace('Aminoacids_swissprot: ', '') for x in r4.columns]
    r4 = r4[
        [
            'amount_measured_amino_acids', 'aromatic', 'basic', 'charged',
            'gravy_ignoring_O_and_U', 'helix_affine',
            'hydrophobic', 'isoelectric_point', 'molecular_weight', 'polar',
            'polar_uncharged', 'sheet_affine', 'turn_affine']


    ]
    r4.columns = ['protein_{}'.format(x) for x in r4.columns]
    return r4


def get_gene(taxon_id):
    r25 = properties.genbank_gene(taxon_id).set_index('gene_ncbi')
    r25.columns = [x.replace('Genbank__gene: ', 'gene_') for x in r25.columns]
    return r25


def get_validated_rna_no_codon_bias(taxon_id):

    r1 = properties.genbank_validated_rna(taxon_id).set_index('gene_ncbi')
    r1.columns = [x.replace('Genbank_validated_RNA: ', '') for x in r1.columns]
    r1 = r1[[
        'cds_SumACGT',
        'full_SumACGT',
        'full_CG'
    ]]
    r1.columns = ['rna_{}'.format(x) for x in r1.columns]
    return r1


def get_validated_rna(taxon_id):

    r1 = properties.genbank_validated_rna(taxon_id).set_index('gene_ncbi')
    r1.columns = [x.replace('Genbank_validated_RNA: ', '') for x in r1.columns]
    r1 = r1[[
        'cds_SumACGT',
        'codon_bias_cdc',
        'codon_bias_novembre2002',
        'codon_bias_rcbs',
        'codon_bias_scuo',
        'codon_bias_sun2013',
        'codon_bias_wright1990',
        'full_SumACGT',
        'full_CG'
    ]]
    r1.columns = ['rna_{}'.format(x) for x in r1.columns]
    return r1


def get_hagr(taxon_id):
    r5 = earlier_studies.hagr_mapped_summary(
        taxon_id)[['gene_ncbi', 'influence']].drop_duplicates()
    r5.loc[:, 'value'] = 1
    r5 = r5.pivot(index='gene_ncbi', columns='influence',
                  values='value').fillna(0)
    return r5


def get_exac(taxon_id):
    r6 = relations.relate_to_homologene(
        taxon_id,
        properties.allelepool_lek_2016_aberration(9606),
        extend_column_names=False).set_index('gene_ncbi')

    r6.columns = [x.replace('Population variability Lek ', 'ExAC_')
                  for x in r6.columns]
    return r6


def get_chromosomes(taxon_id):
    r3 = properties.chromosomes(taxon_id).set_index('gene_ncbi')
    f = r3['chromosome'].isin(['X', 'Y', 'X|Y'])
    r3.loc[f, 'chromosome'] = -1
    f = r3['chromosome'].isin(['MT'])
    r3.loc[f, 'chromosome'] = -2
    f = r3['chromosome'].isin(['Un', '-', ])
    r3.loc[~f, :]
    return r3


def get_transcription_factors(taxon_id):

    if taxon_id == 9606:
        p = inout.get_internal_path(
            "dynamic/tstoeger/190401_map_transcription_factors/human_200_interval.csv")
    elif taxon_id == 10090:
        p = inout.get_internal_path(
            "dynamic/tstoeger/190401_map_transcription_factors/mouse_200_interval.csv")

    tf = pd.read_csv(p, usecols=['gene_ensembl', 'tf_title']).drop_duplicates()
    tf = mapper.gene_ensembl_2_gene_ncbi_unambiguously(
        tf, taxon_id).reset_index()
    tf.loc[:, 'present'] = True
    tf = tf.pivot(index='gene_ncbi', columns='tf_title',
                  values='present').fillna(False)
    tf.loc[:, 'different_tfs'] = tf.sum(1)
    tf.columns = ['tf_{}'.format(x) for x in tf.columns]

    return tf


def get_biomart_transcript_features(taxon_id):

    if taxon_id == 9606:
        p = inout.get_internal_path(
            'datasets/general/resources/biomart/190401_biomart/human_mart_export.txt')
    elif taxon_id == 10090:
        p = inout.get_internal_path(
            'datasets/general/resources/biomart/190401_biomart/mouse_mart_export.txt')

    df = pd.read_table(p)

    df = df[df['Chromosome/scaffold name'].isin(
        [str(x) for x in range(100)] + ['X', 'Y', 'MT']
    )].copy()

    df = df.rename(columns={
        'Gene stable ID': 'gene_ensembl',
        'Transcript stable ID': 'transcript_ensembl'
    })

    exons_per_transcript = df[[
        'gene_ensembl',
        'transcript_ensembl',
        'Exon region start (bp)']].drop_duplicates()

    d = exons_per_transcript.groupby(
        ['gene_ensembl',
         'transcript_ensembl']).size().to_frame('exons').reset_index()

    g = d[['gene_ensembl', 'exons']].groupby(
        'gene_ensembl').agg([max, min, np.median, len])

    g.columns = g.columns.droplevel()

    g = g.rename(columns={
        'max': 'exons_max_ensembl',
        'min': 'exons_min_ensembl',
        'median': 'exons_median_ensembl',
        'len': 'transcripts_ensembl',
    })

    h = mapper.gene_ensembl_2_gene_ncbi_unambiguously(g, taxon_id)
    return h


def get_miRNAs(taxon_id):
    mi = properties.mirdb(taxon_id)
    mi = mi.set_index('gene_ncbi')
    mi.columns = ['miRNA_{}'.format(x) for x in mi.columns]

    mi.loc[:, 'miRNA_total'] = mi.sum(1)
    mi = mi.astype(float)

    return mi
