import pandas as pd

from access_biology_data import meta
from access_science_shared import inout

from cmappython3 import parse


def _get_path(dataset):

    paths = {
        'conditions': inout.get_path(
            'lincs',
            'GSE92742/GSE92742_Broad_LINCS_sig_info.txt'),
        'genes': inout.get_path(
            'lincs',
            'GSE92742/GSE92742_Broad_LINCS_gene_info.txt'),
        'data': inout.get_path(
            'lincs',
                'GSE92742/GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx')
    }
    return paths[dataset]


def get_all_condition_meta():
    """
    loads LINCS meta on experimental conditions
    """
    sig_info = pd.read_csv(
        _get_path('conditions'),
        sep='\t',
        low_memory=False)
    return sig_info


def get_all_gene_meta():
    """
    loads LINCS meta on  genes
    """
    lincs_gene_info = pd.read_csv(
        _get_path('genes'),
        sep="\t")
    return lincs_gene_info


def load_gene_perturbations(genes_exp, genes_pert, pert_type='all'):
    """
    Loads LINCS experiemntal data for gene specific perurbations

    Input:
        genes_exp   list, genes whose expression shall be returned
                        (gene_ncbi identifiers)
        genes_pert  list, genes which should be considered as perturbants
                        (gene_ncbi identifiers)
        pert_type   str, or list, of perturbation that should be
                        considered. Alternatively 'all' to load
                        everything. Recommended options:
                            trt_sh.cgs  consensus signature of knockdown
                            trt_oe      overexpression
                        Full list of options: see documentation of
                            LINCS
                            clue.io/connectopedia/perturbagen_types_and_controls

    Output:
        data        df, rows: conditions, columns: genes whose
                        expression was monitored
        meta_conditions df, experimental conditions
        meta_genes      df, LINCS specific meta on genes

    """

    sig_info = get_all_condition_meta()
    column_order = sig_info.columns

    gi = meta.gene_info(9606, ['gene_ncbi', 'symbol_ncbi'])
    gi = gi.drop_duplicates('symbol_ncbi', keep=False)
    mapper = gi[gi['gene_ncbi'].isin(genes_pert)]

    if pert_type is not 'all':
        if isinstance(pert_type, str):
            pert_type = list(pert_type)
        sig_info = sig_info[sig_info['pert_type'].isin(pert_type)]

    sig_info = pd.merge(
        sig_info,
        mapper.rename(columns={'symbol_ncbi': 'pert_iname'})).drop(
        'pert_iname', 1).rename(columns={'gene_ncbi': 'pert_iname'})

    sig_info = sig_info.loc[:, column_order]

    lincs_gene_info = get_all_gene_meta()
    lincs_gene_info = lincs_gene_info[lincs_gene_info['pr_gene_id'].isin(
        genes_exp)]

    lincs_data = parse.parse(
        _get_path('data'),
        rid=[str(x) for x in lincs_gene_info['pr_gene_id']],
        cid=sig_info['sig_id'].values
    )

    lincs_data.index = [int(x) for x in lincs_data.index]
    lincs_data = lincs_data.sort_index().transpose().sort_index()
    lincs_data.columns.name = 'gene_ncbi'

    if not all(sig_info['sig_id'].isin(lincs_data.index)):
        raise ValueError('Internal error: not all conditions retreived.')

    if not all(lincs_data.index.isin(sig_info['sig_id'])):
        raise ValueError('Internal error: additional conditions retreived.')

    if not all(lincs_gene_info['pr_gene_id'].isin(lincs_data.columns)):
        raise ValueError('Internal error: not all genes retreived.')

    if not all(lincs_data.columns.isin(lincs_gene_info['pr_gene_id'])):
        raise ValueError('Internal error: additional genes retreived.')

    sig_info = sig_info.set_index(
        'sig_id', verify_integrity=True).sort_index().rename_axis('cid')

    lincs_gene_info = lincs_gene_info.sort_values(
        'pr_gene_id').set_index('pr_gene_id', verify_integrity=True)[
        ['pr_is_lm', 'pr_is_bing']
    ].rename_axis('gene_ncbi')

    data = lincs_data
    meta_conditions = sig_info
    meta_genes = lincs_gene_info

    return data, meta_conditions, meta_genes
