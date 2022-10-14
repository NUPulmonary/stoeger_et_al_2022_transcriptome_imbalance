import numpy as np
import pandas as pd

from statsmodels.sandbox.stats.multicomp import multipletests

from scipy.stats import fisher_exact


def jaccard_index(set_x, set_y):
    # Compute jaccard index

    set_x = set(set_x)
    set_y = set(set_y)

    inner = set_x.intersection(set_y)
    outer = set_x.union(set_y)

    jaccard_index = len(inner) / len(outer)

    return jaccard_index


def gini(list_of_values):
    """
    Computes Gini Coefficient,
    ranging between 0 (no inequality) and 1 (absolute inequality)

    Input:
        list_of_values      list with values

    Output:
        gini_coefficient    float, Gini Coefficient,

    """

    # Modified after:
    # http://planspace.org/2013/06/21/
    # how-to-calculate-gini-coefficient-from-raw-data-in-python/

    sorted_list = sorted(list_of_values)
    height, area = 0, 0
    for value in sorted_list:
        height += value
        area += height - value / 2.

    length_of_list = len(list_of_values)

    if length_of_list > 0:
        fair_area = height * length_of_list / 2.
        result = (fair_area - area) / fair_area
    else:
        raise ValueError('Gini Coefficient is not defined for empty lists.')

    return result


def functional_ratio(
        significant_genes,
        background,
        painter,
        test='enrichment'):
    """
    Computes functional enrichments, and does so with directly interpretable
    Fisher test, that would also allow to test for depletion,

    Input
        significant_genes   list of entrez gene IDs found by user
        background       list of entrez gene IDs, background considered
                            by user ("significant" and "non-signficant"
                            genes)
        painter         df, with ['gene_ncbi', 'annotation']
        test            'enrichment' (test only for enrichment),
                        or 'difference' (test for enrichments or depletion)
    """

    if test == 'enrichment':
        alternative = 'greater'
    elif test == 'difference':
        alternative = 'two-sided'
    else:
        raise ValueError('Test only supports enrichment or difference')

    sign_genes = set(significant_genes).intersection(set(painter['gene_ncbi']))
    background = set(background).intersection(set(painter['gene_ncbi']))
    painter = painter.drop_duplicates()

    if len(sign_genes - set(background)) > 0:
        raise ValueError(
            (
                'The list of significant entries contains some'
                'entries absent from the background list.'
            )
        )
    painter = painter[painter['gene_ncbi'].isin(background)]
    non_sign_genes = set(background) - sign_genes

    minimally_required_candidates_in_annotation = 2
    c = painter['annotation'].value_counts()
    painter = painter[painter['annotation'].isin(
        c[c >= minimally_required_candidates_in_annotation].index)]
    painter = painter.set_index('annotation')

    annotations = painter.index.unique()

    collector = pd.DataFrame(
        index=annotations,
        columns=[
            'fold_enrichment',
            'pvalue',
            'bonferroni',
            'benjamini_hochberg',
            'in_significant',
            'in_background'],
    )

    if len(annotations) == 1:
        raise AssertionError('Must have multiple annotations')

    for annotation in annotations:

        genes_with_annotation = set(painter.loc[annotation, 'gene_ncbi'])
        genes_without_annotation = background - genes_with_annotation

        a = len(sign_genes.intersection(genes_with_annotation))
        b = len(non_sign_genes.intersection(genes_with_annotation))

        c = len(sign_genes.intersection(genes_without_annotation))
        d = len(non_sign_genes.intersection(genes_without_annotation))

        # https://www.pathwaycommons.org/guide/primers/statistics/fishers_exact_test/
        _, pvalue = fisher_exact(
            [[a, b], [c, d]],
            alternative=alternative
        )

        in_significant = len(sign_genes.intersection(genes_with_annotation))
        in_background = len(background.intersection(genes_with_annotation))

        if len(sign_genes) > 0:
            fract_of_significant_w_annot = in_significant / len(sign_genes)
            fract_of_background_w_annot = in_background / len(background)

            if fract_of_background_w_annot == 0:
                raise ValueError('something wrong')

            enrichment = np.log2(
                fract_of_significant_w_annot / fract_of_background_w_annot)
        else:
            enrichment = -np.inf

        collector.loc[annotation, 'fold_enrichment'] = enrichment
        collector.loc[annotation, 'pvalue'] = pvalue
        collector.loc[annotation, 'in_significant'] = in_significant
        collector.loc[annotation, 'in_background'] = in_background

    collector.loc[:, 'bonferroni'] = multipletests(
        collector.loc[:, 'pvalue'].values, method='bonferroni')[1]
    collector.loc[:, 'benjamini_hochberg'] = multipletests(
        collector.loc[:, 'pvalue'].values, method='fdr_bh')[1]

    collector = collector.rename_axis('annotation')

    collector = collector.sort_values('pvalue')

    return collector



def functional_ratio_no_filter(
        significant_genes,
        background,
        painter,
        test='enrichment'):
    """
    Computes functional enrichments, and does so with directly interpretable
    Fisher test, that would also allow to test for depletion,

    Input
        significant_genes   list of entrez gene IDs found by user
        background       list of entrez gene IDs, background considered
                            by user ("significant" and "non-signficant"
                            genes)
        painter         df, with ['gene_ncbi', 'annotation']
        test            'enrichment' (test only for enrichment),
                        or 'difference' (test for enrichments or depletion)
    """

    if test == 'enrichment':
        alternative = 'greater'
    elif test == 'difference':
        alternative = 'two-sided'
    else:
        raise ValueError('Test only supports enrichment or difference')

    sign_genes = set(significant_genes)
    background = set(background)
    painter = painter.drop_duplicates()

    if len(sign_genes - set(background)) > 0:
        raise ValueError(
            (
                'The list of significant entries contains some'
                'entries absent from the background list.'
            )
        )
    painter = painter[painter['gene_ncbi'].isin(background)]
    non_sign_genes = set(background) - sign_genes

    minimally_required_candidates_in_annotation = 2
    c = painter['annotation'].value_counts()
    painter = painter[painter['annotation'].isin(
        c[c >= minimally_required_candidates_in_annotation].index)]
    painter = painter.set_index('annotation')

    annotations = painter.index.unique()

    collector = pd.DataFrame(
        index=annotations,
        columns=[
            'fold_enrichment',
            'pvalue',
            'bonferroni',
            'benjamini_hochberg',
            'in_significant',
            'in_background'],
    )

    if len(annotations) == 1:
        raise AssertionError('Must have multiple annotations')

    for annotation in annotations:

        genes_with_annotation = set(painter.loc[annotation, 'gene_ncbi'])
        genes_without_annotation = background - genes_with_annotation

        a = len(sign_genes.intersection(genes_with_annotation))
        b = len(non_sign_genes.intersection(genes_with_annotation))

        c = len(sign_genes.intersection(genes_without_annotation))
        d = len(non_sign_genes.intersection(genes_without_annotation))

        # https://www.pathwaycommons.org/guide/primers/statistics/fishers_exact_test/
        _, pvalue = fisher_exact(
            [[a, b], [c, d]],
            alternative=alternative
        )

        in_significant = len(sign_genes.intersection(genes_with_annotation))
        in_background = len(background.intersection(genes_with_annotation))

        if len(sign_genes) > 0:
            fract_of_significant_w_annot = in_significant / len(sign_genes)
            fract_of_background_w_annot = in_background / len(background)

            if fract_of_background_w_annot == 0:
                raise ValueError('something wrong')

            enrichment = np.log2(
                fract_of_significant_w_annot / fract_of_background_w_annot)
        else:
            enrichment = -np.inf

        collector.loc[annotation, 'fold_enrichment'] = enrichment
        collector.loc[annotation, 'pvalue'] = pvalue
        collector.loc[annotation, 'in_significant'] = in_significant
        collector.loc[annotation, 'in_background'] = in_background

    collector.loc[:, 'bonferroni'] = multipletests(
        collector.loc[:, 'pvalue'].values, method='bonferroni')[1]
    collector.loc[:, 'benjamini_hochberg'] = multipletests(
        collector.loc[:, 'pvalue'].values, method='fdr_bh')[1]

    collector = collector.rename_axis('annotation')

    collector = collector.sort_values('pvalue')

    return collector



def compute_functional_enrichment(
        sign_genes,
        non_sign_genes,
        labelling,
        restrict_to_labelling=False):
    """
    Will compute functional annotation enrichments, simialar
    to DAVID. Note that this function will allow great flexibilty
    about the genes and annotation lists that will be used;

    Input:
        sign_genes      list of entrez gene IDs found by user
        non_sign_genes  list of entrez gene IDs, part of backround,
                            but not found by user
        labelling       df, with
                        ['gene_ncbi', 'annotation_id', 'annotation_name']
        restrict_to_labelling  optional; if True only sign_genes and
                                non_sign_genes within labelling are used
                                (e.g.: use if target label has only
                                been used for a subset of genes)

    Output:
        functional_enrichment df, with several statistics about enrichemnt
            - fold_enrichment
            - oddsratio
            - ease_score (conservative p_value, with one less observation)
            - bonferroni: bonferroni corrected ease_score
            - fdr: false discovery rate corrected ease_score
            - fraction_of_gene_list

    """

    # Settings

    minimally_required_candidates_in_annotation = 1

    # Ensure that only genes, which are present in reference
    # annotation set will be used
    sign_genes = set(sign_genes)
    non_sign_genes = set(non_sign_genes)
    genes_in_labelling = set(labelling['gene_ncbi'])

    if restrict_to_labelling:
        sign_genes = sign_genes.intersection(genes_in_labelling)
        non_sign_genes = non_sign_genes.intersection(genes_in_labelling)

    # ensure that only genes of either cateogy is used
    total_genes = sign_genes.union(non_sign_genes)
    f = labelling['gene_ncbi'].isin(total_genes)
    labelling = labelling.loc[f, :].drop_duplicates()

    # count annotations within user genes & pathway
    f = labelling['gene_ncbi'].isin(sign_genes)
    count_ug_pathway = labelling.loc[
        f, 'annotation_id'].value_counts()

    count_background_pathway = labelling.loc[
        :, 'annotation_id'].value_counts()

    master = pd.concat(
        [
            count_ug_pathway.to_frame('ug_pathway'),
            count_background_pathway.to_frame('background_pathway'),
        ], sort=True, join='outer', axis=1).fillna(0)

    idx = master.index
    res = pd.DataFrame(
        index=idx,
        columns=[
            'fold_enrichment',
            'oddsratio',
            'ease_score',
            'bonferroni',
            'fdr',
            'fraction_of_gene_list'])

    for i in idx:
        a = master.loc[i, 'ug_pathway']
        b = master.loc[i, 'background_pathway']
        c = len(sign_genes) - a
        d = len(total_genes) - b

        if a < minimally_required_candidates_in_annotation:
            res.loc[i, 'oddsratio'] = np.nan
            res.loc[i, 'fold_enrichment'] = np.nan
            res.loc[i, 'ease_score'] = 1  # Set to non-significant
            if len(sign_genes) > 0:
                res.loc[i, 'fraction_of_gene_list'] = a / len(sign_genes)
            else:
                res.loc[i, 'fraction_of_gene_list'] = np.nan
        else:
            res.loc[i, 'fraction_of_gene_list'] = a / len(sign_genes)
            oddsratio, _ = fisher_exact([[a, b], [c, d]])

            # obtain statistics from modified Fisher Exact
            # (EASE score) rather than p-values:
            # Subtract 1 from count of pathway genes in
            # user selection
            # https://david.ncifcrf.gov/content.jsp?file=functional_annotation.html
            e = a - 1
            if e < 0:  # don't allow negative counts
                e = 0
            _, ease = fisher_exact([[e, b], [c, d]])

            res.loc[i, 'fold_enrichment'] = np.log2(
                (a / len(sign_genes)) / (b / len(total_genes)))
            res.loc[i, 'oddsratio'] = oddsratio
            if res.loc[i, 'fold_enrichment'] > 0:
                res.loc[i, 'ease_score'] = ease
            else:  # no enrichment
                res.loc[i, 'ease_score'] = 1

    if res.shape[0] > 1:
        res.loc[:, 'bonferroni'] = multipletests(
            res.loc[:, 'ease_score'].values, method='bonferroni')[1]
        res.loc[:, 'fdr'] = multipletests(
            res.loc[:, 'ease_score'].values, method='fdr_bh')[1]

    functional_enrichment = pd.merge(
        res.sort_values(
            'ease_score', ascending=True),
        labelling[
            ['annotation_id', 'annotation_name']].drop_duplicates(),
        left_index=True,
        right_on='annotation_id',
        how='left')

    functional_enrichment = functional_enrichment.set_index(
        'annotation_id',
        verify_integrity=True)

    return functional_enrichment
