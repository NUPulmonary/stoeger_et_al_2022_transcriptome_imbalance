
import os
import sys

import numpy as np
import pandas as pd

sys.path.append('./../src/')
from aging_tools import inout


def brehme_human_to_mouse():
    """
    Loads chaperome definition from Brehme,
    mapped to mouse homologs
    """

    p = inout.get_internal_path(
        'datasets/general/chaperome/brehme_human_with_mouse_homologs.csv')

    df_chaperome = pd.read_csv(p)

    return df_chaperome


def antalek_human_to_mouse():
    """
    Loads chaperome definition from Mat Antalek,
    mapped to mouse homologs
    """

    p = inout.get_internal_path(
        'datasets/general/chaperome/antalek_human_with_mouse_homologs.csv')

    df_chaperome = pd.read_csv(p)

    return df_chaperome
