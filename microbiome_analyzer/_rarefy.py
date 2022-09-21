# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import glob
import numpy as np
import pandas as pd
from scipy.stats import skew
np.set_printoptions(precision=2, suppress=True)


def get_dat_depths(
        dat: str,
        depths_yml: dict,
        output_folder: str,
        sam_sum: pd.Series) -> tuple:
    """
    Parameters
    ----------
    dat : str
        Dataset name
    depths_yml : dist
        Mapping Dataset nanme -> Depths at which to rarefy
    output_folder : str
        Path to the output folder
    sam_sum : pd.Series
        Sum of reads per sample

    Returns
    -------
    skip : bool
        Whether to skip a rarefaction
    depths_tuple : tuple
        (boolean, Rarefaction depths)
    """
    skip = False
    if not depths_yml:
        depths = get_default_raref_depth(dat, output_folder, sam_sum)
        depths_tuple = (0, depths)
    elif dat in depths_yml:
        skip, depths = get_depths(dat, depths_yml[dat], sam_sum)
        depths_tuple = (1, depths)
    else:
        skip = True
        depths_tuple = []
    return skip, depths_tuple


def get_default_raref_depth(
        dat: str,
        output_folder: str,
        sam_sum: pd.Series):
    """
    Parameters
    ----------
    dat : str
        Dataset name
    output_folder : str
        Path to the output folder
    sam_sum : pd.Series
        Sum of reads per sample

    Returns
    -------
    depths : list
        Rarefaction depths
    """
    raref_files = glob.glob('%s/rarefy/%s/tab_raref*.qza' % (
        output_folder, dat))
    if len(raref_files):
        depths = [x.split('_raref')[-1].split('.tsv')[0] for x in raref_files]
    else:
        second_quantile = sam_sum.quantile(0.2)
        print_skew(dat, sam_sum)
        if second_quantile < 1000:
            depths = []
            print_not_rarefying(dat, sam_sum)
        else:
            nfigure = len(str(int(second_quantile)))
            second_quantile_to_round = second_quantile / (10 ** (nfigure - 2))
            second_quantile_rounded = round(second_quantile_to_round) * (
                        10 ** (nfigure - 2))
            depths = [str(int(second_quantile_rounded))]
            print('[%s] Proposed rarefaction depth: %s '
                  '(second quantile)' % (dat, depths[0]))
    return depths


def get_depths(
        dat: str,
        depths_yml: list,
        sam_sum: pd.Series) -> tuple:
    """

    Parameters
    ----------
    dat : str
        Dataset name
    depths_yml : list
        Depths at which to rarefy
    sam_sum : pd.Series
        Sum of reads per sample

    Returns
    -------
    skip : bool
        Whether to skip a rarefaction
    depths : list
        Depths at which to rarefy
    """
    skip = False
    depths = []
    for depth in depths_yml:
        if depth == 'min' or sum(sam_sum >= int(depth)) > 10:
            depths.append(depth)
    if not depths:
        print('[%s] Min. proposed rarefaction depths would leave <10 samples: '
              '%s (not rarefaction)' % (dat, ', '.join(depths_yml)))
        skip = True
    elif len(depths) != len(depths_yml):
        print('[%s] Proposed rarefaction depths would leave <10 samples: %s ('
              'not rarefied)' % (dat, ', '.join([x for x in depths_yml
                                                 if x not in depths])))
    return skip, depths


def print_skew(
        dat: str,
        tsv_sam_sum: pd.Series) -> None:
    """
    Parameters
    ----------
    dat : str
        Dataset name
    tsv_sam_sum : pd.Series
        Sum of reads per sample
    """
    count, division = np.histogram(tsv_sam_sum)
    skw = skew(count)
    if abs(skw) > 1:
        print()
        print(' ==> Consider rarefying <==')
        print('[%s] Reads-per-sample distribution [skewness=%s] (>1!)' % (
            dat, round(abs(float(skw)), 3)))
        division_std = np.interp(
            count, (min(count), max(count)), (0, 20))
        print('\treadsbin\tsamples\thistogram')
        for ddx, div in enumerate(division_std):
            if div > 1:
                print('\t%s\t%s\t%s' % (
                format(division[ddx], '6.3E'), count[ddx], '-' * int(div)))
            elif div == 0:
                print('\t%s\t%s\t%s' % (
                format(division[ddx], '6.3E'), count[ddx], ''))
            else:
                print('\t%s\t%s\t%s' % (
                format(division[ddx], '6.3E'), count[ddx], '-'))


def print_not_rarefying(
        dat: str,
        sam_sum: pd.Series) -> None:
    """
    Parameters
    ----------
    dat: str
        Dataset name
    sam_sum : pd.Series
        Sum of reads per sample
    """
    print('[%s] Second quantile of the reads-per-sample '
          'distribution is <1000' % dat)
    print('- The sequencing might have failed! Analyze with caution')
    print('- reads-per-sample distribution described:')
    for x, y in sam_sum.describe().to_dict().items():
        print('\t%s: %s' % (x, round(y, 3)))
    print('!!! NOT RAREFYING %s !!!' % dat)


def get_digit_depth(
        depth_: str,
        tsv_sums: pd.Series) -> int:
    """Get the rarefaction depth integer

    Parameters
    ----------
    depth_ : str
        Rarefaction depth
    tsv_sums : pd.Series
        Sum of reads per sample

    Returns
    -------
    depth : int
        Rarefaction depth
    """
    if depth_.isdigit():
        depth = int(depth_)
    else:
        depth = int(np.floor(min(tsv_sums)))
    return depth
