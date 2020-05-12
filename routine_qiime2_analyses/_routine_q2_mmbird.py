# ----------------------------------------------------------------------------
# Copyright (c) 2020, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

# import os
# import itertools
# import pandas as pd
# from os.path import basename, isfile, splitext
# import multiprocessing


def run_mmbird(i_datasets_folder: str, datasets: dict, taxonomies: dict,
               songbird_outputs: dict, mmvec_outputs: dict,
               prjct_nm: str, qiime_env: str, chmod: str, noloc: bool) -> None:
    print()
    print()
    print("i_datasets_folder")
    print(i_datasets_folder)
    print()
    print()
    print("songbird_outputs")
    print(songbird_outputs)
    # {
    #   'vioscreen_foods_consumed_grams_per_day_1800s_noLiquids':
    #      [
    #         [
    #            'vioscreen_foods_consumed_grams_per_day_1800s_noLiquids',
    #            '.../nut/metadata/meta_vioscreen_foods_consumed_grams_per_day_1800s_noLiquids_alphas.tsv',
    #            '.../data/tab_vioscreen_foods_consumed_grams_per_day_1800s_noLiquids.qza',
    #            'age',
    #            'vioscreen_foods_consumed_grams_per_day_1800s_noLiquids/age/filt_f0_s0/2_1e-4_100_05'
    #         ]
    #      ],
    #   'vioscreen_micromacro_qemistree_1800s':
    #      [
    #         [
    #            'vioscreen_micromacro_qemistree_1800s',
    #            '.../nut/metadata/meta_vioscreen_micromacro_qemistree_1800s_alphas.tsv',
    #            .../nut/data/tab_vioscreen_micromacro_qemistree_1800s.qza',
    #            'age',
    #            'vioscreen_micromacro_qemistree_1800s/age/filt_f0_s0/2_1e-4_100_05'
    #         ]
    #      ]
    #   }
    print()
    print()
    print("mmvec_outputs")
    print(mmvec_outputs)
    # {'grams_qemi':
    #   [
    #     [
    #        '.../metadata/grams_qemi/meta_grams_qemi__10_0_1800s.tsv',
    #        '.../data/grams_qemi/tab_vioscreen_foods_consumed_grams_per_day_1800s_noLiquids__grams_qemi__10_0_1800s.qza',
    #        '.../data/grams_qemi/tab_vioscreen_micromacro_qemistree_1800s__grams_qemi__10_0_1800s.qza',
    #        'None',
    #        'b-2_l-1e-4_e-5000_p-05_f-0_d-3_t-None_n-100_gpu-F'
    #     ]
    #   ]
    # }
    pass
