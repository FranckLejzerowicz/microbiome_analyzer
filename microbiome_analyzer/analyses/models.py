# ----------------------------------------------------------------------------
# Copyright (c) 2023, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import re
import sys
from os.path import dirname


def split_formula(formula):
    formula = formula.strip('"').strip("'")
    formula_split = re.split('[+/:*]', formula)
    return formula_split


def parse_formula(formula):
    variables = {}
    for term in split_formula(formula):
        if term[:2] == 'C(':
            term, factors = parse_cat_term(term)
            variables.setdefault(term, []).extend(factors)
        else:
            if term not in variables:
                variables[term] = []
    return variables


def parse_cat_term(term):
    term, ref = [x.strip() for x in term[2:-1].split(',', 1)]
    factors = get_factors(ref, term)
    return term, factors


def get_factors(ref, term):
    if ref.startswith("Treatment"):
        factors = [ref[11:-2]]
    elif ref.startswith("Diff"):
        factors = eval(re.split('levels ?= ?', ref)[-1])
    else:
        print('[songbird] Formatting issue with model term "%s"' % term)
        sys.exit(1)
    return factors


def show_models_issues(self, mess):
    if mess:
        for m in sorted(mess):
            print(m)
    if self.models_issues:
        print('\n## Issues with model (will not run)')
        for model_issue, metas in self.models_issues.items():
            print('# -', model_issue)
            for meta in sorted(metas):
                print('#\t.%s' % meta.replace(dirname(self.dir), ''))
        print('\n')
