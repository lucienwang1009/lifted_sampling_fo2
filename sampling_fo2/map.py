import argparse
from collections import OrderedDict
from decimal import Decimal
from fractions import Fraction
import logging
import math
import logzero
import numpy as np

from symengine.functions import exp
from logzero import logger
from sampling_fo2.utils.polynomial import round_rational
from sampling_fo2.wfomc import Algo, count_distribution
from sampling_fo2.context.wfomc_context import WFOMCContext
from sampling_fo2.fol.sc2 import to_sc2
from sampling_fo2.fol.syntax import AUXILIARY_PRED_NAME, Equivalence, Pred, QuantifiedFormula, top, Universal
from sampling_fo2.fol.utils import new_predicate
from sampling_fo2.network.mln import MLN
from sampling_fo2.parser import mln_parse
from sampling_fo2.problems import WFOMCSProblem
from sampling_fo2.utils import Rational


def parse_args():
    parser = argparse.ArgumentParser(
        description='Compute the maximum and minimum log probability of a given MLN',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('--input', '-i', type=str, required=True,
                        help='mln file')
    parser.add_argument('--debug', action='store_true', default=False)
    args = parser.parse_args()
    return args


def to_wfomcs(mln: MLN) -> tuple[WFOMCSProblem, list[Pred], list[float]]:
    sentence = top
    weightings: dict[Pred, tuple[Rational, Rational]] = dict()
    aux_preds = []
    formula_weight = []
    for weighted_formula in mln:
        formula, weight = weighted_formula.formula, weighted_formula.weight
        free_vars = formula.free_vars()
        if not weighted_formula.is_hard():
            aux_pred = new_predicate(len(free_vars), AUXILIARY_PRED_NAME)
            aux_preds.append(aux_pred)
            formula_weight.append(weight)
            formula = Equivalence(formula, aux_pred(*free_vars))
            exp_weight = Rational(Fraction(math.exp(weight)).numerator,
                                  Fraction(math.exp(weight)).denominator)
            weightings[aux_pred] = (exp_weight, Rational(1, 1))
        for free_var in free_vars:
            formula = QuantifiedFormula(Universal(free_var), formula)
        sentence = sentence & formula
    try:
        sentence = to_sc2(sentence)
    except:
        raise ValueError('Sentence must be a valid SC2 formula.')
    return WFOMCSProblem(sentence, mln.domain, weightings,
                         mln.cardinality_constraint), aux_preds, formula_weight


def main(mln_file: str) -> tuple[Decimal, Decimal]:
    if not mln_file.endswith('.mln'):
        raise RuntimeError(f'Only support .mln file: {input_file}')
    with open(mln_file, 'r') as f:
        mln = mln_parse(f.read())
    problem, aux_preds, formula_weight = to_wfomcs(mln)
    context = WFOMCContext(problem)
    cdist = count_distribution(context, aux_preds, Algo.FASTERv2)
    wfomc = sum(cdist.values())
    wfomc_ln = float(round_rational(wfomc).ln())
    max_w = float('-inf')
    min_w = float('inf')
    for cc, weight in cdist.items():
        if weight == 0:
            continue
        w = np.dot(cc, formula_weight)
        if w > max_w:
            max_w = w
        if w < min_w:
            min_w = w
    max_log_prob = max_w - wfomc_ln
    min_log_prob = min_w - wfomc_ln
    return max_log_prob, min_log_prob


if __name__ == '__main__':
    args = parse_args()
    if args.debug:
        logzero.loglevel(logging.DEBUG)
    else:
        logzero.loglevel(logging.INFO)
    input_file = args.input
    max_prob, min_prob = main(input_file)
    logger.info(f'Maximum probability: {max_prob}')
    logger.info(f'Minimum probability: {min_prob}')
