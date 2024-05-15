import argparse
from collections import OrderedDict
from decimal import Decimal
from fractions import Fraction
import logging
import math
import logzero
import numpy as np

from symengine.functions import log
from logzero import logger
from sampling_fo2.utils.polynomial import round_rational, to_rational
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
        description='Compute the entropy of a given MLN',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('--input', '-i', type=str, required=True,
                        help='mln file')
    parser.add_argument('--debug', action='store_true', default=False)
    args = parser.parse_args()
    return args


def to_wfomcs(mln: MLN) -> tuple[WFOMCSProblem, list[Pred], list[Decimal]]:
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
            formula_weight.append(Decimal(weight))
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


def main(mln_file: str) -> float:
    if not mln_file.endswith('.mln'):
        raise RuntimeError(f'Only support .mln file: {input_file}')
    with open(mln_file, 'r') as f:
        mln = mln_parse(f.read())
    problem, aux_preds, formula_weight = to_wfomcs(mln)
    # reset the weight for every predicate to 1
    problem.weights = dict()
    context = WFOMCContext(problem)
    cdist = count_distribution(context, aux_preds, Algo.FASTERv2)
    wfomc = Decimal(0)
    unnormalized = Decimal(0)
    for cc, num in cdist.items():
        cc_w_dot = sum(
            Decimal(int(n)) * w for n, w in zip(cc, formula_weight)
        )
        unnormalized = unnormalized + (
            cc_w_dot.exp() * Decimal(int(num)) * cc_w_dot
        )
        wfomc = wfomc + cc_w_dot.exp() * Decimal(int(num))
    entropy = wfomc.ln() - unnormalized / wfomc
    return entropy


if __name__ == '__main__':
    args = parse_args()
    if args.debug:
        logzero.loglevel(logging.DEBUG)
    else:
        logzero.loglevel(logging.INFO)
    input_file = args.input
    entropy = main(input_file)
    logger.info(f'Entropy: {entropy}')
