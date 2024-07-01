from __future__ import annotations

from sampling_fo2.fol.sc2 import SC2, to_sc2
from sampling_fo2.fol.syntax import AtomicFormula, Const, Pred, top, AUXILIARY_PRED_NAME, \
    Formula, QuantifiedFormula, Universal, Equivalence
from sampling_fo2.fol.utils import new_predicate
from sampling_fo2.network.constraint import CardinalityConstraint
from sampling_fo2.network.mln import MLN
from sampling_fo2.utils.polynomial import Rational
from fractions import Fraction
import math


class WFOMCSProblem(object):
    """
    A weighted first-order model counting/sampling problem.
    """

    def __init__(self, sentence: SC2,
                 domain: set[Const],
                 weights: dict[Pred, tuple[Rational, Rational]],
                 cardinality_constraint: CardinalityConstraint = None):
        self.domain: set[Const] = domain
        self.sentence: SC2 = sentence
        self.weights: dict[Pred, tuple[Rational, Rational]] = weights
        self.cardinality_constraint: CardinalityConstraint = cardinality_constraint

    def __str__(self) -> str:
        s = ''
        s += 'Domain: \n'
        s += '\t' + str(self.domain) + '\n'
        s += 'Sentence: \n'
        s += '\t' + str(self.sentence) + '\n'
        s += 'Weights: \n'
        s += '\t' + str(self.weights) + '\n'
        if self.cardinality_constraint is not None:
            s += 'Cardinality Constraint: \n'
            s += '\t' + str(self.cardinality_constraint) + '\n'
        return s

    def __repr__(self) -> str:
        return str(self)


def MLN_to_WFOMC(mln: MLN):
    sentence = top
    weightings: dict[Pred, tuple[Rational, Rational]] = dict()
    for weighted_formula in mln:
        weight, formula = weighted_formula.weight, weighted_formula.formula
        free_vars = formula.free_vars()
        if not weighted_formula.is_hard():
            aux_pred = new_predicate(len(free_vars), AUXILIARY_PRED_NAME)
            formula = Equivalence(formula, aux_pred(*free_vars))
            weightings[aux_pred] = (Rational(Fraction(math.exp(weight)).numerator,
                                             Fraction(math.exp(weight)).denominator), Rational(1, 1))
        sentence = sentence & formula

    try:
        sentence = to_sc2(sentence)
    except:
        raise ValueError('Sentence must be a valid SC2 formula.')
    return WFOMCSProblem(sentence, mln.domain, weightings, mln.cardinality_constraint)
