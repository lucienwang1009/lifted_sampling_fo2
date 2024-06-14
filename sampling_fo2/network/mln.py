from __future__ import annotations

from typing import List, Set, FrozenSet

from sampling_fo2.fol.syntax import Formula, QuantifiedFormula, Const, Pred, Var
from sampling_fo2.fol.utils import pad_vars, universally_quantify
from sampling_fo2.network.constraint import CardinalityConstraint


class WeightedFormula(object):
    def __init__(self, formula: Formula, weight: float):
        """
        A formula with a weight.
        :param formula: The formula.
        :param weight: The weight.
        """
        self.formula: Formula = formula
        self.weight: float = weight

        # if it's a hard formula, let's present the formula explicitly,
        # i.e., include the quantifiers in the formula
        if self.is_hard():
            free_vars = self.formula.free_vars()
            if free_vars:
                self.formula = universally_quantify(
                    free_vars, self.formula
                )

    def __str__(self):
        return '{} {}'.format(self.weight, self.formula)

    def is_hard(self) -> bool:
        return self.weight == float('inf')


class MLN(object):
    def __init__(self, weighted_formulas: List[WeightedFormula],
                 domain: Set[Const] = None,
                 cardinality_constraint: CardinalityConstraint = None):
        """
        A Markov Logic Network. Slighly different from the original definition:
        * All hard formulas are quantified explicitly;
        * No predicate definitions are in the input, but infered from the formulas;
        * Domain is included in the MLN.

        :param weighted_formulas: A list of weighted formulas.
        :param domain: The domain of the MLN.
        :param cardinality_constraint: The cardinality constraint in the MLN.
        """
        self.weighted_formulas: List[WeightedFormula] = weighted_formulas
        self.domain = domain
        self.cardinality_constraint: CardinalityConstraint = cardinality_constraint

        # deal with predicate_definition
        self.preds: set[Pred] = set()
        for w_formula in self.weighted_formulas:
            self.preds.update(w_formula.formula.preds())
        self.idx: int

    def append(self, weighted_formula: WeightedFormula):
        self.weighted_formulas.append(weighted_formula)
        self.preds.update(weighted_formula.formula.preds())

    def extend(self, weighted_formulas: List[WeightedFormula]):
        self.weighted_formulas.extend(weighted_formulas)
        for w_formula in weighted_formulas:
            self.preds.update(w_formula.formula.preds())

    def __iter__(self):
        self.idx = 0
        return self

    def __next__(self) -> WeightedFormula:
        if self.idx < self.size():
            ret = self.weighted_formulas[self.idx]
            self.idx += 1
            return ret
        else:
            raise StopIteration

    def vars(self) -> FrozenSet[Var]:
        variables = set()
        for w_formula in self.weighted_formulas:
            variables.update(w_formula.vars())
        return frozenset(variables)

    def size(self) -> int:
        return len(self.weighted_formulas)

    def __str__(self):
        s = ''
        for w_formula in self:
            s += '{}\n'.format(w_formula)
        return s


# class ComplexMLN(MLN):
#     def __init__(self, formulas: List[QuantifiedFormula], weights: List[List[complex]] = None,
#                  domain: Set[Const] = None, predicate_definition: Set[Pred] = None):
#         super().__init__(formulas, weights, domain, predicate_definition)
#
#     def __str__(self):
#         s = ''
#         s += 'domain = {}\n'.format(','.join(
#             str(element) for element in self.domain
#         ))
#         for f, ws in self:
#             s += '{} {}\n'.format(','.join(str(w) for w in ws), f)
#         return s
