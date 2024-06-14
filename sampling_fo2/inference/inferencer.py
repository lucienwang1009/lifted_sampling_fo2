from decimal import Decimal
import logzero
import numpy as np
import logging

from collections import defaultdict
from logzero import logger
from math import exp

from sampling_fo2.context.wfomc_context import WFOMCContext
from sampling_fo2.fol.syntax import AtomicFormula, Const, Existential, QFFormula, \
    QuantifiedFormula, Quantifier, Universal, a, b
from sampling_fo2.inference.db import Database
from sampling_fo2.network.constraint import CardinalityConstraint
from sampling_fo2.network.mln import MLN, WeightedFormula
from sampling_fo2.problems import MLN_to_WFOMC
from sampling_fo2.utils.polynomial import round_rational
from sampling_fo2.wfomc import Algo, wfomc


class Inferencer:
    def __init__(self, model: MLN):
        self.model: MLN = model

        self.weighted_formulas: WeightedFormula = model.weighted_formulas
        self.domain: set[Const] = model.domain
        self.ccs: CardinalityConstraint = model.cardinality_constraint
        self.preds = list(model.preds)

        # extracted grounding of quantified formulas
        # all formulas are in the form of quantifier-free FO2,
        self.ground_formulas: list[list[QFFormula]] = []
        self.all_ground_atoms: list[AtomicFormula] = []
        self.formulas_quantifiers: list[list[Quantifier]] = []
        # NOTE(lucien): use Decimal here, it would be more accurate
        self.partition_function_ln: Decimal = 0.0

        self.preprocess()

    def preprocess(self):
        # compute partition function
        loglevel = logzero._loglevel
        logzero.loglevel(logging.ERROR)
        wfomc_problem = MLN_to_WFOMC(self.model)
        context = WFOMCContext(wfomc_problem)
        # disable logging when computing wfomc
        self.partition_function_ln = round_rational(wfomc(context, Algo.FASTER)).ln()
        logzero.loglevel(loglevel)
        logger.debug(
            'Partition function in ln: %s', self.partition_function_ln
        )

        for weighted_formula in self.model:
            _, formula = weighted_formula.weight, weighted_formula.formula
            if len(formula.vars()) > 2 or len(formula.vars()) == 0:
                raise ValueError('Formula with arity > 2 or = 0 is not supported')
            qf_formula = formula
            quantified_vars = []
            quantifiers = []
            while(isinstance(qf_formula, QuantifiedFormula)):
                quantifier = qf_formula.quantifier_scope
                quantifiers.append(quantifier)
                quantified_vars.append(
                    quantifier.quantified_var
                )
                qf_formula = qf_formula.quantified_formula
            self.formulas_quantifiers.append(quantifiers)
            if not isinstance(qf_formula, QFFormula):
                raise ValueError(
                    'Only support Vx, Ex, VxEy, Ux, Uy, UxUy, ExEy, ExUy, EyUx, '
                    f'but get {qf_formula}'
                )
            all_vars = list(formula.free_vars()) + quantified_vars
            substituted_vars = ([a, b])[:len(all_vars)]
            grounding_formulas = [
                qf_formula.substitute(
                    dict(zip(all_vars, substituted_vars))
                )
            ]
            if len(qf_formula.vars()) == 2:
                grounding_formulas.append(
                    qf_formula.substitute(
                        dict(zip(all_vars, [a, a]))
                    )
                )
            self.ground_formulas.append(grounding_formulas)
        logger.debug(
            'Ground formulas: %s', self.ground_formulas
        )

        for pred in self.preds:
            if pred.arity == 1:
                self.all_ground_atoms.append(pred(a))
                self.all_ground_atoms.append(pred(b))
            elif pred.arity == 2:
                self.all_ground_atoms.append(pred(a, b))
                self.all_ground_atoms.append(pred(b, a))
                self.all_ground_atoms.append(pred(a, a))
                self.all_ground_atoms.append(pred(b, b))
            else:
                raise ValueError('Predicate with arity > 2 or = 0 is not supported')
        logger.debug(
            'All ground atoms in the MLN: %s', self.all_ground_atoms
        )

        if self.ccs:
            self.ccs.build()

    def satisfied_groundings(self, db: Database) -> list[int]:
        """
        Compute the number of groundings of each formula that are satisfied by the database.
        :param db: The database.
        :return: The number of groundings of each formula that are satisfied by the database.
        """
        if db.preds.difference(self.preds):
            raise ValueError('Database contains predicates not in the MLN: %s',
                             db.preds.difference(self.preds))
        if db.domain.difference(self.domain):
            raise ValueError('Database contains constants not in the MLN: %s',
                             db.domain.difference(self.domain))
        truth_substs = dict()
        for i, c1 in enumerate(self.domain):
            for j, c2 in enumerate(self.domain):
                if i >= j:
                    continue
                substitution = {a: c1, b: c2}
                reverse_substitution = {c1: a, c2: b}
                ground_atoms = set(
                    atom.substitute(substitution)
                    for atom in self.all_ground_atoms
                )
                truth_subst = dict()
                for atom in ground_atoms:
                    reverse_atom = atom.substitute(reverse_substitution)
                    if atom in db.facts:
                        truth_subst[reverse_atom] = True
                    else:
                        truth_subst[reverse_atom] = False
                truth_substs[(c1, c2)] = truth_subst
                substitution = {a: c2, b: c1}
                reverse_substitution = {c2: a, c1: b}
                ground_atoms = set(
                    atom.substitute(substitution)
                    for atom in self.all_ground_atoms
                )
                truth_subst = dict()
                for atom in ground_atoms:
                    reverse_atom = atom.substitute(reverse_substitution)
                    if atom in db.facts:
                        truth_subst[reverse_atom] = True
                    else:
                        truth_subst[reverse_atom] = False
                truth_substs[(c2, c1)] = truth_subst
        for c in self.domain:
            substitution = {a: c, b: c}
            reverse_substitution = {c: a}
            ground_atoms = set(
                atom.substitute(substitution)
                for atom in self.all_ground_atoms
            )
            truth_subst = dict()
            for atom in ground_atoms:
                reverse_atom = atom.substitute(reverse_substitution)
                if atom in db.facts:
                    truth_subst[reverse_atom] = True
                else:
                    truth_subst[reverse_atom] = False
            truth_substs[(c, c)] = truth_subst
        # satisfications of formulas on each pair of constants
        sats = []
        for formulas in self.ground_formulas:
            if len(formulas) == 2:
                formula_ab, formula_aa = formulas
                sat = np.zeros((len(self.domain), len(self.domain)), dtype=bool)
                for i, c1 in enumerate(self.domain):
                    for j, c2 in enumerate(self.domain):
                        if i != j:
                            sat[i, j] = formula_ab.eval(truth_substs[(c1, c2)])
                        else:
                            sat[i, j] = formula_aa.eval(truth_substs[(c1, c2)])
                sats.append(sat)
            else:
                formula = formulas[0]
                sat = np.zeros((len(self.domain), ), dtype=bool)
                for i, c in enumerate(self.domain):
                    sat[i] = formula.eval(truth_substs[(c, c)])
                sats.append(sat)
        num_satisfied = []
        for i, weighted_formula in enumerate(self.model):
            formula = weighted_formula.formula
            all_sat = sats[i]
            quantifiers = self.formulas_quantifiers[i]
            if weighted_formula.is_hard():
                continue
            for quantifier in quantifiers[::-1]:
                if isinstance(quantifier, Universal):
                    all_sat = np.all(all_sat, axis=-1)
                if isinstance(quantifier, Existential):
                    all_sat = np.any(all_sat, axis=-1)
            N = all_sat.sum()
            num_satisfied.append(N)
        return num_satisfied

    def log_prob(self, db: Database):
        """
        Compute the log probability of the database given the MLN.
        :param db: The database.
        :return: The log probability.
        """
        if db.preds.difference(self.preds):
            raise ValueError('Database contains predicates not in the MLN: %s',
                             db.preds.difference(self.preds))
        if db.domain.difference(self.domain):
            raise ValueError('Database contains constants not in the MLN: %s',
                             db.domain.difference(self.domain))
        truth_substs = dict()
        cardinalities = defaultdict(set)
        for i, c1 in enumerate(self.domain):
            for j, c2 in enumerate(self.domain):
                if i >= j:
                    continue
                substitution = {a: c1, b: c2}
                reverse_substitution = {c1: a, c2: b}
                ground_atoms = set(
                    atom.substitute(substitution)
                    for atom in self.all_ground_atoms
                )
                truth_subst = dict()
                for atom in ground_atoms:
                    reverse_atom = atom.substitute(reverse_substitution)
                    if atom in db.facts:
                        truth_subst[reverse_atom] = True
                        cardinalities[atom.pred].add(atom)
                    else:
                        truth_subst[reverse_atom] = False
                truth_substs[(c1, c2)] = truth_subst
                substitution = {a: c2, b: c1}
                reverse_substitution = {c2: a, c1: b}
                ground_atoms = set(
                    atom.substitute(substitution)
                    for atom in self.all_ground_atoms
                )
                truth_subst = dict()
                for atom in ground_atoms:
                    reverse_atom = atom.substitute(reverse_substitution)
                    if atom in db.facts:
                        truth_subst[reverse_atom] = True
                        cardinalities[atom.pred].add(atom)
                    else:
                        truth_subst[reverse_atom] = False
                truth_substs[(c2, c1)] = truth_subst
        for c in self.domain:
            substitution = {a: c, b: c}
            reverse_substitution = {c: a}
            ground_atoms = set(
                atom.substitute(substitution)
                for atom in self.all_ground_atoms
            )
            truth_subst = dict()
            for atom in ground_atoms:
                reverse_atom = atom.substitute(reverse_substitution)
                if atom in db.facts:
                    truth_subst[reverse_atom] = True
                    cardinalities[atom.pred].add(atom)
                else:
                    truth_subst[reverse_atom] = False
            truth_substs[(c, c)] = truth_subst
        cardinalities = dict(
            (pred, len(atoms))
            for pred, atoms in cardinalities.items()
        )
        logger.debug(
            'Cardinalities: %s', cardinalities
        )

        # satisfications of formulas on each pair of constants
        sats = []
        for formulas in self.ground_formulas:
            if len(formulas) == 2:
                formula_ab, formula_aa = formulas
                sat = np.zeros((len(self.domain), len(self.domain)), dtype=bool)
                for i, c1 in enumerate(self.domain):
                    for j, c2 in enumerate(self.domain):
                        if i != j:
                            sat[i, j] = formula_ab.eval(truth_substs[(c1, c2)])
                        else:
                            sat[i, j] = formula_aa.eval(truth_substs[(c1, c2)])
                sats.append(sat)
            else:
                formula = formulas[0]
                sat = np.zeros((len(self.domain), ), dtype=bool)
                for i, c in enumerate(self.domain):
                    sat[i] = formula.eval(truth_substs[(c, c)])
                sats.append(sat)

        # deal with cardinality constraints first
        if self.ccs:
            kwargs = dict((pred.name, cc) for pred, cc in cardinalities.items())
            if not eval(self.ccs.validator.format(**dict(kwargs))):
                logger.debug('Cardinality constraints are not satisfied')
                return float('-inf')
        # hard rules must be satisfied, i.e.,
        # the formula must be true on all pairs of constants
        total_weight = 0.0
        for i, weighted_formula in enumerate(self.model):
            weight, formula = weighted_formula.weight, weighted_formula.formula
            all_sat = sats[i]
            quantifiers = self.formulas_quantifiers[i]
            if weighted_formula.is_hard():
                for quantifier in quantifiers[::-1]:
                    if isinstance(quantifier, Universal):
                        all_sat = np.all(all_sat, axis=-1)
                    if isinstance(quantifier, Existential):
                        all_sat = np.any(all_sat, axis=-1)
                assert not all_sat.shape
                if not all_sat:
                    logger.debug(f'Hard formula {formula} is not satisfied')
                    return float('-inf')
            else:
                for quantifier in quantifiers[::-1]:
                    if isinstance(quantifier, Universal):
                        all_sat = np.all(all_sat, axis=-1)
                    if isinstance(quantifier, Existential):
                        all_sat = np.any(all_sat, axis=-1)
                N = all_sat.sum()
                total_weight += weight * N
        logger.debug('Total weight: %s', total_weight)
        return Decimal(total_weight) - self.partition_function_ln
