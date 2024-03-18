from __future__ import annotations
from functools import reduce
import math
from collections import defaultdict

from sampling_fo2.fol.syntax import AUXILIARY_PRED_NAME
from sampling_fo2.network.constraint import PartitionConstraint
from sampling_fo2.utils.multinomial import MultinomialCoefficients
from .syntax import *

PREDICATES = defaultdict(list)

PERMITTED_VAR_NAMES = range(ord('A'), ord('Z') + 1)
def new_var(exclude: frozenset[Var]) -> Var:
    for c in PERMITTED_VAR_NAMES:
        v = Var(chr(c))
        if v not in exclude:
            return v
    raise RuntimeError(
        "No more variables available"
    )


def new_predicate(arity: int, name: str) -> Pred:
    global PREDICATES
    p = Pred('{}{}'.format(name, len(PREDICATES[name])), arity)
    PREDICATES[name].append(p)
    return p


def get_predicates(name: str) -> list[Pred]:
    return PREDICATES[name]


def new_scott_predicate(arity: int) -> Pred:
    return new_predicate(arity, SCOTT_PREDICATE_PREFIX)


def pad_vars(vars: frozenset[Var], arity: int) -> frozenset[Var]:
    if arity > 3:
        raise RuntimeError(
            "Not support arity > 3"
        )
    ret_vars = set(vars)
    default_vars = [X, Y, Z]
    idx = 0
    while(len(ret_vars) < arity):
        ret_vars.add(default_vars[idx])
        idx += 1
    return frozenset(list(ret_vars)[:arity])


def exactly_one_qf(preds: list[Pred]) -> QFFormula:
    if len(preds) == 1:
        return top
    lits = [p(X) for p in preds]
    # p1(x) v p2(x) v ... v pm(x)
    formula = reduce(lambda x, y: x | y, lits)
    for i, l1 in enumerate(lits):
        for j, l2 in enumerate(lits):
            if i < j:
                formula = formula & ~(l1 & l2)
    return formula


def exactly_one(preds: list[Pred]) -> QuantifiedFormula:
    if len(preds) == 1:
        return top
    return QuantifiedFormula(Universal(X), exactly_one_qf(preds))

def convert_counting_formula(formula: QuantifiedFormula, domain: set[Const]):
    """
    Only need to deal with \forall X \exists_{=k} Y: f(X,Y)
    """
    uni_formula = top
    ext_formulas = []

    cnt_quantified_formula = formula.quantified_formula.quantified_formula
    cnt_quantifier = formula.quantified_formula.quantifier_scope
    count_param = cnt_quantifier.count_param

    repeat_factor = (math.factorial(count_param)) ** len(domain)

    # Follow the steps in "A Complexity Upper Bound for
    # Some Weighted First-Order Model Counting Problems With Counting Quantifiers"
    # (2)
    aux_pred = new_predicate(2, AUXILIARY_PRED_NAME)
    aux_atom = aux_pred(X, Y)
    uni_formula = uni_formula & (cnt_quantified_formula.equivalent(aux_atom))
    # (3)
    sub_aux_preds, sub_aux_atoms = [], []
    for i in range(count_param):
        aux_pred_i = new_predicate(2, f'{aux_pred.name}_')
        aux_atom_i = aux_pred_i(X, Y)
        sub_aux_preds.append(aux_pred_i)
        sub_aux_atoms.append(aux_atom_i)
        sub_ext_formula = QuantifiedFormula(Existential(Y), aux_atom_i)
        sub_ext_formula = QuantifiedFormula(Universal(X), sub_ext_formula)
        ext_formulas.append(sub_ext_formula)
    # (4)
    for i in range(count_param):
        for j in range(i):
            uni_formula = uni_formula & (~sub_aux_atoms[i] | ~sub_aux_atoms[j])
    # (5)
    or_sub_aux_atoms = QFFormula(False)
    for atom in sub_aux_atoms:
        or_sub_aux_atoms = or_sub_aux_atoms | atom
    uni_formula = uni_formula & or_sub_aux_atoms.equivalent(aux_atom)
    # (6)
    cardinality_constraint = (aux_pred, '=', len(domain) * count_param)

    return uni_formula, ext_formulas, cardinality_constraint, repeat_factor


def unary_evidence_to_ccs(evidence: set[AtomicFormula], domain: set[Const]) \
        -> tuple[QFFormula, PartitionConstraint]:
    """
    Convert unary evidence to cardinality constraints
    """
    assert(all(len(atom.args) == 1 for atom in evidence)), \
        "Only support unary evidence for now"
    element2lits = defaultdict(list)
    for atom in evidence:
        element2lits[atom.args[0]].append(atom.substitute({atom.args[0]: X}))
    evi_size = defaultdict(int)
    for _, lits in element2lits.items():
        evi_size[frozenset(lits)] += 1
    formula = top
    aux_preds = []
    ccs = list()
    for lits, size in evi_size.items():
        aux_pred = new_predicate(1, AUXILIARY_PRED_NAME)
        aux_preds.append(aux_pred)
        aux_atom = aux_pred(X)
        formula = formula & (
            aux_atom.implies(reduce(lambda x, y: x & y, lits))
        )
        ccs.append((aux_pred, '=', size))
    ns = tuple(n for _, _, n in ccs)
    n_sum = sum(ns)
    repeat_factor = (
        MultinomialCoefficients.coef((n_sum, len(domain) - n_sum))
        * MultinomialCoefficients.coef(ns)
    )
    formula = formula & exactly_one_qf(aux_preds)
    return formula, ccs, repeat_factor


def unary_evidence_to_pc(evidence: set[AtomicFormula]) \
        -> tuple[QFFormula, PartitionConstraint]:
    """
    Convert unary evidence to partition constraint
    """
    assert(all(len(atom.args) == 1 for atom in evidence)), \
        "Only support unary evidence for now"
    element2lits = defaultdict(list)
    for atom in evidence:
        element2lits[atom.args[0]].append(atom.substitute({atom.args[0]: X}))
    evi_size = defaultdict(int)
    for _, lits in element2lits.items():
        evi_size[frozenset(lits)] += 1
    formula = top
    aux_preds = []
    partition = list()
    for lits, size in evi_size.items():
        aux_pred = new_predicate(1, AUXILIARY_PRED_NAME)
        aux_preds.append(aux_pred)
        aux_atom = aux_pred(X)
        formula = formula & (
            aux_atom.implies(reduce(lambda x, y: x & y, lits))
        )
        partition.append((aux_pred, size))
    formula = formula & exactly_one_qf(aux_preds)
    return formula, PartitionConstraint(partition)
