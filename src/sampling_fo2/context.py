from __future__ import annotations
import random

from itertools import product
from collections import defaultdict
from copy import deepcopy
from dataclasses import dataclass, field
from typing import Generator
from loguru import logger
from itertools import product
from wfomc import AtomicFormula, Pred, a, b, multinomial, \
    RingElement, top, SC2, exactly_one_qf, new_predicate, \
    convert_counting_formula, Const, Pred, QFFormula, QuantifiedFormula, \
    CardinalityConstraint, Rational, RingElement, EVIDOM_PRED_NAME, X, \
    Y, AUXILIARY_PRED_NAME, SKOLEM_PRED_NAME
from wfomc.cell_graph import Cell


class BlockType(frozenset):
    pass


@dataclass(frozen=True)
class ExistentialTwoTable(object):
    code_ab: tuple[bool, ...]
    code_ba: tuple[bool, ...]
    preds: tuple[Pred, ...]
    atoms: tuple[tuple[AtomicFormula, AtomicFormula]] = field(
        default=None, repr=False, init=False, hash=False, compare=False
    ) # pyright: ignore[reportAssignmentType]

    def __post_init__(self):
        atoms = list()
        for p in self.preds:
            args_ab = [a, b]
            args_ba = [b, a]
            atoms.append((
                p(*args_ab),
                p(*args_ba)
            ))
        object.__setattr__(self, 'atoms', tuple(atoms))

    def get_evidences(self, reverse=False) -> frozenset[AtomicFormula]:
        evidences = set()
        for ab, ba, (atom_ab, atom_ba) in zip(self.code_ab, self.code_ba, self.atoms):
            if not reverse:
                evidences.add(atom_ab if ab else ~atom_ab)
                evidences.add(atom_ba if ba else ~atom_ba)
            else:
                evidences.add(atom_ab if ba else ~atom_ab)
                evidences.add(atom_ba if ab else ~atom_ba)
        return frozenset(evidences)

    def is_positive(self, pred: Pred) -> tuple[bool, bool]:
        idx = self.preds.index(pred)
        return self.code_ab[idx], self.code_ba[idx]

    def __str__(self):
        evidences = self.get_evidences()
        return str(evidences)

    def __repr__(self):
        return str(self)


class ExistentialContext(object):
    def __init__(self, cell_assignment: list[Cell],
                 ext_preds: list[Pred]):
        self.cell_assignment = cell_assignment
        self.cell_elements: dict[Cell, set[int]] = defaultdict(set)
        self.cell_config: dict[Cell, int] = defaultdict(lambda: 0)
        for idx, cell in enumerate(cell_assignment):
            self.cell_elements[cell].add(idx)
            self.cell_config[cell] += 1
        self.cells = list(self.cell_config.keys())
        self.ext_preds = ext_preds

        self.cb_config: dict[tuple[Cell, BlockType],
                             int] = defaultdict(lambda: 0)
        self.cb_elements: dict[
            tuple[Cell, BlockType], set[int]
        ] = defaultdict(set)
        # initialize cb_config and cb_elements
        for cell, num in self.cell_config.items():
            block_type = BlockType()
            for ext in self.ext_preds:
                if not cell.is_positive(ext):
                    block_type = BlockType(block_type | {ext})
            key = (cell, block_type)
            self.cb_config[key] = num
            self.cb_elements[key].update(self.cell_elements[cell])
        logger.debug('initial cell-block config: {}', self.cb_elements)

    def all_satisfied(self) -> bool:
        return all(self.cb_config[(cell, BlockType())] == self.cell_config[cell]
                   for cell in self.cells)

    def select_cell_block_type(self) -> tuple[Cell, BlockType]: # pyright: ignore[reportReturnType]
        if self.all_satisfied():
            raise RuntimeError('all satisfied, cannot select more')
        for (cell, block), num in self.cb_config.items():
            if len(block) != 0 and num > 0:
                return cell, block

    def reduce_element(self, cell: Cell, block: BlockType) -> int:
        self.cell_config[cell] -= 1
        self.cb_config[(cell, block)] -= 1
        element = self.cb_elements[(cell, block)].pop()
        return element

    def reduce_block_type(self, block: BlockType,
                          etable: ExistentialTwoTable,
                          ab_or_ba: int = 0) -> BlockType:
        return BlockType(
            block.difference(set(
            ext_pred for ext_pred in block
            if etable.is_positive(ext_pred)[ab_or_ba]
            ))
        )

    def satisfied(self, block: BlockType,
                  overall_etable_config: dict[ExistentialTwoTable, int]) -> bool:
        if len(block) == 0:
            return True
        for etable, num in overall_etable_config.items():
            if num > 0:
                # ab_p
                block = self.reduce_block_type(block, etable)
            if len(block) == 0:
                return True
        return len(block) == 0

    def iter_etable_config(
        self,
        etable_weights: dict[Cell, dict[ExistentialTwoTable, RingElement]]
    ) -> Generator[dict[tuple[Cell, BlockType], dict[ExistentialTwoTable, int]]]:
        for raw_config in product(
                *list(
                    multinomial(len(etable_weights[cell]), num)
                    for (cell, _), num in self.cb_config.items()
                )
        ):
            etable_config = defaultdict(dict)
            for i, (cell, block) in enumerate(self.cb_config.keys()):
                config = raw_config[i]
                for j, etable in enumerate(etable_weights[cell].keys()):
                    etable_config[(cell, block)][etable] = config[j]
            yield etable_config

    def reduce_cb_config(self, etable_config:
                         dict[tuple[Cell, BlockType], dict[ExistentialTwoTable, int]]) \
            -> frozenset[tuple[Cell, BlockType, int]]:
        reduced_cb_config = deepcopy(self.cb_config)
        for (cell, block), config in etable_config.items():
            for etable, num in config.items():
                reduced_block = self.reduce_block_type(
                    block, etable, ab_or_ba=1)
                reduced_cb_config[(cell, block)] -= num
                reduced_cb_config[(cell, reduced_block)] += num
        return frozenset(
            (*k, v) for k, v in reduced_cb_config.items() if v > 0
        )

    def sample_and_update(self, etable_config:
                          dict[tuple[Cell, BlockType], dict[ExistentialTwoTable, int]]) \
            -> dict[ExistentialTwoTable, list[int]]:
        etable_to_elements: dict[ExistentialTwoTable,
                                 list[int]] = defaultdict(list)
        updated_elements = dict()
        # sample
        for cb_type, config in etable_config.items():
            idx = 0
            elements = list(self.cb_elements[cb_type])
            # NOTE: we need to shuffle it again!
            random.shuffle(elements)
            updated_elements[cb_type] = defaultdict(list)
            for etable, num in config.items():
                etable_to_elements[etable] += elements[idx:(idx + num)]
                reduced_block = self.reduce_block_type(
                    cb_type[1], etable, ab_or_ba=1
                )
                sampled_elements = elements[idx:(idx + num)]
                updated_elements[cb_type][reduced_block].extend(
                    sampled_elements)
                idx += num

        # update
        for cb_type, config in updated_elements.items():
            for reduced_block, elements in config.items():
                self.cb_elements[cb_type].difference_update(elements)
                self.cb_elements[
                    (cb_type[0], reduced_block)
                ].update(elements)
                num = len(elements)
                self.cb_config[cb_type] -= num
                self.cb_config[(cb_type[0], reduced_block)] += num
        logger.debug('update eu_config: {}', self.cb_elements)
        return etable_to_elements

    def __str__(self):
        s = ''
        for (cell, block), num in self.cb_config.items():
            s += 'Cell {}, {}: {}\n'.format(cell, list(block), num)
        return s

    def __repr__(self):
        return str(self)


class WFOMSContext(object):
    """
    Context for WFOMS algorithm
    """

    def __init__(self, problem):
        self.domain: set[Const] = problem.domain
        self.sentence: SC2 = problem.sentence
        self.weights: dict[Pred, tuple[RingElement, RingElement]] = problem.weights
        self.cardinality_constraint: CardinalityConstraint = problem.cardinality_constraint

        logger.info('sentence: \n{}', self.sentence)
        logger.info('domain: \n{}', self.domain)
        logger.info('weights:')
        for pred, w in self.weights.items():
            logger.info('{}: {}', pred, w)
        logger.info('cardinality constraint: {}', self.cardinality_constraint)

        # For existantial quantifiers existential quantified predicates,
        # they should be some converted tseitin predicates
        self.binary_ext_preds: list[Pred] = list()
        self.other_ext_preds: list[Pred] = list()
        # self.uni_var_indices: list[int] = list()
        # self.tseitin_preds: list[Pred] = list()
        # self.tseitin_to_extpred: dict[Pred, Pred] = dict()
        # self.tseitin_to_skolem: dict[Pred, Pred] = dict()
        self.domain_preds: list[Pred] = list()
        self.domain_to_block_type: dict[Pred, BlockType] = dict()

        self.formula: QFFormula   # skolemized sentence
        self.uni_formula: QFFormula  # universally quantified sentence
        self._build()
        logger.info('skolemized sentence: {}', self.formula)
        logger.info(
            'universally quantified sentence: {}', self.uni_formula)
        self.block_encoded_formula: QFFormula = self.uni_formula & \
            self._encode_block_types() # pyright: ignore[reportAttributeAccessIssue]
        logger.info('block encoded sentence: {}',
                    self.block_encoded_formula)

        logger.info('weights for WFOMC:')
        for pred, w in self.weights.items():
            logger.info('{}: {}', pred, w)
        # build etables
        self.etables: list[ExistentialTwoTable] = self._build_etables()

    def contain_cardinality_constraint(self) -> bool:
        return self.cardinality_constraint is not None

    def contain_existential_quantifier(self) -> bool:
        return self.sentence.contain_existential_quantifier() or \
            self.sentence.contain_counting_quantifier()

    def get_weight(self, pred: Pred) -> tuple[RingElement, RingElement]:
        default = Rational(1, 1)
        if pred in self.weights:
            return self.weights[pred]
        return (default, default)

    def decode_result(self, res: RingElement):
        if not self.contain_cardinality_constraint():
            return res
        return self.cardinality_constraint.decode_poly(res) # pyright: ignore[reportArgumentType]

    def _skolemize_one_formula(self, formula: QuantifiedFormula):
        """
        Only need to deal with \\forall X \\exists Y: f(X,Y) or \\exists X: f(X,Y)
        """
        quantified_formula = formula.quantified_formula
        quantifier_num = 1
        while(isinstance(quantified_formula, QuantifiedFormula)):
            quantified_formula = quantified_formula.quantified_formula
            quantifier_num += 1

        # always introduce auxiliary predicate
        aux_pred = new_predicate(quantifier_num, AUXILIARY_PRED_NAME)
        aux_atom = aux_pred(X, Y) if quantifier_num == 2 else aux_pred(X)
        self.formula = self.formula & (quantified_formula.equivalent(aux_atom)) # pyright: ignore[reportAttributeAccessIssue]
        self.uni_formula = self.uni_formula & (
            quantified_formula.equivalent(aux_atom)) # pyright: ignore[reportAttributeAccessIssue]
        ext_formula = aux_atom
        ext_pred = aux_pred

        if ext_pred.arity == 2:
            self.binary_ext_preds.append(ext_pred)
        else:
            self.other_ext_preds.append(ext_pred)

        if quantifier_num == 2:
            skolem_pred = new_predicate(1, SKOLEM_PRED_NAME)
            skolem_atom = skolem_pred(X)
        elif quantifier_num == 1:
            skolem_pred = new_predicate(0, SKOLEM_PRED_NAME)
            skolem_atom = skolem_pred()
        self.formula = self.formula & (skolem_atom | ~ext_formula) # type: ignore
        self.weights[skolem_pred] = (Rational(1, 1), Rational(-1, 1)) # pyright: ignore[reportPossiblyUnboundVariable]

    def _build(self):
        formula = self.sentence.uni_formula
        while(isinstance(formula, QuantifiedFormula)):
            formula = formula.quantified_formula
        self.formula = formula # pyright: ignore[reportAttributeAccessIssue]

        self.ext_formulas = self.sentence.ext_formulas
        if self.sentence.contain_counting_quantifier():
            logger.info('translate SC2 to SNF')
            if not self.contain_cardinality_constraint():
                self.cardinality_constraint = CardinalityConstraint()
            for cnt_formula in self.sentence.cnt_formulas:
                uni_formula, ext_formulas, cardinality_constraint, _ = \
                    convert_counting_formula(cnt_formula, self.domain)
                self.formula = self.formula & uni_formula
                self.ext_formulas = self.ext_formulas + ext_formulas
                self.cardinality_constraint.add_simple_constraint(*cardinality_constraint)

        if self.contain_cardinality_constraint():
            self.cardinality_constraint.build()

        self.uni_formula = self.formula
        for ext_formula in self.ext_formulas:
            self._skolemize_one_formula(ext_formula)
        if self.contain_cardinality_constraint():
            self.weights.update(
                self.cardinality_constraint.transform_weighting(
                    self.get_weight, # pyright: ignore[reportArgumentType]
                )
            )

    def _encode_block_types(self) -> QFFormula:
        # Encode block type, every block type can be seen as a set of unary facts
        ext_atoms = []
        for extp in self.binary_ext_preds:
            ext_atom = extp(X, Y)
            ext_atoms.append(ext_atom)

        evidence_sentence = top
        for flag in product(*([[True, False]] * len(ext_atoms))):
            domain_pred = new_predicate(1, EVIDOM_PRED_NAME)
            domain_atom = domain_pred(X)
            if any(flag):
                for idx, f in enumerate(flag):
                    if not f:
                        continue
                    evidence = ext_atoms[idx]
                    skolem_pred = new_predicate(1, SKOLEM_PRED_NAME)
                    self.weights[skolem_pred] = (
                        Rational(1, 1), Rational(-1, 1))
                    skolem_lit = skolem_pred(X)
                    evidence_sentence = evidence_sentence & (
                        domain_atom | skolem_lit
                    )
                    evidence_sentence = evidence_sentence & (
                        skolem_lit | (~evidence)
                    )
                block_type = BlockType(
                    pred for idx, pred in enumerate(self.binary_ext_preds) if flag[idx]
                )
            else:
                block_type = BlockType()
            self.domain_preds.append(domain_pred)
            self.domain_to_block_type[domain_pred] = block_type
        sentence = evidence_sentence & exactly_one_qf(self.domain_preds)
        return sentence # pyright: ignore[reportReturnType]

    def _build_etables(self) -> list[ExistentialTwoTable]:
        etables = list()
        n_ext_preds = len(self.binary_ext_preds)
        for i in product(*([[True, False]] * n_ext_preds)):
            for j in product(*([[True, False]] * n_ext_preds)):
                etables.append(
                    ExistentialTwoTable(i, j, tuple(self.binary_ext_preds))
                )
        return etables
