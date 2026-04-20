# tests/test_model_satisfaction.py
from pathlib import Path
import pytest
from wfomc.parser.wfomcs_grammar import grammar
from wfomc.parser.wfomcs_parser import WFOMSTransformer
from wfomc.fol.syntax import (
    Formula, Pred, Const, Top, Bot, QFFormula,
    Var, QuantifiedFormula, Universal, Existential, Counting,
    Conjunction, Disjunction, Negation, Implication, Equivalence,
)
from wfomc.problems import WFOMCProblem
from lark import Lark
from sampling_fo2.sampler import Sampler
from sampling_fo2.context import WFOMSContext


class _RawFormulaTransformer(WFOMSTransformer):
    def __init__(self):
        super().__init__()
        self.raw_formula: Formula | None = None

    def wfomcs(self, args):
        self.raw_formula = args[0]
        return super().wfomcs(args)


def parse_wfomcs_raw(filepath: str) -> tuple[Formula, WFOMCProblem]:
    """Return (raw_formula, WFOMCProblem) for a .wfomcs file."""
    text = Path(filepath).read_text()
    lark_parser = Lark(grammar, start='wfomcs')
    tree = lark_parser.parse(text)
    transformer = _RawFormulaTransformer()
    sentence, domain, weightings, cc, unary_evidence = transformer.transform(tree)
    pred_weightings = {k: v for k, v in {sentence.pred_by_name(p): w for p, w in weightings.items()}.items() if k is not None}
    problem = WFOMCProblem(sentence, domain, pred_weightings, cc, unary_evidence)
    return transformer.raw_formula, problem


def test_parse_raw_returns_formula():
    models_dir = Path(__file__).parent.parent / 'models'
    raw, problem = parse_wfomcs_raw(str(models_dir / '2-colored-graph.wfomcs'))
    assert raw is not None
    assert len(problem.domain) == 10



def eval_ground_qf(formula, true_atoms: set) -> bool:
    """Evaluate a ground (variable-free) QFFormula against a set of true atoms."""
    import sympy
    if isinstance(formula, Top):
        return True
    if isinstance(formula, Bot):
        return False
    pos_atoms = formula.atoms()  # always returns positive AtomicFormulas
    subst = {atom: (atom in true_atoms) for atom in pos_atoms}
    result = formula.sub_nullary_atoms(subst)
    return bool(sympy.simplify(result.expr))


def test_eval_ground_qf_basic():
    P = Pred('P', 1)
    a = Const('a')
    atom = P(a)          # P(a)
    neg_atom = ~P(a)     # ~P(a)

    assert eval_ground_qf(atom, {atom}) is True
    assert eval_ground_qf(atom, set()) is False
    assert eval_ground_qf(neg_atom, set()) is True
    assert eval_ground_qf(neg_atom, {atom}) is False

    Q = Pred('Q', 1)
    b = Const('b')
    q_b = Q(b)
    conj = atom & q_b
    assert eval_ground_qf(conj, {atom, q_b}) is True
    assert eval_ground_qf(conj, {atom}) is False


def eval_formula(formula, domain: set, true_atoms: set, assignment: dict = None) -> bool:
    """Recursively evaluate a FOL formula against a Herbrand interpretation."""
    if assignment is None:
        assignment = {}

    if isinstance(formula, Top):
        return True
    if isinstance(formula, Bot):
        return False
    if isinstance(formula, QFFormula):
        grounded = formula.substitute(assignment)
        return eval_ground_qf(grounded, true_atoms)
    if isinstance(formula, QuantifiedFormula):
        q = formula.quantifier_scope
        var = formula.quantified_var
        body = formula.quantified_formula
        if isinstance(q, Universal):
            return all(
                eval_formula(body, domain, true_atoms, {**assignment, var: c})
                for c in domain
            )
        elif isinstance(q, Existential):
            return any(
                eval_formula(body, domain, true_atoms, {**assignment, var: c})
                for c in domain
            )
        elif isinstance(q, Counting):
            count = sum(
                1 for c in domain
                if eval_formula(body, domain, true_atoms, {**assignment, var: c})
            )
            cmp = q.comparator
            k = q.count_param
            if cmp == '=':    return count == k
            if cmp == '!=':   return count != k
            if cmp == '<':    return count < k
            if cmp == '>':    return count > k
            if cmp == '<=':   return count <= k
            if cmp == '>=':   return count >= k
            if cmp == 'mod':
                r, m = k
                return count % m == r
    if isinstance(formula, Conjunction):
        return (eval_formula(formula.left_formula, domain, true_atoms, assignment)
                and eval_formula(formula.right_formula, domain, true_atoms, assignment))
    if isinstance(formula, Disjunction):
        return (eval_formula(formula.left_formula, domain, true_atoms, assignment)
                or eval_formula(formula.right_formula, domain, true_atoms, assignment))
    if isinstance(formula, Negation):
        return not eval_formula(formula.sub_formula, domain, true_atoms, assignment)
    if isinstance(formula, Implication):
        return (not eval_formula(formula.left_formula, domain, true_atoms, assignment)
                or eval_formula(formula.right_formula, domain, true_atoms, assignment))
    if isinstance(formula, Equivalence):
        l = eval_formula(formula.left_formula, domain, true_atoms, assignment)
        r = eval_formula(formula.right_formula, domain, true_atoms, assignment)
        return l == r
    raise TypeError(f"Unsupported formula type: {type(formula)}")


def test_eval_formula_universal():
    P = Pred('TestP', 1)
    X = Var('X')
    domain = {Const('a1'), Const('b1')}
    formula = QuantifiedFormula(Universal(X), P(X))
    a1, b1 = Const('a1'), Const('b1')
    assert eval_formula(formula, domain, {P(a1), P(b1)}) is True
    assert eval_formula(formula, domain, {P(a1)}) is False


def test_eval_formula_existential():
    P = Pred('TestQ', 1)
    X = Var('X')
    domain = {Const('c1'), Const('d1')}
    formula = QuantifiedFormula(Existential(X), P(X))
    c1 = Const('c1')
    assert eval_formula(formula, domain, {P(c1)}) is True
    assert eval_formula(formula, domain, set()) is False


def test_eval_formula_counting():
    P = Pred('TestR', 1)
    X = Var('X')
    domain = {Const('e1'), Const('f1'), Const('g1')}
    formula = QuantifiedFormula(Counting(X, '=', 2), P(X))
    e1, f1, g1 = Const('e1'), Const('f1'), Const('g1')
    assert eval_formula(formula, domain, {P(e1), P(f1)}) is True
    assert eval_formula(formula, domain, {P(e1)}) is False
    assert eval_formula(formula, domain, {P(e1), P(f1), P(g1)}) is False


def _user_preds(formula) -> set:
    """Collect all predicates appearing in formula (excludes auxiliary @-prefixed ones)."""
    if isinstance(formula, QFFormula):
        return {a.pred for a in formula.atoms()
                if not a.pred.name.startswith('@')}
    if isinstance(formula, QuantifiedFormula):
        return _user_preds(formula.quantified_formula)
    if isinstance(formula, (Conjunction, Disjunction, Implication, Equivalence)):
        return _user_preds(formula.left_formula) | _user_preds(formula.right_formula)
    if isinstance(formula, Negation):
        return _user_preds(formula.sub_formula)
    return set()


def _check_cardinality(cc, sample: set) -> bool:
    """Check CardinalityConstraint against sample (set of positive AtomicFormula)."""
    if cc is None or cc.empty():
        return True
    for expr, cmp, k in cc.constraints:
        count = sum(
            coef * sum(1 for atom in sample if atom.positive and atom.pred == pred)
            for pred, coef in expr.items()
        )
        if cmp == '=':    ok = count == k
        elif cmp == '!=': ok = count != k
        elif cmp == '<':  ok = count < k
        elif cmp == '>':  ok = count > k
        elif cmp == '<=': ok = count <= k
        elif cmp == '>=': ok = count >= k
        else:             ok = True
        if not ok:
            return False
    return True


def verify_sample(raw_formula, domain: set, sample: set, cc=None) -> bool:
    """Return True iff sample satisfies raw_formula and optional cardinality constraint."""
    user_preds = _user_preds(raw_formula)
    filtered = {atom for atom in sample if atom.pred in user_preds}
    if not eval_formula(raw_formula, domain, filtered):
        return False
    return _check_cardinality(cc, filtered)


def test_verify_sample_colored_graph():
    """2-node colored graph: check formula E(X,Y)->~(B(X)&R(Y))&~(B(Y)&R(X))."""
    E = Pred('E', 2)
    R = Pred('R', 1)
    B = Pred('B', 1)
    a, b = Const('a'), Const('b')
    domain = {a, b}
    models_dir = Path(__file__).parent.parent / 'models'
    raw, problem = parse_wfomcs_raw(str(models_dir / '2-colored-graph.wfomcs'))
    # valid: R(a), R(b), no edges — same color, no edges satisfies all constraints
    valid_sample = {R(a), R(b)}
    assert verify_sample(raw, domain, valid_sample) is True
    # invalid: one node with both R and B — violates ~R(X) | ~B(X)
    invalid_sample = {R(a), B(a), R(b)}  # a has both colors
    assert verify_sample(raw, domain, invalid_sample) is False


_MODELS_DIR = Path(__file__).parent.parent / 'models'
_WFOMCS_FILES = sorted(_MODELS_DIR.glob('*.wfomcs'))

_SKIP_MODELS = {
    'linear_order_perm.wfomcs',  # LEQ linear order (special semantics)
    'head-middle-tail.wfomcs',   # LEQ linear order
}

N_SAMPLES = 10


@pytest.mark.parametrize(
    'model_file',
    [str(f) for f in _WFOMCS_FILES if f.name not in _SKIP_MODELS],
    ids=[f.name for f in _WFOMCS_FILES if f.name not in _SKIP_MODELS],
)
def test_samples_satisfy_sentence(model_file: str):
    raw_formula, problem = parse_wfomcs_raw(model_file)
    context = WFOMSContext(problem)
    sampler = Sampler(context)
    samples = sampler.sample(N_SAMPLES)
    failures = [
        i for i, s in enumerate(samples)
        if not verify_sample(raw_formula, problem.domain, s, problem.cardinality_constraint)
    ]
    assert failures == [], (
        f"{len(failures)}/{N_SAMPLES} samples violated the sentence in {Path(model_file).name}. "
        f"First failing sample index: {failures[0]}, "
        f"sample: {sorted(str(a) for a in samples[failures[0]])}"
    )
