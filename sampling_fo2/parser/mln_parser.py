from lark import Lark
from sampling_fo2.network.mln import MLN, WeightedFormula
from sampling_fo2.fol.syntax import Const, Pred
from sampling_fo2.network.constraint import CardinalityConstraint
from sampling_fo2.parser.cardinality_constraints_parser import CCTransfomer
from sampling_fo2.parser.mln_grammar import grammar
from sampling_fo2.utils import Rational


from sampling_fo2.parser.fol_parser import FOLTransformer

from sampling_fo2.fol.syntax import *

class MLNTransformer(FOLTransformer, CCTransfomer):

    def domain_elements(self, args):
        return list(args)

    def int_domain(self, args):
        return int(args[0])

    def element(self, args):
        return args[0].value

    def set_domain(self, args):
        return set(args[0])

    def domain_name(self, args):
        return args[0].value

    def domain(self, args):
        domain_name, domain_spec = args
        if isinstance(domain_spec, int):
            domain_spec = set(
                Const(f'{domain_name}{i}') for i in range(domain_spec)
            )
        return (domain_name, domain_spec)

    def weight(self, args):
        return float(args[0])

    def rules(self, args):
        return args

    def rule(self, args):
        w, r = args[0]
        return WeightedFormula(r, w)

    def hard_rule(self, args):
        return float('inf'), args[0]

    def soft_rule(self, args):
        return args[0], args[1]

    def mln(self, args):
        weighted_formulas = args[0]
        domain = args[1][1] # Only one definition domain is supported
        cardinality_constraints = args[2]

        ccs: list[tuple[dict[Pred, float], str, float]] = list()
        if len(cardinality_constraints) > 0:
            for cc in cardinality_constraints:
                new_expr = dict()
                expr, comparator, param = cc
                for pred_name, coef in expr.items():
                    pred = self.name2pred.get(pred_name, None)
                    if not pred:
                        raise ValueError(f'Predicate {pred_name} not found')
                    new_expr[pred] = coef
                ccs.append((new_expr, comparator, param))
            cardinality_constraint = CardinalityConstraint(ccs)
        else:
            cardinality_constraint = None

        mln = MLN(weighted_formulas, domain, cardinality_constraint)
        return mln

def parse(text: str) -> MLN:
    mln_parser = Lark(grammar,
                        start='mln')
    tree = mln_parser.parse(text)
    mln = MLNTransformer().transform(tree)

    return mln
