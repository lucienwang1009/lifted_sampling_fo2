import argparse
from copy import deepcopy
import logging

from wfomc import parse_input, Const, AtomicFormula, Pred, wfomc, WFOMCProblem, round_rational, to_sc2
from wfomc.parser.fol_parser import parse
from collections import Counter
import matplotlib.pyplot as plt

from sampling_fo2.sampler import Sampler
from sampling_fo2.context import WFOMSContext

def parse_args():
    parser = argparse.ArgumentParser(
        description='Perform weighted sampling for FO2 with counting quantifiers.'
    )
    parser.add_argument(
        '--input_file',
        '-i',
        type=str,
        help='Path to the input WFOMS file.'
    )
    parser.add_argument(
        '--n_samples',
        type=int,
        help='Number of samples to generate.'
    )
    parser.add_argument(
        '--domain-size',
        type=int,
        required=True,
        help='Size of the domain to sample from.'
    )
    parser.add_argument(
        '--stats_pred',
        type=str,
        required=True,
        help='Predicate name to compute statistics on.'
    )
    return parser.parse_args()


def statistics(samples: list[set[AtomicFormula]], pred_name: str) -> None:
    stats_counter = Counter()
    for sample in samples:
        cnt = len([atom for atom in sample if atom.pred.name == pred_name])
        stats_counter[cnt] += 1
    # plot histogram
    x = list(stats_counter.keys())
    y = [stats_counter[k] for k in x]
    plt.bar(x, y)
    plt.xlabel(f'Number of true atoms for predicate {pred_name}')
    plt.ylabel('Frequency')
    plt.title(f'Statistics for predicate {pred_name}')
    plt.show()


def asymptotic_statistics(problem: WFOMCProblem, pred_name: str, max_domain_size: int) -> None:
    probs = []
    for domain_size in range(1, max_domain_size + 1, 10):
        problem.domain = set(Const(f'{i}') for i in range(1, domain_size + 1))
        context = WFOMSContext(problem)
        total_w = wfomc(problem)
        print(problem)
        conditional_problem = deepcopy(problem)
        conditional_problem.sentence.uni_formula = to_sc2(conditional_problem.sentence.uni_formula & \
            parse(f'\\forall X: ({pred_name}(X))')).uni_formula
        print(conditional_problem)
        true_w = wfomc(conditional_problem)
        probs.append(round_rational(true_w / total_w))
    print(probs)


def main():
    args = parse_args()
    problem = parse_input(args.input_file)
    asymptotic_statistics(problem, args.stats_pred, args.domain_size)
    context = WFOMSContext(problem)
    sampler = Sampler(context)
    samples = sampler.sample(args.n_samples)
    statistics(samples, args.stats_pred)


if __name__ == '__main__':
    main()
