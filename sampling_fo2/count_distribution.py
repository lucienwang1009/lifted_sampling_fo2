import argparse
import logging
import logzero

from sampling_fo2.context import WFOMCContext
from sampling_fo2.parser import mln_parse
from sampling_fo2.problems import MLN_to_WFOMC
from sampling_fo2.wfomc import Algo
from sampling_fo2.wfomc import count_distribution


def parse_args():
    parser = argparse.ArgumentParser(
        description='Compute the count distribution of a given MLN',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('--input', '-i', type=str, required=True,
                        help='mln file')
    parser.add_argument('--debug', action='store_true', default=False)
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = parse_args()
    if args.debug:
        logzero.loglevel(logging.DEBUG)
    else:
        logzero.loglevel(logging.INFO)
    input_file = args.input
    if not input_file.endswith('.mln'):
        raise RuntimeError(f'Only support .mln file: {input_file}')

    with open(input_file, 'r') as f:
        input_content = f.read()
    mln = mln_parse(input_content)
    preds = mln.preds()
    wfomcs_problem = MLN_to_WFOMC(mln)
    context = WFOMCContext(wfomcs_problem)
    count_dist = count_distribution(
        context, preds, Algo.FASTERv2
    )
    print(count_dist)
