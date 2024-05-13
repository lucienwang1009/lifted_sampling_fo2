


import argparse
import logging
import logzero
from logzero import logger
from sampling_fo2.context.wfomc_context import WFOMCContext

from sampling_fo2.inference.inferencer import Inferencer
from sampling_fo2.inference.db import Database
from sampling_fo2.network.mln import MLN
from sampling_fo2.parser import parse_db, parse_mln
from sampling_fo2.problems import MLN_to_WFOMC
from sampling_fo2.utils.polynomial import round_rational
from sampling_fo2.wfomc import Algo, wfomc


def get_log_prob(mln: MLN, db: Database) -> float:
    """
    Get the probability of the MLN given the database.
    :param mln: The MLN.
    :param db: The database.
    :return: The probability.
    """
    inferencer = Inferencer(mln)
    return inferencer.log_prob(db)


def parse_args():
    parser = argparse.ArgumentParser(
        description='Inference for Markov Logic Networks',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument('--input', '-i', type=str, required=True,
                        help='mln file')
    parser.add_argument('--db', type=str, required=True,
                        help='database file')
    parser.add_argument('--prob', action='store_true', default=False,
                        help='output probability of the query')
    parser.add_argument('--debug', action='store_true', default=False)
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = parse_args()
    if args.debug:
        logzero.loglevel(logging.DEBUG)
    else:
        logzero.loglevel(logging.INFO)

    mln = parse_mln(args.input)
    db = parse_db(args.db)
    logger.debug(f'MLN: {mln}')
    logger.debug(f'Database: {db}')
    if args.prob:
        prob = get_log_prob(mln, db)
        print(prob)
    else:
        raise NotImplementedError('Only probability inference is supported')
