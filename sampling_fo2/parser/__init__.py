from sampling_fo2.network.mln import MLN
from sampling_fo2.problems import WFOMCSProblem, MLN_to_WFOMC
from .wfomcs_parser import parse as wfomcs_parse
from .mln_parser import parse as mln_parse
from .db_parser import DBParser


def parse_mln(input_file: str) -> MLN:
    with open(input_file, 'r') as f:
        input_content = f.read()
    return mln_parse(input_content)


def parse_wfomcs(input_file: str) -> WFOMCSProblem:
    with open(input_file, 'r') as f:
        input_content = f.read()
    return wfomcs_parse(input_content)


def parse_input(input_file: str) -> WFOMCSProblem:
    if input_file.endswith('.mln'):
        mln = parse_mln(input_file)
        wfomcs_problem = MLN_to_WFOMC(mln)
        return wfomcs_problem
    elif input_file.endswith('.wfomcs'):
        return parse_wfomcs(input_file)
    else:
        raise RuntimeError(f'Unknown input file type: {input_file}')


def parse_db(db_file: str) -> set:
    with open(db_file, 'r') as f:
        db_content = f.read()
    return DBParser().parse(db_content)
