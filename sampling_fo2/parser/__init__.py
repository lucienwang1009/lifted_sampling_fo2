from typing import Union

from sampling_fo2.problems import WFOMCSProblem

from .wfomcs_parser import parse as wfomcs_parse
from .mln_parser import parse as mln_parse


def parse_input(input_file: str) -> WFOMCSProblem:
    if input_file.endswith('.mln'):
        raise NotImplementedError
    elif input_file.endswith('.wfomcs'):
        with open(input_file, 'r') as f:
            input_content = f.read()
        return wfomcs_parse(input_content)
    else:
        raise RuntimeError(f'Unknown input file type: {input_file}')
