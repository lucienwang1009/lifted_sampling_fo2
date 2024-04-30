from __future__ import annotations

from sampling_fo2.fol.syntax import Pred
from sampling_fo2.inference.db import Database


class DBParser(object):
    """
    A parser for database files.
    The database file format is that one ground atom per line.
    """
    def parse(self, text: str):
        """
        Parse a database file.
        :param text: The text of the database file.
        :return: A set of ground atoms.
        """
        text = text.strip().split('\n')
        facts = set()
        for line in text:
            pred_name, args = line.split('(')
            args = args[:-1].split(',')
            pred = Pred(pred_name, len(args))
            facts.add(pred(*args))
        return Database(facts)


if __name__ == '__main__':
    db_parser = DBParser()
    db = db_parser.parse('A(a,b)\nB(a)')
    print(db)
    # Expected output:
    # {A(a, b), B(a)}
