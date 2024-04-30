from __future__ import annotations

from sampling_fo2.fol.syntax import AtomicFormula


class Database(object):
    def __init__(self, facts: frozenset[AtomicFormula]):
        """
        Database class, i.e., a set of facts.
        """
        self.facts = facts
        self.preds = set()
        self.domain = set()
        for fact in facts:
            self.preds.add(fact.pred)
            self.domain.update(fact.args)

    def __str__(self):
        return str(self.facts)
