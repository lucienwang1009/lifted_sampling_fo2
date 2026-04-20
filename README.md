# Exact Lifted Sampler for Two-Variable Logic

This tool is for sampling instances or combinatorial structures from the two-variable fragment of first-order logic.

### Installation

Install UV via:
[github](https://github.com/astral-sh/uv) or
```
pip install uv
```

Sync the dependencies:
```
uv sync
```

### How to use
```
$ uv run sampler -i [input] -k [N]
```
where
- `input` is the input file with the suffix `.wfomcs` or `.mln`
- `N` is the number of samples to generate

Find more arguments:
```
$ uv run sampler -h
```

## Input format

The input file with the suffix `.wfomcs` contains the following information **in order**:
1. First-order sentence with at most two logic variables (must be in capital letters, e.g., `X`, `Y`, `Z`, etc.), see [fol_grammar.py](src/sampling_fo2/parser/fol_grammar.py) for details, e.g.,
  * `\forall X: (\forall Y: (R(X, Y) <-> Z(X, Y)))`
  * `\forall X: (\exists Y: (R(X, Y)))`
  * `\exists X: (F(X) -> \forall Y: (R(X, Y)))`
  * ..., even more complex sentence...
2. Domain:
  * `domain=3` or
  * `domain={p1, p2, p3}`, where `p1`, `p2`, `p3` are the constants in the domain (must start with a lowercase letter).
3. Weighting (optional): `positive_weight negative_weight predicate`
4. Cardinality constraint (optional):
  * `|P| = k`
  * `|P| > k`
  * `|P| >= k`
  * `|P| < k`
  * `|P| <= k`
5. Unary evidence (optional):
  * `P(p1), ~P(p3)`

### Example input file

- 2 colored graphs:
```
\forall X: (\forall Y: ((E(X,Y) -> E(Y,X)) &
                        (R(X) | B(X)) &
                        (~R(X) | ~B(X)) &
                        (E(X,Y) -> ~(R(X) & R(Y)) & ~(B(X) & B(Y)))))

V = 10
```

- 2 regular graphs:
```
\forall X: (~E(X,X)) &
\forall X: (\forall Y: ((E(X,Y) -> E(Y,X)) &
                        (E(X,Y) <-> (F1(X,Y) | F2(X,Y))) &
                        (~F1(X, Y) | ~F2(X,Y)))) &
\forall X: (\exists Y: (F1(X,Y))) & 
\forall X: (\exists Y: (F2(X,Y)))

V = 6
|E| = 12
```

- 2 regular graphs using counting quantifier (`\exists_{=2} Y: (E(X,Y))` means there are exactly 2 edges from each node):
```
\forall X: (~E(X,X)) &
\forall X: (\forall Y: (E(X,Y) -> E(Y,X))) &
\forall X: (\exists_{=2} Y: (E(X,Y)))

V = 6
```

- Sampling possible worlds from `friends-smokes` MLN:
```
\forall X: (~fr(X,X)) &
\forall X: (\forall Y: (fr(X,Y) -> fr(Y,X))) &
\forall X: (\forall Y: (aux(X,Y) <-> (fr(X,Y) & sm(X) -> sm(Y)))) &
\forall X: (\exists Y: (fr(X,Y)))

person = 10
2.7 1 aux
```

> **Note: You can also directly input the MLN in the form defined in [mln_grammar.py](src/sampling_fo2/parser/mln_grammar.py)**
```
~friends(X,X).
friends(X,Y) -> friends(Y,X).
2.7 friends(X,Y) & smokes(X) -> smokes(Y)
\forall X: (\exists Y: (friends(X,Y))).

person = 10
```

> Add unary evidence:
```
~friends(X,X).
friends(X,Y) -> friends(Y,X).
2.7 friends(X,Y) & smokes(X) -> smokes(Y)
\forall X: (\exists Y: (friends(X,Y))).

person = {alice, bob, charlie, david, eve}

smokes(alice), ~smokes(bob)
```

More examples are in [models](models/)

## References

```
@article{DBLP:journals/ai/WangPWK24,
  author       = {Yuanhong Wang and
                  Juhua Pu and
                  Yuyi Wang and
                  Ondrej Kuzelka},
  title        = {Lifted algorithms for symmetric weighted first-order model sampling},
  journal      = {Artif. Intell.},
  volume       = {331},
  pages        = {104114},
  year         = {2024},
  url          = {https://doi.org/10.1016/j.artint.2024.104114},
  doi          = {10.1016/J.ARTINT.2024.104114},
  timestamp    = {Fri, 31 May 2024 21:06:28 +0200},
  biburl       = {https://dblp.org/rec/journals/ai/WangPWK24.bib},
  bibsource    = {dblp computer science bibliography, https://dblp.org}
}
```
