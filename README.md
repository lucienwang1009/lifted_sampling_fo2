# Exact Lifted Sampler for Two-Variable Logic

This tool is for sampling instances or combinatorical structures from the two-variable fragment of first-order logic.



## Input format

1. First-order sentence with at most two logic variables, see [fol_grammer.py](sampling_fo2/parser/fol_grammer.py) for details, e.g.,
  * `\forall X: (\forall Y: (R(X, Y) <-> Z(X, Y)))`
  * `\forall X: (\exists Y: (R(X, Y)))`
  * `\exists X: (F(X) -> \forall Y: (R(X, Y)))`
  * ..., even more complex sentence...
2. Domain: 
  * `domain=3` or
  * `domain={p1, p2, p3}`
3. Weighting (optional): `positve_weight negative_weight predicate`
4. Cardinality constraint (optional): 
  * `|P| = k`
  * `|P| > k`
  * `|P| >= k`
  * `|P| < k`
  * `|P| <= k`
  * ...


### Example input file

2 colored graphs:
```
\forall X: (\forall Y: ((E(X,Y) -> E(Y,X)) &
                        (R(X) | B(X)) &
                        (~R(X) | ~B(X)) &
                        (E(X,Y) -> ~(R(X) & R(Y)) & ~(B(X) & B(Y)))))

V = 10
```


2 regular graphs:
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
> Note: You need to convert SC2 sentence into FO2 sentence with cardinality constraints by yourself.

Sampling possible worlds from `friends-smokes` MLN:
```
\forall X: (~fr(X,X)) &
\forall X: (\forall Y: (fr(X,Y) -> fr(Y,X))) &
\forall X: (\forall Y: (aux(X,Y) <-> (fr(X,Y) & sm(X) -> sm(Y)))) &
\forall X: (\exists Y: (fr(X,Y)))

person = 10
2.7 1 aux
```

More examples are in [models](models/)


### Installation
Install the package:
```
$ pip install -e .
```


### How to use
Run the following command:
```
$ python sampling_fo2/sampler.py -i models/friendsmoker.mln -k 10 -s
```
Find more arguments: 
```
$ python sampling_fo2/sampler.py -h
```

## References

```
@inproceedings{DBLP:conf/lics/WangP0K23,
  author       = {Yuanhong Wang and
                  Juhua Pu and
                  Yuyi Wang and
                  Ondrej Kuzelka},
  title        = {On Exact Sampling in the Two-Variable Fragment of First-Order Logic},
  booktitle    = {{LICS}},
  pages        = {1--13},
  year         = {2023},
  url          = {https://doi.org/10.1109/LICS56636.2023.10175742},
  doi          = {10.1109/LICS56636.2023.10175742},
  timestamp    = {Thu, 20 Jul 2023 11:32:59 +0200},
  biburl       = {https://dblp.org/rec/conf/lics/WangP0K23.bib},
  bibsource    = {dblp computer science bibliography, https://dblp.org}
}
```
