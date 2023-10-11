# Exact Lifted Sampler for Two-Variable Logic

This tool is for sampling instances or combinatorical structures from the two-variable fragment of first-order logic.


## Input format

1. First-order sentence with at most two logic variables, see [sampling_fo2/parser/fol_grammer.py](fol_grammar.py) for details, e.g.,
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

2 regular graphs (you need to convert C2 to FO2 + cardinality constraint):
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
Install requirements:
```
$ pip install -r requirements.txt
```
Add path to your PYTHONPATH:
```
$ export PYTHONPATH=$(pwd)/sampling_fo2:$PYTHONPATH
```


### How to use
Run the following command:
```
$ python sampling_ufo2/sampler.py -i models/friendsmoker.mln -k 10
```
Find more arguments: 
```
$ python sampling_ufo2/sampler.py -h
```

## References
