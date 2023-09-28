# Exact Sampler for FO2 sentence and SC2

This tool is for sampling instances or combinatorical structures from FO2 and SC2 theory.


## Input format

- Markov Logic Network (MLN) format, also see [Pracmln](http://www.pracmln.org/mln_syntax.html)
- Cardinality constraint: `|P| = k`


### Example input file

A `friends-smokes` MLN:

```
person = {C1, C2, C3, C4, C5, C6}
friends(person,person)
smokes(person)

3 smokes(x)
1 !friends(x,y) v !smokes(x) v smokes(y) # i.e., friends(x,y) ^ smokes(x) => smokes(y). NOTE: only support CNF for now
```

2-regular graphs:

```
vertex = 10
E(vertex, vertex)
F1(vertex, vertex)
F2(vertex, vertex)

!E(x,x).
!E(x,y) v E(y,x).
!E(x,y) v F1(x,y) v F2(x,y).
!F1(x,y) v E(x,y).
!F2(x,y) v E(x,y).
!F1(x,y) v !F2(x,y).
Exist y F1(x,y).
Exist y F2(x,y).
|E| = 20
```
> Note: You need to convert SC2 sentence into FO2 sentence with cardinality constraints by yourself.

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
$ python sampling_fo2/sampler.py -i models/friendsmoker.mln -k 10 -s
```
Find more arguments: 
```
$ python sampling_fo2/sampler.py -h
```

## References
