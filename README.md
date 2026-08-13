# SBGP: Size-Bin Genetic Programming

![Photo of a magpie, by Zdeněk Macháček on Unsplash](img/logo.jpg)

Photo by <a href="https://unsplash.com/@zmachacek?utm_content=creditCopyText&utm_medium=referral&utm_source=unsplash">Zdeněk Macháček</a> on <a href="https://unsplash.com/photos/black-and-white-bird-on-grass-field-mOKHZYMhnQA?utm_content=creditCopyText&utm_medium=referral&utm_source=unsplash">Unsplash</a>

"SBGP" describes the algorithm and "SB" can also stand for *snag breac*, the Irish name for Magpie.
This project was previously called Magpie, but that name clashed with several
existing projects in neighbouring topics, including
[bloa/magpie](https://github.com/bloa/magpie) for genetic improvement and
[MAgPIE](https://www.pik-potsdam.de/en/institute/departments/activities/land-use-modelling/magpie)
for land-use modelling.

# Introduction

Genetic Programming is a form of program synthesis by evolutionary algorithms,
[invented by Richard Forsyth in 1981](https://www.emerald.com/insight/content/doi/10.1108/eb005587/full/html).
Symbolic regression means non-linear, free-form regression (also known as equation search), often
done using GP. GE (O'Neill & Ryan) is a well-known form of GP, where the GP
language is defined by a grammar and individual genomes are just
lists of ints. Operator Equalisation is a form of GP which constrains the size of individuals
in the population to match a target distribution.

This repo contains SBGP, a very early / experimental implementation of an approach 
to GP symbolic regression. Like the magpie, it borrows shiny ideas from many sources:

**SBGP = Grammatical Evolution + Caching + Steady-State + pseudo-Pareto Front (via Operator Equalisation) + Interval Arithmetic + Constant
Optimisation.**



# Algorithm

1. Create a set of bins, one per integer size, up to some maximum size. Each bin has a max capacity.
2. Create N individuals, ie GE genomes sampled uniformly, and convert by GE derivation to regression equations. For each individual, evaluate fitness and add it to the appropriate bin based on its size, measured as number of active codons used by the GE derivation.
3. While our budget is not exhausted: select an individual at random from the bins and mutate (or select two and crossover). Evaluate fitness and try-add it to the appropriate bin based on its size. (There is no fixed pop size or number of generations.) If a bin is at capacity, try-add only adds it if it's better than at least one (and then the worst is removed).
4. By the way, keep a cache of individuals seen. If step 2 or step 3 creates a duplicate at genotype or phenotype level, discard it before adding.
5. By the way, evaluation for fitness may include evaluation for invalid computations, using interval arithmetic.




# Installation

## Just to use it
```
pip install git+https://github.com/jmmcd/SBGP.git
python -m SBGP.test # check installation worked
```

## To install in order to hack on it

```
git clone https://github.com/jmmcd/SBGP.git
cd SBGP
pip install -e .
python -m SBGP.test # check installation worked
# that test is also available directly in src/SBGP/test.py
```

## Compatibility

Tested on MacOS. On Windows, equations may be not fully simplified as Sympy can hang when attempting simplification and Windows lacks a timeout mechanism to catch the hang.

# Example usage

```
from sklearn.datasets import load_diabetes
from sklearn.model_selection import train_test_split
from SBGP import SBGPRegressor
import pandas as pd
X, y = load_diabetes(return_X_y=True)
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=42)
mr = SBGPRegressor(maxevals=20000)
mr.fit(X_train, y_train)
pd.set_option('display.float_format', '{:.2f}'.format)
pd.set_option('display.max_colwidth', None)
print(mr.equations_[["size", "loss", "loss_validation", "equation"]])
print(f"The R^2 on test data of the eqn with lowest validation MSE: {mr.score(X_test, y_test):.2f}")
```

Result:

```
   size  loss  loss_validation                                                                             equation
0     2  0.41             0.33                                                                                 X[2]
1     6  0.35             0.27                                                                        (X[2] + X[8])
2    10  0.35             0.27                                                               (X[8] + (X[2] + X[2]))
3    18  0.31             0.26                                             (X[8] + ((X[2] - X[6]) + (X[2] + X[3])))
4    27  0.31             0.26             (X[8] + ((X[3] - (-150.8749 - (X[2] + log(364.6661)))) - (X[6] - X[2])))
5    30  0.35             0.26  ((((527.9763 + X[9]) + (0.5623 / X[0])) * X[8]) - ((X[2] * -765.4451) + -151.6155))
The R^2 on test data of the eqn with lowest validation MSE: 0.49
```



# Authors

James McDermott (2022-2026).

# Copyright

See LICENSE.

# Citation

A paper is in review (May 2026).

# Ongoing research and suggestions for student projects

See `TODO.md`.

