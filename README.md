# Magpie: multi-objective archive genetic programming (from Ireland)

![Photo of a magpie, by Zdeněk Macháček on Unsplash](img/logo.jpg)

Photo by <a href="https://unsplash.com/@zmachacek?utm_content=creditCopyText&utm_medium=referral&utm_source=unsplash">Zdeněk Macháček</a> on <a href="https://unsplash.com/photos/black-and-white-bird-on-grass-field-mOKHZYMhnQA?utm_content=creditCopyText&utm_medium=referral&utm_source=unsplash">Unsplash</a>

# Introduction

Genetic Programming is a form of program synthesis by evolutionary algorithms,
[invented by Richard Forsyth in 1981](https://www.emerald.com/insight/content/doi/10.1108/eb005587/full/html).


GE (O'Neill & Ryan) is a well-known form of GP, where the GP
language is defined by a grammar and individual genomes are just
lists of ints. 

This repo contains Magpie, a very early / experimental implementation of an approach 
to GE-like symbolic regression which borrows shiny ideas from many sources.

Magpie = Grammatical Evolution + Steady-State + pseudo-Pareto Front (via Operator Equalisation) + Interval Arithmetic + Constant
Optimisation.



# Algorithm

1. Create a set of bins, one per integer size, up to some maximum size. Each bin has a max capacity.
2. Create N individuals, ie regression equations. For each individual, evaluate fitness and add it to the appropriate bin based on its size.
3. While our budget is not exhausted:
4. Select an individual at random from the bins and mutate (or select two and crossover). Evaluate fitness and add it to the appropriate bin based on its size. (There is no fixed pop size or number of generations.)
5. When a bin is over-capacity, remove the worst from that bin based on fitness.
6. By the way, keep a cache of individuals seen. If step 2 or step 4 creates a duplicate at genotype or phenotype level, discard it before adding.
7. By the way, evaluation for fitness may include evaluation for invalid computations, using interval arithmetic.




# Installation

`pip install -e .`

# Example usage

```
from Magpie import MagpieRegressor
from sklearn.datasets import load_diabetes
from sklearn.model_selection import train_test_split
X, y = load_diabetes(return_X_y=True)
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=42)
mr = MagpieRegressor()
mr.fit(X_train, y_train)
print(f"R^2 on test data: {mr.score(X_test, y_test):.2f}")
```



# Authors

James McDermott (2022-2025)

# Copyright

See LICENSE




# Many TODOs and DONE items and IDEAS

In GE there is a genotype to phenotype mapping, which is
just grammatical derivation. The production choice at each step is
given by the next int in the genome. 

Improvements over GE are certainly possible (Rothlauf; Whigham; LTGE, by Moraglio and McDermott, unpublished) but we will not pursue them for now.

Pareto GP is due to Vladislavleva and colleagues, also other researchers. The idea is to use multi-objective optimisation to minimise error and minimise equation complexity (often measured just by tree size). NOT NEEDED.

Operator Equalisation (OE) is a neat idea due to Silva, Dignum &
Vanneschi. The idea is the population is a sequence of bins, each
bin storing individuals of certain lengths. Bins are not allowed to
go over some maximum size. In OE there is a proper population, with
selection and replacement etc. OE is not seen as a Pareto method,
because the final result is a single tree rather than a front. NOT
NEEDED.

Silva Reassembling OpEq paper: an unexpected feature of 
Operator Equalisation is revealed, one that may be the true 
responsible for its success: a nearly flat length distribution 
target.

We measure length as *effective genome length*. This is very simple and enables a quick/easy population structure.

See "A secret revealed" https://dl.acm.org/doi/abs/10.1145/2043118.2043120

Our method can be seen as a slight modification of OE as a replacement for Pareto GP. That is, we propose to use OE-like bins, one bin for every tree size up to some maximum, say 30. For each bin we store the best `maxcohortlen` individuals ever encountered. Thus the data structure is purely elitist, as is common in multi-objective algorithms. We do not have a true Pareto front, because our best individual of length 27 might be worse (in error) than our best individual of length 26. 
However, the data structure will tend to approximate a Pareto front, and at the end of the run we can quickly extract the Pareto front from the current population. DONE.

We also have the option to calculate the Pareto front at the end of the run and discard individuals not on the Pareto front. TODO.

Idea: we could have a Pareto front within each bin, seeing the test cases as objectives. Or a many-objective approach such as lexicase selection.

Insight: multi-objective algorithms can be factored into
multi-objective search per se, diversity preservation and protection
of weak individuals during the run, and return of multi-objective
front. NSGA2 does all of these using the same mechanism and then
returns the front at the end. Methods like ALPS and niching are for
diversity and protection of weak individuals. 

Our algorithm is a steady-state algorithm, as we create one
individual at a time, not a generation. DONE.

Another main idea in GP is Geometric Semantic GP and more broadly,
the use of multiple trees, with linear combination, to approximate
the target. This seems to lead to good performance, but the result
is a huge, non-interpretable set of trees. WILL NOT DO

The other main idea in GP these days is lexicase selection, or 
epsilon-lexicase in the case of regression. How can we incorporate this?

Perhaps each cohort is actually unlimited in size, but every selection
is just a lexicase selection within a randomly-chosen cohort. Perhaps
we sometimes discard ones which have not been selected for a long time,
if we wish to save memory. TODO.

As mentioned above, we have a maximum size, eg 30. This is expressed
as the number of codons in the genome, which corresponds to the
number of nodes in the derivation tree (except for the case where a
grammar LHS has only a single RHS, as in this case some GE
implementations don't use a codon). We choose this size based on:
how large a tree would we feel comfortable calling it "readable", eg
by printing the equation in our paper or showing it to a domain
expert in the area our data is collected from? DONE

A possible concern is that this maximum size might limit
evolutionary dynamics, eg we might find that even if the best
solution is less than 30 nodes, imposing this maximum makes it
harder to *find* that best solution. I think this is unlikely to be a
genuine problem. If we are concerned about this, we can just
increase the maximum size, eg to 40, to allow evolution to search in
this area of the space, but then discard all over-length individuals
(greater than 30) after the run. This could be part of a more
general post-processing to produce a Pareto front (mentioned
above). TODO.

We have an unusual initialisation also. We use a large proportion of
the run to create random individuals. Most runs create a random
initial population, but we allow the number of initial random
individuals to be much larger. They are "protected" from each other
by the population structure. If we make the number of initial
randoms small, we notice a loss of diversity early in the search, eg
all cohorts dominated by x[10] + something. DONE

After the initial random individuals, we create individuals by init,
mut or xover. Those created by init are "protected" by the
population structure. DONE

Interval arithmetic was proposed by Keijzer. Here the idea is to
avoid singularities, such as divide-by-zero, log of zero or negative
numbers, square root of negative numbers, and some others. Koza
defined protected operators which prevent this type of problem, but
they still introduce discontinuities, and are less readable in
applications. Another approach is to simply discard individuals
which give errors during training, but the problem then is that they
may give no errors during training but may still give errors on
unseen data. Keijzer's solution was to use interval arithmetic,
which calculates an upper and lower bound on every calculation in
the tree, given bounds on the variables. For example, if x0 is in
[-1, 1], then x0**2 is in [0, 1]. We may know bounds on variables
from domain knowledge, or we can derive them with a margin of error
from the training data. If a calculation has the possibility of
error, then the individual is discarded, for example if x0 is in
[-1, 1] then 1 / x0 could give divide-by-zero. DONE.

There is no point in evaluating the same genome twice. Hemberg
investigated the several forms of redundancy in GE. TODO: just
discard a duplicate genome. DONE?

Duplicate phenotypes (trees) is a huge problem in GE, and to a
lesser extent in other forms of GP. There is no point in evaluating
the same phenotype twice. More than 75% of phenotypes produced by a
typical GE grammar could be duplicates, as shown by Nicolau. We do
not design any mechanisms to avoid this, instead we simply store a
cache of previously seen phenotypes and throw away duplicates before
evaluation. DONE.

There are further de-duplication steps which could be
incorporated. Worm & Chiu reduce each tree to a canonical form, eg
both x0 + c0 and c0 + x0 are really the same. By reducing to a
canonical form, we can identify further effective duplicates and
avoid evaluating them. Sympy can do this. TODO.

Other pointless formulae can also be identified and discarded before
evaluation, eg log(c0). TODO.

PySR (Cramner) incorporates special complexity penalties which
prevent non-idiomatic functions such as sin(sin(sin(x0))). ITEA (Olivetti) achieved similar in a different way. TODO.

Sympy could also be used to simplify formulae before
evaluation. Rockett showed that simplification at the end helps more
than during, I think. TODO.

Optimisation of constants has been used by many authors. The idea is
that each equation is represented by something like X0 + sin(X1 *
C0) / C1. The best values [C0, C1] can be determined by a
pseudo-gradient algorithm, eg provided by scipy. DONE.

We use training data for optimisation of constants, and then
validation data to determine which individuals go into the
bins. This helps to favour individuals which are not over-fitted. We
use the unusual setting of train size 0.3 and validation size 0.7,
in order to make the optimisation of constants much faster, and to
strongly prevent over-fitting by the constants. DONE.

Rockett investigated both optimisation of constants during a run,
and after a run, and did not find a large advantage in favour of
using it during a run. Since it's so much slower, there is an
argument for using it only after. TODO.

Rockett also investigated and gave insight into different algorithms
for optimisation of constants, eg LM versus SLSQP. We are using
scipy fmin, which uses LBFGS. TODO.

The grammar is specified by a .bnf file. Eg we have power **, and
log, and 3 constants. These can be changed. If adding new functions,
might need to add them to fn_mappings in fitness.evaluate().


No need to discard old individuals. Just insert each ind at the
appropriate idx in its cohort. Select with a triangular distribution
with one hyperparameter. This assumes we have plenty of RAM to store every individual.

TODO: save processing by caching subtrees' semantics.

TODO: Use namedtuples instead of plain tuples to store individuals.

TODO: compare against PySR, FFX, PyOperon, QD-GP or MAP-Elites GP

TODO: Use GE_RANGE for consts as well as for vars
<constidx> ::= GE_RANGE:n_consts
C[<constidx>]

TODO: Set column_names_in_ carefully, generate x0 etc where needed, use varnames (add 'x' if starts with a number, replace everything with alphanumeric)

TODO: export final equations to a csv, to sympy, and to latex (and simplify/beautify them) compatible with PySR.


TODO: Sensitivity analysis. One approach would be, apply a 95%, 90% confidence bound on each variable using interval, ie x0 = [0.5, 0.7], and see what are the possible outputs of the model. Another approach: for each variable, choose the *mean* for every other variable, and vary this variable from min to max of its range, to see the effect on the output. Another approach: for each training point, for each variable x (holding others constant), calculate f(x) - f(x - 1SD) and calculate f(x + 1SD) - f(x), that is calculate the effect of a 1SD increase, but in two scenarios. If either x - 1SD or x + 1SD goes outside the range of y, consider discarding it? Maybe use 0.1 SD instead. Then for each variable, we can see the distribution of the effect on f of a 1 SD increase in x.

Maybe this distribution tells us something interesting about f. If it includes both pos and neg, does that tell us the model is non-linear in a useful way?


interval arithmetic example:

```
from interval import interval, inf, imath
s = "3 * x / y"
x = interval([0, 1])
y = interval([0, 1])
eval(s) # Out[9]: interval([-inf, inf])
y = interval([0.1, 1])
eval(s) # Out[11]: interval([0.0, 30.0])

inf in s_int # Out[15]: False 
```



