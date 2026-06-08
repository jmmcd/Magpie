# Suggestions for student projects

* **Hybrid model selection** is a set of methods invented by Ramlan et al (https://link.springer.com/chapter/10.1007/978-3-032-23005-8_12). The main idea is: we have a training dataset (X, y) on which to evaluate our models' fitness at model selection time, but nothing is stopping us from looking at other X points, interpolated or slightly extrapolated outside the training data. We don't know the correct y for those new X points, so we don't know *exactly* how our model should behave. But we know *something*. It should not divide by zero - this is already captured in Magpie by interval arithmetic. It should not return enormous values, far outside the range of y observed in the training data, and should not be too sensitive to small changes in X. This part is captured by three model selection metrics proposed by Ramlan et al. It would be interesting to add these to Magpie and evaluate. The experimental design could be something like: compare the recommended model selection methods and lambda hyperparameters from the Ramlan et al paper, versus interval arithmetic, versus neither, versus both. That would be a 2x2 design.

* **Better simplification**: use a fast simplification algorithm (https://icml.cc/virtual/2026/poster/64015) and investigate whether it's better to use simplification on every equation during the run (and at which point in the evaluation cycle), or only on the final Pareto front, or something in between.

* **Assortative mating** is a concept in evolutionary biology: individuals are more likely to mate with other individuals which are *different* from them. That is, opposites attract. Does assortative mating help in Magpie? Or does the opposite help? We could measure "different" in several ways:
  * Similar **size** (so, crossover works by choose a random individual, then choose another in a similar-sized bin)
  * Similar **semantics**

  For each of these, we could try:
  
  * Assortative mating, ie bias crossover to make the second parent highly different
  * The opposite, ie bias crossover to make the second parent highly similar
  * The default, ie no such bias.

  This gives a 2 x 3 experimental design. Run the whole thing on 10-20 datasets with say 10 runs each. For all other hyperparameters, use defaults. Read the Magpie paper for background.

* The **Hypervolume** of a Pareto front is the volume of objective space covered by the set of individuals of the front. It's a standard measure of how well a multi-objective algorithm has performed. We should add it to our calculations at the end of a run, and test the implementation carefully, including thinking about the choice of reference point. Then we should run a set of experiments, starting from the Magpie reference config (default hyperparameters) and trying out changes "one factor at a time" (OFAT) to see if any improve the reported hypervolume. 


* **Explicit constants versus a placeholder versus multiple placeholders**. We currently have numbered constant placeholders `C[0]`, `C[1]`, etc up to a max determined by `prop_consts`. We set their values by numerical optimisation. We also have a `n_num_optimisations=0` setting, where there are instead explicit constants in the grammar. But maybe these multiple placeholders are a mistake. The effect is that sometimes the same constant `C[i]` will appear in more than one place in the equation, and of course all occurrences of same constant must have the same value after optimisation. Pros: simpler equations in the first place (from an information theory point of view) and more opportunities for simplification. Cons: complex code and a hyperparameter to tune (prop_consts). Could just allow `C` to occur in the grammar, and then convert every `C` to `C[i]` when walking the string. This would give a simple experimental design with just two setups: a single placeholder in the grammar versus multiple. However, it might be important to sweep the `n_num_optimisations` parameter at the same time as it could have different effects in the two cases.



# Many other TODOs

* Improvements over plain GE are certainly possible (Rothlauf; Whigham; LTGE, by Moraglio and McDermott, unpublished). Or a GP-tree approach.

* We could have a Pareto front within each bin, seeing the test cases as objectives. Or a many-objective approach such as lexicase selection.

* One main idea in GP in recent years is lexicase selection, or epsilon-lexicase in the case of regression. How can we incorporate this? Perhaps each cohort is actually unlimited in size, but every selection is just a lexicase selection within a randomly-chosen cohort. Perhaps we sometimes discard ones which have not been selected for a long time, if we wish to save memory.

* We have a limit on individual size. 
A possible concern is that this maximum size might limit
evolutionary dynamics, eg we might find that even if the best
solution is less than 30 nodes, imposing this maximum makes it
harder to *find* that best solution. I think this is unlikely to be a
genuine problem. If we are concerned about this, we can just
increase the maximum size, eg to 40, to allow evolution to search in
this area of the space, but then discard all over-length individuals
(greater than 30) after the run. This could be part of a more
general post-processing to produce a Pareto front (mentioned
above).

* We have an unusual initialisation also. We use a large proportion of
the run to create random individuals. Most runs create a random
initial population, but we allow the number of initial random
individuals to be much larger. They are "protected" from each other
by the population structure. If we make the number of initial
randoms small, we notice a loss of diversity early in the search, eg
all cohorts dominated by x[10] + something.

* We already de-duplicate using a cache, There are further de-duplication steps which could be
incorporated. Worm & Chiu reduce each tree to a canonical form, eg
both x0 + c0 and c0 + x0 are really the same. By reducing to a
canonical form, we can identify further effective duplicates and
avoid evaluating them. Rockett showed that simplification at the end helps more
than during, I think. Sympy can do this. DONE but not tested yet

* Other pointless formulae should also be identified and discarded before
evaluation, eg log(c0). TODO.

* PySR (Cramner) incorporates special complexity penalties which
prevent non-idiomatic functions such as sin(sin(sin(x0))). ITEA (Olivetti) achieved similar in a different way. TODO.

* Rockett investigated both optimisation of constants during a run,
and after a run, and did not find a large advantage in favour of
using it during a run. Since it's so much slower, there is an
argument for using it only after. TODO.

* Rockett also investigated and gave insight into different algorithms
for optimisation of constants, eg LM versus SLSQP. We are using
scipy `curve_fit`, which uses partial least squares. We could investigate this more.

* No need to discard old individuals. Just insert each ind at the
appropriate idx in its cohort. Select with a triangular distribution
with one hyperparameter. This assumes we have plenty of RAM to store every individual.

* TODO: save processing by caching subtrees' semantics.

* TODO: compare against PySR (Done), FFX (Done), PyOperon, QD-GP or MAP-Elites GP

* TODO: Set column_names_in_ carefully, generate x0 etc where needed, use varnames (add 'x' if starts with a number, replace everything with alphanumeric)

* TODO: export final equations to a csv, to sympy, and to latex (and simplify/beautify them) compatible with PySR. Some of this is DONE

* TODO: Sensitivity analysis. One approach would be, apply a 95%, 90% confidence bound on each variable using interval, ie x0 = [0.5, 0.7], and see what are the possible outputs of the model. Another approach: for each variable, choose the *mean* for every other variable, and vary this variable from min to max of its range, to see the effect on the output. Another approach: for each training point, for each variable x (holding others constant), calculate f(x) - f(x - 1SD) and calculate f(x + 1SD) - f(x), that is calculate the effect of a 1SD increase, but in two scenarios. If either x - 1SD or x + 1SD goes outside the range of y, consider discarding it? Maybe use 0.1 SD instead. Then for each variable, we can see the distribution of the effect on f of a 1 SD increase in x. Maybe this distribution tells us something interesting about f. If it includes both pos and neg, does that tell us the model is non-linear in a useful way?


* Threshold on accuracy to avoid too-simple models contributing too often to crossover/mutation

* How to visualise a large factorial experiment? 

* Need to rethink validation data for constant optimisation and for model selection. Maybe it is ok, but maybe re-unify
them for hypervolume calculation? Maybe use *both* for model selection?

* Have implemented hypervolume and simplification but need to test both.

* Simplification before fitting constants would be nice as it would save curve_fit time, but it requires some hacking to deal with (c[0] + c[1]) because sympy doesn't know they are constants that can be absorbed. Simplification before printing the final equation would be easier because the constants are now constants. Simplification before the final *latex* is what we have right now.

* Consider counting tree complexity using sympy `len(list(sp.preorder_traversal(expr)))` either during the run, or after simplification at the end.

* Consider changing to MSE or R^2 instead of ugly 1-R^2.

* We should seed the initial population with some fitted equations:
  * Constant
  * Single linear regressions of each single variable
  * Multiple linear regressions of some pairs of triples of variables
  * Same for polynomials up to degree 4
  * Same for rational functions up to degree 4
  * Same for some log and exp templates
  * (Maybe just run FFX and use it to seed the population?)




