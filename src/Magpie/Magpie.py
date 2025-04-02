# Magpie: multi-objective archive genetic programming (from Ireland)

from pathlib import Path
import random
import numpy as np
import pandas as pd
np.seterr(all='raise')
from sklearn.model_selection import train_test_split
from sklearn.base import BaseEstimator, RegressorMixin

from .grammar import Grammar, derive_string
from .exceptions import *
from .fitness import evaluate, one_m_r2, latex_eqn
from .interval import Interval, generate_bounds
from .pareto import is_pareto_efficient

class MagpieRegressor(BaseEstimator, RegressorMixin):
    def __init__(self,
                 maxevals=10000, 
                 initevals=3000, 
                 mutprob=0.8,
                 maxgenomelen=30, 
                 maxcohortlen=7,
                 gramfile="symbolic_regression.bnf",
                 valsize=0.2, 
                 initprob=0.0):
        self.maxevals = maxevals
        self.initevals = initevals
        self.mutprob = mutprob
        self.initprob = initprob # xoverprob = 1 - initprob - mutprob
        self.maxgenomelen = maxgenomelen
        self.maxcohortlen = maxcohortlen
        if Path(gramfile).exists():
            self.gramfile = gramfile
        else:
            from importlib.resources import files
            self.gramfile = str(files(__package__) / 'grammars' / gramfile)
        self.valsize = valsize
        self.n_consts = self.maxgenomelen // 2

        assert initprob + mutprob <= 1.0
        assert initevals <= maxevals

    def fit(self, X, y=None, X_bounds=None):

        n_vars = X.shape[1]
        print('n_vars', n_vars)
        self.gram = Grammar(self.gramfile, n_vars=n_vars, n_consts=self.n_consts)
        try:
            self.column_names_in_ = X.columns
            X = X.values
            y = y.values
        except:
            self.column_names_in_ = [f'X{i}' for i in range(n_vars)]



        if self.valsize > 0:
            X_train, X_val, y_train, y_val = \
                train_test_split(X, y, test_size=self.valsize)
        else:
            X_train, X_val, y_train, y_val = \
                X, None, y, None
        
        if X_bounds is None:
            self.X_bounds = None # don't use bounds
        elif type(self.X_bounds) == float:
            assert 0.0 <= X_bounds <= 1.0 # X_bounds is a margin of error
            self.X_bounds = generate_bounds(X_train, X_bounds) # add the margin of error
        else:
            assert len(X_bounds) == n_vars
            self.X_bounds = X_bounds

            
        f_cache = {} # a cache to store results for all individuals ever seen

        # make the main population structure, empty for now...
        # maxgenomelen is the max genome length, but it effectively
        # tells us the number of cohorts also because a cohort
        # is defined by the effective genome length (ie number of used codons)
        pop = EvoLengthPop(self.maxgenomelen, self.maxcohortlen)

        # main loop. at every step we make one new genome, evaluate
        # it, and add it to the population structure.  evaluation
        # means several checks plus evaluation on training data.
        # adding to the population structure really means trying to
        # add: the population can reject it. at the start of the run,
        # we generate the genome at random. later, it's by
        # mutation. (it would be easy to add crossover also.)  the
        # philosophy of the method is: we keep only pretty good
        # individuals of each length (call this a cohort). we select
        # randomly (selection is not biased by fitness). but the new
        # offspring is added to the population structure only if it's
        # good (so replacement is controlled by fitness).
        
        evals = 0
        while evals < self.maxevals:

            if evals < self.initevals:
                # we're in the early part of the run: make random genome
                g = self.random_genome()
            else:
                # we've finished the early part, so make a new
                # individual by init, mut, or xover
                r = random.random()
                if r < self.initprob:
                    g = self.random_genome()
                elif r < self.initprob + self.mutprob:
                    # take a random genome from the elite and mutate
                    ind = pop.random() # get random ind from pop. genome is at idx -1
                    uc, g = ind[0], ind[-1]
                    g = self.mutate_genome(uc, g)
                else:
                    # xover
                    ind0 = pop.random()
                    ind1 = pop.random()
                    uc0, g0 = ind0[0], ind0[-1]
                    uc1, g1 = ind1[0], ind1[-1]
                    g = self.xover_genome(uc0, g0, uc1, g1)

            try:
                # make phenotype and check if valid and non-duplicate
                ps, used_codons = derive_string(self.gram, g)
                if ps is None: raise InvalidIndividualException
                if ps in f_cache:
                    raise DuplicateIndividualException
                # store ps, because we have seen it. We should never revisit the same in future
                f_cache[ps] = used_codons, None, None, None, None, None, None

                # evaluate, optimise, simplify
                evals += 1
                res = evaluate(ps, X_train, y_train, X_val, y_val, X_bounds)
                if not(np.isfinite(res[0]) and np.isfinite(res[1])):
                    # nan or inf as train or val fitness
                    raise SingularityException
                       
                f_cache[ps] = used_codons, *res # update ps fitness value

                
                # store the individual in the pop structure
                if pop.add((used_codons, *res, g)):
                    pass 
                    #print(f"evals {evals} / {self.maxevals}")
                    #print(pop) # print pop only if there was a change
                    #print("\n\n\n\n\n\n")
                
            except (InvalidIndividualException, DuplicateIndividualException,
                    SingularityException, FailedOptimisationError,
                    FloatingPointError, ZeroDivisionError):
                # we are not using protected operators, eg AQ. instead we use raw operators like
                # division and log, and we catch exceptions here. We also check intervals using
                # Interval bounds-checking. For any of these "expected" exceptions, we just react
                # as follows: don't add the individual to our population because it's bad, just
                # continue to the next iteration of the loop.
                continue 
            if evals % 1000 == 0: print(evals)

        # TODO may need to check on self.equation_ as it is used in .predict()
        # is it ok to store self.equation_ as a df?
        self.equations_, self.equation_ = pop.extract_best_eqns(self.column_names_in_)
        return self

    def random_genome(self):
        genome = np.random.randint(self.gram.production_lcm, size=self.maxgenomelen)
        return genome

    def mutate_genome(self, uc, genome):
        genome = genome.copy()
        idx = random.randrange(uc)
        gi = genome[idx]
        while genome[idx] == gi:
            genome[idx] = random.randrange(self.gram.production_lcm)
        return genome

    def xover_genome(self, uc0, g0, uc1, g1):
        g0 = g0.copy()
        # TODO investigate use of uc (used_codons) here
        idx0 = random.randrange(len(g0) - 2)
        idx1 = random.randrange(idx0, len(g0))
        g0[idx0:idx1] = g1[idx0:idx1]
        return g0
    
    def predict(self, X):
        # item at index -3 of self.equation_ 
        # is a function which expects X data (rows, columns)
        # TODO: all these indices in the individuals are ugly!
        try:
            X = X.values # if X is a DF, extract the values TODO should handle DFs as sklearn does
        except:
            pass
        # access the equation in the right format
        # self.equation_ is a DataFrame, so we can use .at[] to get the function
        return self.equation_.at[0, "equation_fn_transpose"](X)

class EvoLengthPop:
    """The population data structure. It consists of cohorts, one for each
    integer length L, from 0 up to maxlen. Eg maxlen=30. A cohort is a
    list of individuals of good individuals of length L. When we add
    an individual of length L to the population, it's checked against
    cohort L. If it's better than any individual in the cohort, it
    replaces that individual.

    """
    
    def __init__(self, maxlen, maxcohortlen):
        assert maxlen >= 20
        assert maxcohortlen >= 1
        self.maxlen = maxlen
        # the maximum size of each cohort is given by maxcohortlen.
        # however, for cohorts above L=20, we reduce it.
        self.sizes = [maxcohortlen if L < 20 else maxcohortlen // 2
                      for L in range(maxlen + 1)]
        self.total_size = sum(self.sizes)
        self.inds = [[] for _ in range(maxlen + 1)]

    def add(self, ind):
        # add <ind> to the data structure, according to logic already
        # described. return True if we actually add
        L, fit, *rest = ind
        if L > self.maxlen: return False # ind is too long, don't use
        cohort = self.inds[L] # which cohort does this ind fit in?
        if len(cohort) < self.sizes[L]:
            cohort.append(ind) # the cohort is not yet full, so add
            return True
        else:
            cohort_fits = [ind[1] for ind in cohort] # get the train fit values of the cohort
            idx = np.argmax(cohort_fits) # get the worst element of the cohort
            if fit <= cohort[idx][1]: # (L, cost, cost_test, ...)
                cohort[idx] = ind # replace it
                return True
        return False

    def random(self):
        # get a random element of the population. notice we select the
        # cohort (length) first, then choose a random element of that.
        # this is fast, but arguably a different slower method would
        # be better.. different distribution over sizes.
        cohort = random.choice(self.inds)
        # use a loop, because no guarantee the first cohort we pick
        # will be non-empty
        while not len(cohort): 
            cohort = random.choice(self.inds)
        return random.choice(cohort)

    def simplify(self):
        # TODO: try sympy for simplification. currently not used.
        ps_simplified = simplify(p, n_vars)
        p_simplified = eval("lambda X, C: " + ps_simplified)
        def simplify(p, n_vars):
            X = sympy.symbols(f"X:{n_vars}")
            C = sympy.symbols(f"C:{n_vars}")
            s = p(X, C)
            zs = s.simplify()
            return str(zs)
        

    def extract_best_eqns(self, colnames):
        # take the Pareto Front based on validation cost
        # TODO use simplify to return nicer-looking equations.
        eqns, best = [], None
        for L in range(len(self.sizes)): # we are always *increasing* by size, so PF should be *decreasing* by MSE
            cohort = self.inds[L]
            if not len(cohort): continue # empty bin
            current = min(cohort, key=lambda x: x[2]) # (L, cost, cost_val, <various fns>, g)
            if best is None or current[2] < best[2]: # this is the best so far by validation MSE
                best = current[:] # have to copy or we will modify badly (add latex twice)
                eqns.append(current) # append only if current equation is the best yet

        # add latex equation as the last element of each equation
        # element[4] is psc, ie an equation in text with numerical constants
        eqns = [eqn + (latex_eqn(eqn[4], colnames),) for eqn in eqns] 
        best = best + (latex_eqn(best[4], colnames),)
        print(len(eqns[0]), len(best))
        columns = ['size', 'loss', 'loss_validation', 'equation_no_consts', 'equation', 
                   'equation_fn', 'consts_for_equation_fn_no_consts', 'equation_fn_transpose_no_consts', 'equation_fn_transpose', 'genome', 'latex']
        # for x, y, z in zip(eqns[0], best, columns):
        #     print(z, x, y)
        eqns = pd.DataFrame(eqns, columns=columns) 
        # print(best)
        best = pd.DataFrame([best], columns=columns)
        return eqns, best

    def __str__(self): # redo this or remove 
        s = "L trainfit valfit eqn\n"
        for L in range(len(self.sizes)):
            cohort = self.inds[L]
            if not len(cohort): continue # empty bin
            ind = min(cohort, key=lambda x: x[1]) # (L, cost, cost_test, <various fns>, g)
            s += f"{L} {ind[1]:.3f} {ind[2]:.3f} {ind[4]}\n"
        return s


if __name__ == '__main__':
    from sklearn.datasets import load_diabetes
    from sklearn.model_selection import train_test_split
    X, y = load_diabetes(return_X_y=True)
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=42)
    mr = MagpieRegressor()
    mr.fit(X_train, y_train)
    print(f"R^2 on test data: {mr.score(X_test, y_test):.2f}")
