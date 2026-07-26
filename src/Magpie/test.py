print("Running a very quick test to check installation was ok.\n")

import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split
from Magpie import MagpieRegressor

rng = np.random.default_rng(42)

n = 200
v0 = rng.uniform(-1, 1, n)
v1 = rng.uniform(0.01, 1, n)  # positive to keep log(0.5 * v1) defined
noise = rng.normal(0, 0.05, n)

# use a DataFrame/Series with custom names so we can check that Magpie
# picks up the column names instead of falling back to X0, X1, ...
X = pd.DataFrame({"v0": v0, "zoo": v1})
y = pd.Series(v0**3 + np.log(0.5 * v1) + noise, name="y")

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=0)

# 3000 evals (still ~2s) so that the Pareto front contains equations which
# actually use the variables, ie where the names are visible
mr = MagpieRegressor(maxevals=3000, init_prop=0.2,
                     mutprob=0.5,
                     gramfile="base_consts.bnf",
                     prop_consts=0.0,
                     n_num_optimisations=0,
                     random_state=42)
mr.fit(X_train, y_train)

eqs = mr.equations_
print(f"Column names seen by Magpie: {list(mr.column_names_in_)}")
print(f"Found {len(eqs)} equations on the Pareto front\n")

for i, row in eqs.iterrows():
    print(f"--- Equation {i} (size={row['size']}, train_loss={row['loss']:.4f}) ---")
    print(f"  str:   {row['equation']}")
    print(f"  latex: {row['latex']}")
    # the raw equation function wants a plain array, not a DataFrame
    preds = row["equation_fn_transpose"](X_test.values)
    mse = np.mean((preds - y_test) ** 2)
    print(f"  test MSE: {mse:.4f}")
    print()
