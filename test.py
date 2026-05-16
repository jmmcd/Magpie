import numpy as np
from sklearn.model_selection import train_test_split
from Magpie import MagpieRegressor

rng = np.random.default_rng(42)

n = 200
x0 = rng.uniform(-1, 1, n)
x1 = rng.uniform(0.01, 1, n)  # positive to keep log(0.5 * x1) defined
noise = rng.normal(0, 0.05, n)
y = x0**3 + np.log(0.5 * x1) + noise

X = np.column_stack([x0, x1])
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=0)

mr = MagpieRegressor(maxevals=500, init_prop=0.2,
                     gramfile="base_const.bnf",
                     prop_consts=0.0,
                     n_num_optimisations=0)
mr.fit(X_train, y_train)

eqs = mr.equations_
print(f"Found {len(eqs)} equations on the Pareto front\n")

for i, row in eqs.iterrows():
    print(f"--- Equation {i} (size={row['size']}, train_loss={row['loss']:.4f}) ---")
    print(f"  str:   {row['equation']}")
    print(f"  latex: {row['latex']}")
    preds = row["equation_fn_transpose"](X_test)
    mse = np.mean((preds - y_test) ** 2)
    print(f"  test MSE: {mse:.4f}")
    print()
