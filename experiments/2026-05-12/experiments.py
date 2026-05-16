import os, sys
from datetime import datetime
from pathlib import Path
import random
import time
import importlib.resources as pkg_resources
from pmlb import fetch_data, regression_dataset_names
from joblib import Parallel, delayed

import numpy as np
import pandas as pd
pd.set_option('display.float_format', '{:.2f}'.format)
pd.set_option('display.max_colwidth', None)

from sklearn.model_selection import ParameterGrid, train_test_split

# from pysr import PySRRegressor
# from ffx import FFXRegressor

from Magpie import MagpieRegressor

import sklearn, pysr
print(sklearn.__version__)
print(pysr.__version__)


## load datasets

# pmlb 1.x doesn't expose dataset_stats directly, but bundles the TSV
tsv_path = str(pkg_resources.files("pmlb") / "all_summary_stats.tsv")
stats = pd.read_csv(tsv_path, sep="\t")

# Filter: regression, 3–8 features, <2000 rows, exclude synthetic fri_c* datasets
reg_stats = stats[stats["dataset"].isin(regression_dataset_names)]
filtered = reg_stats[
    (reg_stats["n_features"] >= 3) &
    (reg_stats["n_features"] <= 8) &
    (reg_stats["n_instances"] < 2000) &
    (~reg_stats["dataset"].str.contains("fri_c"))
].sort_values("n_instances")

# Pick 10 spread evenly across the size range
step = max(1, len(filtered) // 10)
selected_names = filtered.iloc[::step].head(10)["dataset"].tolist()

print(f"Found {len(filtered)} matching datasets, selected:")
print(filtered[filtered["dataset"].isin(selected_names)][["dataset", "n_instances", "n_features"]].to_string(index=False))

# Load each as (X, y) arrays
srbench_datasets = {}
for name in selected_names:
    df = fetch_data(name, local_cache_dir=str(Path(__file__).parent.parent / "data"))
    X = df.drop("target", axis=1).values
    y = df["target"].values
    srbench_datasets[name] = (X, y)
    print(f"  {name}: X{X.shape}")

## create train-test splits

from sklearn.model_selection import train_test_split
srbench_splits = {}
for name, (X, y) in srbench_datasets.items():
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=0)
    srbench_splits[name] = (X_train, X_test, y_train, y_test)

## run_model

def run_model(model, rep, X_train, X_test, y_train, y_test):
    random.seed(rep)
    np.random.seed(rep)
    start = time.time()
    model.fit(X_train, y_train)
    end = time.time()
    elapsed = end - start
    rsq_test = model.score(X_test, y_test)
    rsq_train = model.score(X_train, y_train)
    eq_row = model.equation_.iloc[0]
    size = eq_row['size']
    n_vars_used = eq_row['n_vars_used']
    n_consts_used = eq_row['n_consts_used']
    equation = eq_row['equation']
    r = (elapsed, rsq_train, rsq_test, size, n_vars_used, n_consts_used, equation)
    return r

## Run the large experiment
# 
# * Magpie we want to investigate hyperparameters, so `n_reps` runs per dataset per hyperparameter config.
# Every run gives a row which we put in a Pandas DataFrame.

budget = 25000
nreps = 3 # 10 datasets x 3 reps is plenty.

shared = {
    'maxgenomelen': [30, 50],
    'maxcohortlen': [2, 4, 8, 15],
    'X_bounds': ['none', 'train', 'test'],
    'valsize': [0.0, 0.25, 0.5],
    'mutprob': [0.6, 1.0],
    #'prop_consts': [0.5, 0.75, 1.0],
    #'initevals': [10, 30, 80],
    #'X_bounds_margin': [0.05]
}
param_grid = [
    {**shared, 'gramfile': ['base_placeholders.bnf',
                            'extended_placeholders.bnf'],
               'n_num_optimisations': [1, 3, 5],
               'prop_consts': [0.75]},
    {**shared, 'gramfile': ['base_consts.bnf',
                            'extended_consts.bnf'],
               'n_num_optimisations': [0],
               'prop_consts': [0.0]},
]

def count_grid(g):
    result = 1
    for k in g:
        result *= len(g[k])
    return result
est_time_per_ind = 50 / 25000 # in seconds, estimated from some typical runs of 25000 individuals which took 50s
ncores = 30
print(f"we will run {count_grid(param_grid)} configurations")
print(f"with {nreps} nreps, {len(srbench_datasets)} datasets, and budget {budget}")
print(f"we have {ncores} cores")
print(f"estimate: {nreps * count_grid(param_grid) * len(srbench_datasets) * budget * est_time_per_ind / (ncores * 60 * 60)} hours")
print(f"we will have {nreps * count_grid(param_grid) * len(srbench_datasets)} rows in experiments_checkpoint.csv")
param_keys = list(param_grid.keys())
sys.exit()

import json
with open(f"param_grid_{datetime.now().strftime('%Y_%m_%d_%H%M')}.json", "w") as f:
    json.dump(param_grid, f, indent=2)

def run_single(datasetname, rep, params):
    print(f"DEBUG start run_single dataset={datasetname} params={params} rep={rep} {datetime.now()}", flush=True)
    dataset_n, dataset_k = srbench_datasets[datasetname][0].shape
    X_train, X_test, y_train, y_test = srbench_splits[datasetname]
    mr = MagpieRegressor(maxevals=budget,
                         X_bounds_test_data=X_test, # usually unused
                         **params)
    r = run_model(mr, rep, X_train, X_test, y_train, y_test)
    param_vals = tuple(params[k] for k in param_keys)
    r = (datasetname, dataset_n, dataset_k, rep, "Magpie", *param_vals, *r)
    print(r, str(datetime.now()))
    return r


## Checkpointing, so if a run hangs, we can restart this whole script and it will
# skip runs that are already finished

CHECKPOINT_CSV = "experiments_checkpoint.csv"
cols = ("Dataset", "N", "K", "Rep", "Model", *param_keys, "Time", "Rsq_train", "Rsq_test", "Size", "N_vars_used", "N_consts_used", "Equation")

def task_key(datasetname, rep, params):
    return (datasetname,) + tuple(params[k] for k in param_keys) + (rep,)

if os.path.exists(CHECKPOINT_CSV):
    checkpoint_df = pd.read_csv(CHECKPOINT_CSV, index_col=0)
    for k in param_keys:
        checkpoint_df[k] = checkpoint_df[k].astype(type(param_grid[k][0]))
    done = set(
        task_key(row["Dataset"], int(row["Rep"]), {k: row[k] for k in param_keys})
        for _, row in checkpoint_df.iterrows()
    )
    results = [tuple(row) for row in checkpoint_df.itertuples(index=False)]
else:
    done = set()
    results = []


## run the runs in parallel

all_tasks = [
    (datasetname, rep, params)
    for datasetname in srbench_splits
    for params in ParameterGrid(param_grid)
    for rep in range(nreps)
]
pending_tasks = [t for t in all_tasks if task_key(*t) not in done]
print(f"{len(done)} tasks already complete, {len(pending_tasks)} remaining")

for r in Parallel(n_jobs=ncores, return_as='generator_unordered')(
    delayed(run_single)(datasetname, rep, params)
    for datasetname, rep, params in pending_tasks
):
    results.append(r)
    pd.DataFrame(results, columns=cols).to_csv(CHECKPOINT_CSV)

results_df = pd.DataFrame(results, columns=cols)
results_df.to_csv(f"experiments_results_{datetime.now().strftime('%Y_%m_%d_%H%M')}.csv")
