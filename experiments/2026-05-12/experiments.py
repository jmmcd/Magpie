import os, sys
from datetime import datetime
from pathlib import Path
import time
import importlib.resources as pkg_resources
from pmlb import fetch_data, regression_dataset_names
from joblib import Parallel, delayed

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

def run_model(model, X_train, X_test, y_train, y_test):
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
    'init_prop': [0.12, 0.20],
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

import argparse
import json

parser = argparse.ArgumentParser()
parser.add_argument('--no-dry-run', action='store_true', dest='run',
                    help='Actually run experiments (default: dry-run only)')
args = parser.parse_args()
dry_run = not args.run

# param_keys: all keys appearing in any sub-grid (order-preserving, deduped)
param_keys = list(dict.fromkeys(k for subgrid in param_grid for k in subgrid.keys()))

all_configs = list(ParameterGrid(param_grid))
n_configs = len(all_configs)
est_time_per_ind = 50 / 25000
ncores = 30
print(f"Param grid: {n_configs} configurations")
print(f"  nreps={nreps}, datasets={len(srbench_datasets)}, budget={budget}, cores={ncores}")
print(f"  Total rows if complete: {nreps * n_configs * len(srbench_datasets)}")
print(f"  Estimated time: {nreps * n_configs * len(srbench_datasets) * budget * est_time_per_ind / (ncores * 3600):.1f} hours")


## Checkpointing

CHECKPOINT_CSV = "experiments_checkpoint.csv"
cols = ("Dataset", "N", "K", "Rep", "Model", *param_keys, "Time", "Rsq_train", "Rsq_test", "Size", "N_vars_used", "N_consts_used", "Equation")

def task_key(datasetname, rep, params):
    return (datasetname,) + tuple(params[k] for k in param_keys) + (rep,)

first_config = all_configs[0]

if os.path.exists(CHECKPOINT_CSV):
    checkpoint_df = pd.read_csv(CHECKPOINT_CSV, index_col=0)
    for k in param_keys:
        checkpoint_df[k] = checkpoint_df[k].astype(type(first_config[k]))
    done = set(
        task_key(row["Dataset"], int(row["Rep"]), {k: row[k] for k in param_keys})
        for _, row in checkpoint_df.iterrows()
    )
    results = [tuple(row) for row in checkpoint_df.itertuples(index=False)]
else:
    checkpoint_df = None
    done = set()
    results = []

all_tasks = [
    (datasetname, rep, params)
    for datasetname in srbench_splits
    for params in all_configs
    for rep in range(nreps)
]
all_task_keys = {task_key(*t) for t in all_tasks}
pending_tasks = [t for t in all_tasks if task_key(*t) not in done]

## Dry-run report

print(f"\nCheckpoint: {len(done)} tasks done, {len(pending_tasks)} pending")

if checkpoint_df is not None:
    stale_mask = [
        task_key(row["Dataset"], int(row["Rep"]), {k: row[k] for k in param_keys})
        not in all_task_keys
        for _, row in checkpoint_df.iterrows()
    ]
    stale_df = checkpoint_df[stale_mask]
    if len(stale_df):
        print(f"\nWARNING: {len(stale_df)} checkpoint rows do not match the current param grid:")
        print(stale_df[["Dataset", "Rep", *param_keys]].to_string(index=False))
    else:
        print("No stale checkpoint rows.")

print(f"\nPending tasks ({len(pending_tasks)}):")
i = 0
for datasetname, rep, params in pending_tasks:
    print(f"  dataset={datasetname} rep={rep} " + " ".join(f"{k}={v}" for k, v in params.items()))
    i += 1
    if i > 10:
        print("...")
        break

if dry_run:
    print("\nDry run complete. To run experiments: python experiments.py --no-dry-run")
    sys.exit(0)

## Live run

with open(f"param_grid_{datetime.now().strftime('%Y_%m_%d_%H%M')}.json", "w") as f:
    json.dump(param_grid, f, indent=2)

def run_single(datasetname, rep, params):
    print(f"DEBUG start run_single dataset={datasetname} params={params} rep={rep} {datetime.now()}", flush=True)
    try:
        dataset_n, dataset_k = srbench_datasets[datasetname][0].shape
        X_train, X_test, y_train, y_test = srbench_splits[datasetname]
        mr = MagpieRegressor(maxevals=budget,
                             X_bounds_test_data=X_test,
                             random_state=rep,
                             **params)
        r = run_model(mr, X_train, X_test, y_train, y_test)
        param_vals = tuple(params[k] for k in param_keys)
        r = (datasetname, dataset_n, dataset_k, rep, "Magpie", *param_vals, *r)
        print(r, str(datetime.now()))
        return r
    except Exception as e:
        print(f"ERROR dataset={datasetname} rep={rep} params={params}: {e}", flush=True)
        return None

for r in Parallel(n_jobs=ncores, return_as='generator_unordered')(
    delayed(run_single)(datasetname, rep, params)
    for datasetname, rep, params in pending_tasks
):
    if r is None:
        continue
    results.append(r)
    pd.DataFrame(results, columns=cols).to_csv(CHECKPOINT_CSV)

results_df = pd.DataFrame(results, columns=cols)
results_df.to_csv(f"experiments_results_{datetime.now().strftime('%Y_%m_%d_%H%M')}.csv")
