from sklearn.model_selection import train_test_split
from Magpie import MagpieRegressor, generate_bounds
from pysr import PySRRegressor
from ffx import FFXRegressor
import pandas as pd
import numpy as np
import random
import time
pd.set_option('display.float_format', '{:.2f}'.format)
pd.set_option('display.max_colwidth', None)
import sklearn, pysr
print(sklearn.__version__)
print(pysr.__version__)
from datetime import datetime

import importlib.resources as pkg_resources
from pmlb import fetch_data, regression_dataset_names


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
    df = fetch_data(name)
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
    try:
        rsq_test = model.score(X_test, y_test)
    except:
        rsq_test = -1e6
    try:
        rsq_train = model.score(X_train, y_train)
    except:
        rsq_train = -1e6
    r = (elapsed, rsq_train, rsq_test)
    print(r)
    return r

## Run the large experiment
# 
# * FFX is deterministic so 1 run per dataset
# * PySR we just use its default setup so `n_reps` runs per dataset
# * Magpie we want to investigate hyperparameters, so `n_reps` runs per dataset per hyperparameter config.
# * For now, we just have: 
#   * 3 ways of doing interval-checking
#   * 3 values for `initevals`
#   * 3 values for `maxcohortlen`
# 
# Every run gives a row which we put in a Pandas DataFrame.

nreps = 5

from joblib import Parallel, delayed
from sklearn.base import clone

def run_single(datasetname, maxgenomelen, maxcohortlen, rep):
    print(f"DEBUG start run_single dataset={datasetname} maxgenomelen={maxgenomelen} maxcohortlen={maxcohortlen} rep={rep} {datetime.now()}", flush=True)
    dataset_n, dataset_k = srbench_datasets[datasetname][0].shape
    X_train, X_test, y_train, y_test = srbench_splits[datasetname]
    initevals = 5000
    X_bounds_key = 'test'
    X_bounds = generate_bounds(X_test)
    mr = MagpieRegressor(maxevals=25000,
                         mutprob=0.5,
                         maxgenomelen=maxgenomelen,
                         maxcohortlen=maxcohortlen,
                         initevals=initevals,
                         X_bounds=X_bounds)
    r = run_model(mr, rep, X_train, X_test, y_train, y_test)
    r = ("Magpie", datasetname, dataset_n, dataset_k, rep, X_bounds_key, initevals, maxgenomelen, maxcohortlen, *r)
    print(r, str(datetime.now()))
    return r

import os

## Checkpointing, so if a run hangs, we can restart this whole script and it will
# skip runs that are already finished

CHECKPOINT_CSV = "experiments_checkpoint.csv"
cols = ("Model", "Dataset", "N", "K", "Rep", "X_bounds", "initevals", "maxgenomelen", "maxcohortlen", "Time", "Rsq_train", "Rsq_test")

if os.path.exists(CHECKPOINT_CSV):
    checkpoint_df = pd.read_csv(CHECKPOINT_CSV, index_col=0)
    done = set(zip(
        checkpoint_df["Dataset"],
        checkpoint_df["maxgenomelen"].astype(int),
        checkpoint_df["maxcohortlen"].astype(int),
        checkpoint_df["Rep"].astype(int),
    ))
    results = [tuple(row) for row in checkpoint_df.itertuples(index=False)]
else:
    done = set()
    results = []


## run the runs in parallel

all_tasks = [
    (datasetname, maxgenomelen, maxcohortlen, rep)
    for datasetname in srbench_splits
    for maxgenomelen in (30,)
    for maxcohortlen in (3, 5, 7)
    for rep in range(nreps)
]
pending_tasks = [t for t in all_tasks if (t[0], t[1], t[2], t[3]) not in done]
print(f"{len(done)} tasks already complete, {len(pending_tasks)} remaining")

for r in Parallel(n_jobs=7, return_as='generator_unordered')(
    delayed(run_single)(datasetname, maxgenomelen, maxcohortlen, rep)
    for datasetname, maxgenomelen, maxcohortlen, rep in pending_tasks
):
    results.append(r)
    pd.DataFrame(results, columns=cols).to_csv(CHECKPOINT_CSV)

results_df = pd.DataFrame(results, columns=cols)
results_df.to_csv(f"experiments_results_{datetime.now().strftime('%Y_%m_%d_%H%M')}.csv")
