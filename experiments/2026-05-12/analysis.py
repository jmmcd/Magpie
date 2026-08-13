import json
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
_TAB10 = ['#1f77b4','#ff7f0e','#2ca02c','#d62728','#9467bd',
          '#8c564b','#e377c2','#7f7f7f','#bcbd22','#17becf']
from pathlib import Path
import statsmodels.formula.api as smf
from statsmodels.stats.anova import anova_lm
from statsmodels.formula.api import mixedlm

HERE = Path(__file__).parent
CSV = HERE / "experiments_results_2026_05_15_0258.csv"
GRID = HERE / "param_grid_2026_05_14_2204.json"

df = pd.read_csv(CSV, index_col=0)
with open(GRID) as f:
    param_grid = json.load(f)

df["Rsq_test_clamped"] = df["Rsq_test"].clip(lower=-1)
df["uops"] = df["gramfile"].map(lambda g: "extended" if "extended" in g else "base")
df = df.drop(columns=["gramfile", "prop_consts"])
df["n_opt"] = df["n_num_optimisations"]
df = df.drop(columns=["n_num_optimisations"])

# --- Diagnostic: negative Rsq_test values after filtering ---
_neg = df[df["Rsq_test"] < 0]
print(f"Negative Rsq_test: {len(_neg)} of {len(df)} rows ({100 * len(_neg) / len(df):.1f}%)")
if len(_neg) > 0:
    print(_neg.groupby("Dataset")["Rsq_test"].agg(["count", "min"]).to_string())

    # 1. Marginal negative rate by factor level
    print("\nNegative rate by factor level (marginal over all other factors):")
    for _f in ["uops", "n_opt", "maxgenomelen", "maxcohortlen", "X_bounds", "valsize", "mutprob"]:
        if _f not in df.columns:
            continue
        _rates = (
            df.groupby(_f)["Rsq_test"]
            .apply(lambda s: (s < 0).sum())
            / df.groupby(_f)["Rsq_test"].count()
        ).rename("neg_rate")
        print(f"\n  {_f}:")
        print(_rates.to_string())

    # 2. Negative rate by dataset × X_bounds
    print("\nNegative rate by dataset × X_bounds:")
    _pivot = df.assign(_is_neg=df["Rsq_test"] < 0).pivot_table(
        index="Dataset", columns="X_bounds", values="_is_neg",
        aggfunc="mean"
    )[["none", "train", "test"]]
    print(_pivot.round(3).to_string())

factors = [
    "uops" if f == "gramfile" else
    "n_opt" if f == "n_num_optimisations" else
    f for f in dict.fromkeys(k for d in param_grid for k in d.keys())
    if f not in ("prop_consts",)
]

_LABEL_MAP = {
    "X_bounds": {"none": "none", "train": "train", "test": "test"},
}

_LEVEL_ORDER = {
    "X_bounds": ["none", "train", "test"],
}

def short_label(factor, value):
    return _LABEL_MAP.get(factor, {}).get(value, str(value))

ref_config = {
    "maxgenomelen": 30,
    "maxcohortlen": 4,
    "X_bounds": "train",
    "uops": "base",
    "valsize": 0.0,
    "mutprob": 1.0,
    "n_opt": 1,
    "init_prop": 0.2,
}

# --- Diagnostic: negatives in OFAT subsets (ref config + one-factor variants) ---
_ofat_neg = []
for _factor in factors:
    _levels = _LEVEL_ORDER.get(_factor, sorted(df[_factor].unique()))
    for _level in _levels:
        _mask = pd.Series(True, index=df.index)
        for _f, _ref_val in ref_config.items():
            if _f != _factor:
                _mask &= (df[_f] == _ref_val)
        _mask &= (df[_factor] == _level)
        _n_neg = (df[_mask]["Rsq_test"] < 0).sum()
        _n_total = _mask.sum()
        if _n_neg > 0:
            _ofat_neg.append(f"  {_factor}={_level}: {_n_neg}/{_n_total} negative")
if _ofat_neg:
    print("\nNegatives found in OFAT plot subsets:")
    for _line in _ofat_neg:
        print(_line)
else:
    print("\nOFAT plot subsets: no negative Rsq_test values.")

_ds_sorted = sorted(df["Dataset"].unique())
datasets = [ds for ds in _ds_sorted if "1096" not in str(ds)] + [ds for ds in _ds_sorted if "1096" in str(ds)]
dataset_colors = {ds: _TAB10[i] for i, ds in enumerate(datasets)}

n = len(factors)
ncols = 4
nrows = (n + ncols - 1) // ncols
YMIN, YMAX = -0.5, 1.0

fig, axes = plt.subplots(nrows, ncols, figsize=(7, 3.5 * nrows), sharey=True)
axes = axes.flatten()
rng = np.random.default_rng(42)

for i, (ax, factor) in enumerate(zip(axes, factors)):
    levels = _LEVEL_ORDER.get(factor, sorted(df[factor].unique()))
    labels = [short_label(factor, v) + ("*" if ref_config.get(factor) == v else "") for v in levels]

    for xi, level in enumerate(levels):
        # OFAT: fix all other factors at reference, vary only this one
        mask = pd.Series(True, index=df.index)
        for f, ref_val in ref_config.items():
            if f != factor:
                mask &= (df[f] == ref_val)
        mask &= (df[factor] == level)
        sub = df[mask]
        median_y = np.clip(sub["Rsq_test"].median(), YMIN, YMAX)
        ax.hlines(median_y, xi - 0.35, xi + 0.35, colors="black", linewidth=1.5, zorder=1)
        for _, row in sub.iterrows():
            jitter = rng.uniform(-0.2, 0.2)
            ax.scatter(xi + jitter, np.clip(row["Rsq_test"], YMIN, YMAX),
                       color=dataset_colors[row["Dataset"]], alpha=0.7, s=18, linewidths=0, zorder=2)

    ax.set_xticks(range(len(levels)))
    ax.set_xticklabels(labels)
    ax.set_title(factor)
    ax.axhline(0, color="gray", linewidth=0.6, linestyle="--")
    if i % ncols == 0:
        ax.set_ylabel("R² (test)")
        ax.yaxis.set_major_formatter(ticker.FormatStrFormatter("%.1f"))
    else:
        ax.tick_params(labelleft=False)

axes[0].set_ylim(-0.6, YMAX)

for ax in axes[n:]:
    ax.set_visible(False)

fig.tight_layout()
fig.subplots_adjust(wspace=0.05)  # tight_layout resets wspace, so re-apply
plt.savefig(HERE / "main_effects.png", dpi=150)
plt.savefig(HERE / "main_effects.pdf")
# plt.show()

# --- Linear model (ANOVA) ---
# Treat all factors as categorical; Dataset is a blocking variable.
# Numeric factors wrapped in C() to force treatment contrast coding.
factor_terms = " + ".join(f"C({f})" for f in factors)
formula = f"Rsq_test_clamped ~ C(Dataset) + {factor_terms}"
model = smf.ols(formula, data=df).fit()
aov = anova_lm(model, typ=2)
print(model.summary())
print("\nType II ANOVA table:")
print(aov)

# --- Mixed effects model: Dataset as random intercept ---
# Equivalent to a paired analysis: dataset variance is partitioned out as a
# random effect rather than consumed as fixed coefficients.
mixed_formula = f"Rsq_test_clamped ~ {factor_terms}"
mixed_model = mixedlm(mixed_formula, data=df, groups=df["Dataset"]).fit()
print("\nMixed model summary:")
print(mixed_model.summary())

# --- Interaction model ---
df_sub = df.copy()

# Mixed model with all pairwise interactions among all factors.
terms = " + ".join(f"C({f})" for f in factors)
interaction_formula = f"Rsq_test_clamped ~ ({terms}) ** 2"
mixed_model2 = mixedlm(interaction_formula, data=df_sub, groups=df_sub["Dataset"]).fit()
print("\nMixed model (all pairwise interactions):")
print(mixed_model2.summary())

# --- LaTeX table: all main effects + significant interactions ---
def _clean_term(name):
    name = re.sub(r"C\((\w+)\)\[T\.([^\]]+)\]", lambda m: f"\\mi|{m.group(1)}={m.group(2)}|", name)
    return name.replace(":", " $\\times$ ")

def _fmt_coef(x):
    return f"{x:.2f}"

def _fmt_pval(x):
    return f"{x:.3f}"[1:]  # strip leading zero: 0.003 -> .003

def _post_process_latex(tex):
    return tex.replace("\\begin{table}\n", "\\begin{table}\n\\centering\n")

# _ci = mixed_model2.conf_int()
_res = pd.DataFrame({
    "Coef.": mixed_model2.fe_params,
    # "Std. Err.": mixed_model2.bse,
    # "z":         mixed_model2.tvalues,
    "p":     mixed_model2.pvalues,
    # "CI 2.5%":  _ci[0],
    # "CI 97.5%": _ci[1],
})
def _row_sort_key(name):
    # Operates on raw statsmodels names e.g. C(maxcohortlen)[T.15]
    m = re.match(r"C\((\w+)\)\[T\.([^\]]+)\]", name)
    if m:
        factor, level = m.group(1), m.group(2)
        try:
            return (factor, float(level), "")
        except ValueError:
            return (factor, 0.0, level)
    return ("", 0.0, name)

_main = _res[~_res.index.str.contains(":")].copy()
_main = _main.loc[sorted(_main.index, key=_row_sort_key)]  # sort before cleaning
_main.index = pd.Index([_clean_term(t) for t in _main.index])

_sig_int = _res[_res.index.str.contains(":") & (_res["p"] < 0.05)].copy()
_sig_int.index = pd.Index([_clean_term(t) for t in _sig_int.index])

_table = pd.concat([_main, _sig_int])
_latex = (
    _table.style
    .format({"Coef.": _fmt_coef, "p": _fmt_pval})
    .to_latex(
        caption="Main effects and significant pairwise interactions from mixed-effects model (p < 0.05).",
        label="tab:sig_interactions",
        hrules=True,
    )
)
(HERE / "sig_interactions.tex").write_text(_post_process_latex(_latex))
print(f"\nLaTeX table written to {HERE / 'sig_interactions.tex'}")

# --- Best configuration by average within-dataset rank ---
# Rank configs within each (Dataset, Rep) group; rank 1 = highest Rsq_test_clamped.
# Averaging across groups rewards consistency across datasets, not raw magnitude.
df["_rank"] = df.groupby(["Dataset", "Rep"])["Rsq_test_clamped"].rank(
    ascending=False, method="average"
)
avg_ranks = df.groupby(factors)["_rank"].mean().reset_index()
avg_ranks = avg_ranks.rename(columns={"_rank": "mean_rank"}).sort_values("mean_rank")
print("\nTop 10 configurations by average rank:")
print(avg_ranks.head(10).to_string(index=False))
print("\nBest configuration:")
print(avg_ranks.iloc[0].to_string())

_top10 = avg_ranks.head(10).copy()
_top10.insert(0, "rank", range(1, 11))
_float_factors = [f for f in factors if avg_ranks[f].dtype == float]
_fmt = {f"\\mi|{f}|": "{:.1f}" for f in _float_factors} | {"mean rank": "{:.2f}"}
_top10 = _top10.rename(
    columns={f: f"\\mi|{f}|" for f in factors} | {"mean_rank": "mean rank"}
)
_latex_top10 = (
    _top10.style
    .format(_fmt)
    .hide(axis="index")
    .to_latex(
        caption="Top 10 configurations by average within-dataset rank (rank 1 = best).",
        label="tab:top10",
        hrules=True,
    )
)
_latex_top10_processed = (
    _post_process_latex(_latex_top10)
    .replace("\\begin{table}", "\\begin{table*}")
    .replace("\\end{table}", "\\end{table*}")
)
(HERE / "top10_configs.tex").write_text(_latex_top10_processed)
print(f"\nTop-10 LaTeX table written to {HERE / 'top10_configs.tex'}")

# --- OFAT table: median Time and Size for each hyperparameter level ---
_ofat_rows = []
for factor in factors:
    _levels = _LEVEL_ORDER.get(factor, sorted(df[factor].unique()))
    for level in _levels:
        mask = pd.Series(True, index=df.index)
        for f, ref_val in ref_config.items():
            if f != factor:
                mask &= (df[f] == ref_val)
        mask &= (df[factor] == level)
        sub = df[mask]
        _is_ref = ref_config.get(factor) == level
        _ofat_rows.append({
            "Factor": f"\\mi|{factor}|",
            "Level": f"{level}{'*' if _is_ref else ''}",
            "Time (s)": sub["Time"].median(),
            "Size": sub["Size"].median(),
        })

_ofat_df = pd.DataFrame(_ofat_rows).set_index(["Factor", "Level"])
print("\nOFAT table: Time and Size (* = reference level):")
print(_ofat_df.to_string())

_latex_ofat = (
    _ofat_df.style
    .format({"Time (s)": "{:.0f}", "Size": "{:.0f}"})
    .to_latex(
        caption="OFAT analysis: median elapsed time and equation size per hyperparameter level. * denotes the reference configuration.",
        label="tab:ofat_time_size",
        hrules=True,
        multirow_align="t",
    )
)
(HERE / "ofat_time_size.tex").write_text(_post_process_latex(_latex_ofat))
print(f"\nOFAT LaTeX table written to {HERE / 'ofat_time_size.tex'}")

# --- Comparison: SBGP (ref config) vs PySR ---
_baselines = pd.read_csv(HERE / "baselines.csv")
print("\nMean time (s) for baseline models:")
print(_baselines.groupby("Model")["Time"].mean().round(1).to_string())

_pysr = _baselines[_baselines["Model"] == "PySR"][["Dataset", "Rep", "Rsq_test"]].copy()
_pysr["Rsq_test_clamped"] = _pysr["Rsq_test"].clip(lower=-1)
_pysr["Model"] = "PySR"

_ref_mask = pd.Series(True, index=df.index)
for _f, _ref_val in ref_config.items():
    _ref_mask &= (df[_f] == _ref_val)
_sbgp = df[_ref_mask][["Dataset", "Rep", "Rsq_test", "Rsq_test_clamped"]].copy()
_sbgp["Model"] = "SBGP"

_compare = pd.concat([_sbgp, _pysr], ignore_index=True)
print("\nSBGP (ref config) vs PySR — observations per model:")
print(_compare.groupby("Model")["Rsq_test"].agg(["count", "mean", "median"]).to_string())

_cmp_model = mixedlm("Rsq_test_clamped ~ C(Model)", data=_compare, groups=_compare["Dataset"]).fit()
print("\nMixed model: SBGP (ref config) vs PySR (Dataset as random intercept):")
print(_cmp_model.summary())

# --- Comparison plot: FFX vs PySR vs SBGP (ref config) ---
_ffx = _baselines[_baselines["Model"] == "FFX"][["Dataset", "Rep", "Rsq_test"]].copy()
_ffx["Model"] = "FFX"

_cmp_plot_data = {"FFX": _ffx, "PySR": _pysr, "SBGP": _sbgp}
_cmp_labels = ["FFX", "PySR", "SBGP"]
_rng_cmp = np.random.default_rng(42)

fig_cmp, ax_cmp = plt.subplots(figsize=(3.8, 3.5))
for xi, _mname in enumerate(_cmp_labels):
    _sub = _cmp_plot_data[_mname]
    _median_y = np.clip(_sub["Rsq_test"].median(), YMIN, YMAX)
    ax_cmp.hlines(_median_y, xi - 0.35, xi + 0.35, colors="black", linewidth=1.5, zorder=1)
    for _, _row in _sub.iterrows():
        _jitter = _rng_cmp.uniform(-0.2, 0.2)
        ax_cmp.scatter(xi + _jitter, np.clip(_row["Rsq_test"], YMIN, YMAX),
                       color=dataset_colors[_row["Dataset"]], alpha=0.7, s=18, linewidths=0, zorder=2)

ax_cmp.set_xticks(range(3))
ax_cmp.set_xticklabels(_cmp_labels)
ax_cmp.set_ylabel("R² (test)")
ax_cmp.set_ylim(-0.6, YMAX)
ax_cmp.axhline(0, color="gray", linewidth=0.6, linestyle="--")
ax_cmp.yaxis.set_major_formatter(ticker.FormatStrFormatter("%.1f"))
_cmp_legend_handles = [
    ax_cmp.scatter([], [], color=dataset_colors[ds], s=18, label=ds.split("_")[0])
    for ds in datasets
]
ax_cmp.legend(handles=_cmp_legend_handles, title="Dataset",
              bbox_to_anchor=(1.02, 1), loc="upper left", frameon=False, fontsize=7)
fig_cmp.tight_layout()
fig_cmp.savefig(HERE / "comparison.png", dpi=150)
fig_cmp.savefig(HERE / "comparison.pdf")
print(f"\nComparison plot written to {HERE / 'comparison.png'}")
