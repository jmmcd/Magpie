import pandas as pd
from pathlib import Path

HERE = Path(__file__).parent
SRC = HERE / "../2026-05-11/experiments_results_2026_05_11.csv"
DST = HERE / "baselines.csv"

df = pd.read_csv(SRC)
df = df[df["Model"].isin(["FFX", "PySR"])]
df.to_csv(DST, index=False)
print(f"Written {len(df)} rows to {DST}")
