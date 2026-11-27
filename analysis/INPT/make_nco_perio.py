#!/usr/bin/env python3
# Part of the public analysis code for resident inpatient days in Korean districts, 2011-2024 (MIT licence). Inputs are read from environment variables; see README.md.
import os, sqlite3, unicodedata
import pandas as pd
N = lambda s: unicodedata.normalize("NFC", s)
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(os.path.dirname(HERE))
DB = os.environ["ST_REF_DATASET"]
con = sqlite3.connect(f"file:{DB}?mode=ro", uri=True)
for t in ["PERIO"]:
    d = pd.read_sql(f"select year, fbd_region, cases, E from disease_panel where tag='{t}'", con); d["fbd_region"] = d.fbd_region.map(N)
    assert len(d) == 3206 and d.cases.notna().all() and (d.E > 0).all()
    os.makedirs(os.path.join(ROOT, "04_분석", t), exist_ok=True)
    d[["year", "fbd_region", "cases"]].to_csv(os.path.join(ROOT, "04_분석", t, f"panel_{t}.csv"), index=False)
    d.rename(columns={"E": "E_std"})[["year", "fbd_region", "E_std"]].to_csv(os.path.join(ROOT, "02_패널", f"{t}_Estd.csv"), index=False)
    print(t, len(d))
