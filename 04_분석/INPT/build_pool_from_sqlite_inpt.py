#!/usr/bin/env python3
# Part of the public analysis code for resident inpatient days in Korean districts, 2011-2024 (MIT licence). Inputs are read from environment variables; see README.md.
import os, hashlib, sqlite3, unicodedata
import numpy as np, pandas as pd
NFC = lambda s: unicodedata.normalize("NFC", str(s))
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(os.path.dirname(HERE)); PAN = os.path.join(ROOT, "02_패널"); GT = os.path.join(HERE, "gate")
DB = os.environ["ST_REF_DATASET"]
sha = lambda p: hashlib.sha256(open(p, "rb").read()).hexdigest()
h0 = sha(DB); m0 = os.stat(DB).st_mtime
con = sqlite3.connect(f"file:{DB}?mode=ro", uri=True)
M = {"walking_rate": "walking_practice_rate", "moderate_pa_rate": "moderate_pa_rate", "public_sport_fac_per100k": "public_sports_facilities_per_100k", "village_sport_fac_per100k": "village_sports_facilities_per_100k",
     "obesity_rate": "obesity_rate", "high_risk_drinking": "high_risk_drinking_rate", "smoking_rate": "current_smoking_rate", "screening_uptake": "screening_uptake",
     "basic_livelihood_rate": "basic_livelihood_rate", "urban_population_share": "urban_population_share", "physicians_per_1000": "physicians_per_1000", "sex_ratio": "sex_ratio", "fiscal_independence": "fiscal_independence",
     "pop80_share": "pop80_share", "hit_and_run_rate": "hit_and_run_rate", "traffic_culture_index": "traffic_culture_index", "fires_per_10000": "fires_per_10000", "tertiary_educated_share": "tertiary_educated_share",
     "living_alone_elderly": "older_single_household_share",
     "care_institutions_per100k": "care_institutions_per_100k", "ltc_hospitals_per100k_65": "ltc_hospitals_per_100k_65"}
cols = ",".join(f'"{v}_raw"' for v in M.values())
D = pd.read_sql(f'select tag, year, fbd_region, sigungu_cd, population, {cols} from disease_panel where year between 2011 and 2024', con); D["fbd_region"] = D.fbd_region.map(NFC)
D["sigungu_cd"] = D.sigungu_cd.astype(str).str.replace(r"\.0$", "", regex=True); key = D[D.tag == "ALZ"][["year", "fbd_region", "sigungu_cd", "population"]].sort_values(["year", "fbd_region"]).reset_index(drop=True)
assert len(key) == 229 * 14 and key.fbd_region.nunique() == 229
P = key.copy(); log = []
for k, v in M.items():
    c = v + "_raw"; ok = D.groupby("tag")[c].apply(lambda s: s.notna().sum() == 3206); tags = [t for t in ok[ok].index if not (k == "public_sport_fac_per100k" and t == "BRCA")]
    assert tags, (k, "229×14 모두 있는 태그 없음")
    W = D[D.tag.isin(tags)].pivot_table(index=["year", "fbd_region"], columns="tag", values=c)
    spread = float((W.max(1) - W.min(1)).abs().max()); assert spread < 1e-9, (k, spread)
    s = W.iloc[:, 0].rename(k); P = P.merge(s.reset_index(), on=["year", "fbd_region"], how="left", validate="one_to_one"); log.append((k, c, len(tags), tags[0], spread))
V = pd.read_csv(os.path.join(GT, "tables", "T07b_village_facilities_corrected.csv"), dtype={"sigungu_cd": str}, encoding="utf-8-sig")
P = P.merge(V, on=["sigungu_cd", "year"], how="left", validate="one_to_one"); assert P.loc[P.year != 2015, "vil_n_sqlite"].notna().all() and P.loc[P.year == 2015, "vil_n_sqlite"].isna().all()
chk = (P.vil_n_sqlite / P.population * 1e5 - P.village_sport_fac_per100k).abs().max(); assert chk < 1e-6, chk
v15 = P[P.year == 2015].village_sport_fac_per100k.copy()
P["village_sport_fac_per100k_sqlite"] = P.village_sport_fac_per100k
P["village_sport_fac_per100k"] = (P.vil_n_v2 / P.population * 1e5)
P["village_sport_fac_per100k"] = P.groupby("fbd_region").village_sport_fac_per100k.transform(lambda s: s.interpolate(method="linear", limit_area="inside"))
assert P.village_sport_fac_per100k.notna().all(), P[P.village_sport_fac_per100k.isna()][["fbd_region", "year"]]
ch = P[(P.village_sport_fac_per100k - P.village_sport_fac_per100k_sqlite).abs() > 1e-9]
ch[["year", "fbd_region", "sigungu_cd", "population", "vil_n_sqlite", "vil_n_v2", "village_sport_fac_per100k_sqlite", "village_sport_fac_per100k"]].assign(var="village_sport_fac_per100k", rule="마을체육시설 정정(06 사본) + 결측 칸 시군구 안 선형 보간").to_csv(
    os.path.join(GT, "tables", "G00_village_changed_cells.csv"), index=False, encoding="utf-8-sig")
print("마을체육시설 정정 칸(패널 값 변경):", len(ch), "군위 2023 =", round(float(P.loc[P.fbd_region.str.contains("군위") & (P.year == 2023), "village_sport_fac_per100k"].iloc[0]), 1))
E = pd.read_csv(os.path.join(PAN, "INPT_Estd3.csv"), encoding="utf-8-sig"); E["fbd_region"] = E.fbd_region.map(NFC)
P = P.merge(E[["year", "fbd_region", "cases", "offset", "pop65_share"]].rename(columns={"cases": "inpt_cases", "offset": "inpt_pop"}), on=["year", "fbd_region"], how="left", validate="one_to_one")
assert P.inpt_cases.notna().all() and (P.inpt_pop > 0).all()
dp = (P.inpt_pop / P.population - 1).abs(); print("주민등록 성×연령 합 / SQLite population 최대 상대차", round(float(dp.max()), 4))
miss = P.drop(columns=["vil_n_sqlite", "vil_n_v2"]).isna().sum(); assert miss.sum() == 0, miss[miss > 0]
P.drop(columns=["vil_n_sqlite", "vil_n_v2"]).to_csv(os.path.join(PAN, "EXT_canonical_pool_INPT_base.csv"), index=False, encoding="utf-8-sig")
con.close(); h1 = sha(DB); m1 = os.stat(DB).st_mtime
open(os.path.join(GT, "G00_base_hash.txt"), "w").write(f"기준 SQLite {DB}\n전 sha256={h0} mtime={m0}\n후 sha256={h1} mtime={m1}\n같음={h0 == h1 and m0 == m1}\n" + "\n".join(f"{a}<-{b} 태그 {c}개(대표 {d}) 태그 간 최대차 {e}" for a, b, c, d, e in log) + "\n")
print("풀", P.shape, "기준 파일 해시 전후 같음:", h0 == h1 and m0 == m1)
