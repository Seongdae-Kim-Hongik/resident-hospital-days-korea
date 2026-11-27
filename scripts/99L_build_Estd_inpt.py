#!/usr/bin/env python3
# Part of the public analysis code for resident inpatient days in Korean districts, 2011-2024 (MIT licence). Inputs are read from environment variables; see README.md.
import os, sys, json, glob, sqlite3, importlib.util, unicodedata
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE); PAN = os.path.join(ROOT, "02_패널"); SEL = os.environ.get("ST_OUT_SEL", os.path.join(ROOT, "04_분석", "INPT", "gate", "E")); os.makedirs(SEL, exist_ok=True)
RAW = os.environ.get("ST_KOSIS_RAW", "")
DB = os.environ["ST_REF_DATASET"]
def _load(nm, f):
    s = importlib.util.spec_from_file_location(nm, os.path.join(HERE, f)); m = importlib.util.module_from_spec(s); s.loader.exec_module(m); return m
FBDZIP = os.environ.get("ST_FBDZIP"); cw = _load("cw", "40_region_crosswalk.py"); m50 = _load("m50", "50_build_master_panel.py")
NFC = lambda s: unicodedata.normalize("NFC", str(s))
assert os.environ.get("ST_FBDZIP"), "ST_FBDZIP(창원 2023 교정 인구 폴더)를 지정"
YEARS = list(range(2011, 2025)); SEXES = ["남자", "여자"]
con = sqlite3.connect(f"file:{DB}?mode=ro", uri=True)
base = pd.read_sql("select year, key, fbd_region from disease_panel where tag='PERIO'", con); base["fbd_region"] = base.fbd_region.map(NFC); base["key"] = base.key.map(NFC)
PKEYS = set(base.key); assert len(base) == 3206 and len(PKEYS) == 229
K2R = dict(zip(base.key, base.fbd_region))
BANDS18 = [f"{lo}-{lo+4}" for lo in range(0, 85, 5)] + ["85+"]; OLD = ["65-69", "70-74", "75-79", "80-84", "85+"]
SIDO_HDR = {k for k in cw.SIDO_STD if len(k) >= 3 and (k.endswith("시") or k.endswith("도"))}
FIVE = [f"{lo}-{lo+4}세" for lo in range(0, 85, 5)]; OLD85 = ["85-89세", "90-94세", "95-99세", "100세 이상"]
EXTRA = {("충청북도", "청원군"): "충청북도|청주시"}; rows = []; NMISS85 = []
for y in YEARS:
    d = pd.read_csv(os.path.join(FBDZIP, f"{y}.csv"), encoding="cp949", dtype=str)
    d["행정구역별"] = d["행정구역별"].map(NFC).str.strip(); d["연령별"] = d["연령별"].map(NFC).str.strip()
    for s in SEXES: d[s] = pd.to_numeric(d[s].str.replace(",", ""), errors="coerce")
    sido = None; city = None; blocks = (d["행정구역별"] != d["행정구역별"].shift()).cumsum()
    for _, g in d.groupby(blocks, sort=False):
        nm = g["행정구역별"].iat[0]
        if nm == "전국": continue
        if nm in SIDO_HDR:
            sido = nm; city = None
            if nm != "세종특별자치시": continue
            g = g.drop_duplicates("연령별")
        sg = nm.split("-")[-1]
        if sg.endswith("시"): city = sg
        elif sg.endswith("구") and city is not None and cw.sido_std(sido) not in ("서울시", "부산시", "대구시", "인천시", "광주시", "대전시", "울산시"): continue
        else: city = None
        key = EXTRA.get((cw.sido_std(sido), sg), cw.make_key(sido, sg))
        if key not in PKEYS and key.endswith("군") and key[:-1] + "시" in PKEYS: key = key[:-1] + "시"
        for s in SEXES:
            a = g.set_index("연령별")[s]
            if pd.isna(a.get("계", np.nan)): continue
            f = a.reindex(FIVE); o = a.reindex(OLD85); assert f.notna().all(), (y, nm, s)
            NMISS85.append(int(o.isna().sum())); assert abs(a["계"] - f.sum() - o.fillna(0).sum()) < 0.5, ("85+ 결측이 실제 0명이 아님", y, nm, s, a["계"], f.sum(), o.sum())
            o = o.where(o.notna(), 0.0)
            rows.append({"year": y, "sex": s, "key": NFC(key), **{b: f[b + "세"] for b in BANDS18[:-1]}, "85+": o.sum()})
print(f"[99S72L_r2] 85+ 구간 결측 칸 {sum(NMISS85)}개 / 블록 {len(NMISS85)}개 — 모두 계 합 assert 통과(실제 0명)")
POP = pd.DataFrame(rows).groupby(["year", "key", "sex"], as_index=False).sum(numeric_only=True)
assert set(POP.key) >= PKEYS, sorted(PKEYS - set(POP.key))
cwn = POP[(POP.key == NFC("경상남도|창원시")) & (POP.year == 2023)][BANDS18].sum().sum()
print(f"[99S72] 창원 2023 인구 {cwn:,.0f} (정답 1,009,038)"); assert abs(cwn - 1009038) < 1
NATPOP = POP[POP.key.isin(PKEYS)].groupby(["year", "sex"])[BANDS18].sum()

cn = {"in": "nhis__관내입원__입내원일수", "out": "nhis__관외입원__입내원일수", "tin": "nhis__관내전체__입내원일수", "oin": "nhis__관내외래__입내원일수", "outn": "nhis__관외입원__진료실인원수"}
R = pd.read_sql("select year, key, " + ",".join(f'"{v}" as "{k}"' for k, v in cn.items()) + " from region_year where year between 2011 and 2024", con); R["key"] = R.key.map(NFC)
assert len(R) == 3206 and set(R.key) == PKEYS
miss = R["in"].isna(); z = (R.tin - R.oin)[miss]
assert (z == 0).all(), "관내입원 결측인데 관내전체−관내외래 ≠ 0"
assert R.out.notna().all() and (R.out > 0).all()
R.loc[miss, ["year", "key"]].assign(관내전체_관내외래=z.values, 처리="구조적 0(시군구 안 입원 시설 없음)").to_csv(os.path.join(SEL, "E_관내결측_구조적0_INPT.csv"), index=False, encoding="utf-8-sig")
assert (((R.tin - R.oin) - R["in"])[~miss].abs().max()) == 0
R["cases"] = R["in"].fillna(0) + R.out
R["out_dpp_ratio"] = (R.out / R.outn) / (R.out / R.outn).groupby(R.key).transform("median")
BAD = R.out_dpp_ratio > 1.5
R.loc[BAD, ["year", "key", "in", "out", "outn", "out_dpp_ratio", "cases"]].assign(rule="관외 1인당 입내원일수 > 1.5×시군구 중앙값 → 결측 → 시군구 안 선형 보간").to_csv(os.path.join(SEL, "E_결과정정칸_INPT.csv"), index=False, encoding="utf-8-sig")
R.loc[BAD, "cases"] = np.nan
R = R.sort_values(["key", "year"]); R["cases"] = R.groupby("key").cases.transform(lambda x: x.interpolate(method="linear", limit_area="inside")).round(0);
assert (R.cases % 1 == 0).all(), "정수 assert(셈 우도)"
assert R.cases.notna().all() and int(BAD.sum()) == 2, ("결과 정정 칸 수", int(BAD.sum()))
print(f"[99L INPT] 결과 원표 결함 정정 {int(BAD.sum())}칸: " + "; ".join(f"{k} {y} ratio {r:.2f}" for k, y, r in R.assign(b=BAD.values if False else False).pipe(lambda d: [(a, b, c) for a, b, c in zip(d.key, d.year, d.out_dpp_ratio) if c > 1.5])))
print(f"[99L INPT] 관내입원 결측 {int(miss.sum())}칸 → 구조적 0, 키 {R[miss].key.nunique()}개: {sorted(R[miss].key.unique())}")
NT = pd.read_csv(os.environ["ST_NHIS_NATIONAL"], encoding="utf-8-sig")
NT = NT[(NT.ITM_NM == "입내원일수") & NT.C2_NM.isin(SEXES)].copy()
def band(a):
    a = a.replace(" ", "")
    if a in ("0세", "1-4세"): return "0-4"
    if a == "85세이상": return "85+"
    if a == "합계": return None
    return a[:-1]
NT["band"] = NT.C1_NM.map(band)
tot = NT[NT.C1_NM == "합계"].set_index(["PRD_DE", "C2_NM"]).DT
NT = NT[NT.band.notna()]; NAT = NT.groupby(["PRD_DE", "C2_NM", "band"]).DT.sum().unstack()[BANDS18]; NAT.index.names = ["year", "sex"]
assert (NAT.sum(1) - tot.reindex(NAT.index).values).abs().max() == 0, "연령 합 ≠ 합계"
RATE = NAT / NATPOP[BANDS18].reindex(NAT.index)
RATE.reset_index().to_csv(os.path.join(SEL, "E_전국율_INPT.csv"), index=False, encoding="utf-8-sig")
P = base.merge(R[["year", "key", "cases"]], on=["year", "key"], validate="one_to_one")
E = np.zeros(len(P)); Esex = np.zeros(len(P)); Eage = np.zeros(len(P)); Efix = np.zeros(len(P)); off = np.zeros(len(P)); p65 = np.zeros(len(P))
agerate = (NAT.groupby(level=0).sum() / NATPOP[BANDS18].groupby(level=0).sum())
sexrate = NAT.sum(1) / NATPOP[BANDS18].sum(1)
for s in SEXES:
    X = P[["year", "key"]].merge(POP[POP.sex == s], on=["year", "key"], how="left", validate="one_to_one"); assert X[BANDS18].notna().all().all()
    pv = X[BANDS18].values
    E += (pv * RATE.xs(s, level="sex").loc[X.year.values, BANDS18].values).sum(1)
    Efix += (pv * RATE.xs(s, level="sex").loc[2015, BANDS18].values.astype(float)).sum(1)
    Eage += (pv * agerate.loc[X.year.values, BANDS18].values).sum(1)
    Esex += pv.sum(1) * sexrate.xs(s, level="sex").loc[X.year.values].values
    off += pv.sum(1); p65 += X[OLD].values.sum(1)
rs = lambda v: v * P.year.map(P.groupby("year").cases.sum()).values / pd.Series(v).groupby(P.year.values).transform("sum").values
P["offset"] = off; P["E_std"] = rs(E); P["E_sex"] = rs(Esex); P["E_std_age"] = rs(Eage); P["E_fix2015"] = Efix
P["E_crude"] = P.offset * P.year.map(P.groupby("year").cases.sum() / P.groupby("year").offset.sum()); P["pop65_share"] = p65 / off * 100
for c in ["E_std", "E_sex", "E_std_age", "E_crude"]:
    assert (P.groupby("year")[c].sum() - P.groupby("year").cases.sum()).abs().max() < 1e-6 * P.cases.sum() and (P[c] > 0).all(), c
NC = pd.DataFrame({"시군구합": P.groupby("year").cases.sum(), "전국표_건보": tot.groupby(level=0).sum(), "인구합": P.groupby("year").offset.sum()}); NC["비"] = NC.시군구합 / NC.전국표_건보
NC.to_csv(os.path.join(SEL, "E_전국합대조_INPT.csv"), encoding="utf-8-sig"); print(NC.round(4).to_string())
print(f"[99L INPT] 행 {len(P)}, 합 {int(P.cases.sum()):,}, 0칸 {(P.cases == 0).sum()}, SIR 범위 {np.round(np.percentile(P.cases / P.E_std, [0, 50, 100]), 3)}, E_crude/E_std {np.round(np.percentile(P.E_crude / P.E_std, [0, 50, 100]), 3)}")
P[["year", "fbd_region", "key", "offset", "cases", "E_std", "E_sex", "E_std_age", "E_crude", "E_fix2015", "pop65_share"]].to_csv(os.environ.get("ST_OUT_PANEL", os.path.join(PAN, "INPT_Estd3.csv")), index=False, encoding="utf-8-sig")
