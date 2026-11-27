#!/usr/bin/env python3
# Part of the public analysis code for resident inpatient days in Korean districts, 2011-2024 (MIT licence). Inputs are read from environment variables; see README.md.
import os, json, hashlib
import pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(os.path.dirname(HERE)); PAN = os.path.join(ROOT, "02_패널")
POOL = "EXT_canonical_pool_INPT.csv"
P = pd.read_csv(os.path.join(PAN, POOL), encoding="utf-8-sig"); assert len(P) == 229 * 14 and P.drop(columns=["vil_n_sqlite", "vil_n_v2"], errors="ignore").notna().all().all()
FORCED = ["urban_population_share", "physicians_per_1000", "sex_ratio", "fiscal_independence", "pop80_share"]
NEG = ["hit_and_run_rate", "traffic_culture_index", "fires_per_10000"]
SENS = ["beds_per_1000", "ltc_hospitals_per100k_65", "hospitals_per100k", "tertiary_educated_share", "living_alone_elderly", "mort_allcause_asr"]
HYP = ["walking_rate", "moderate_pa_rate", "public_sport_fac_per100k", "village_sport_fac_per100k",
       "obesity_rate", "smoking_rate", "high_risk_drinking", "basic_livelihood_rate", "screening_uptake"]
X = lambda L: ["ext__" + v for v in L]; XK = lambda v: "ext__" + v
EXPECT = {
 "walking_rate": "- (걷기 실천율: 습관적 신체활동 → 심혈관·대사·근골격 질환과 낙상 감소 → 입원 일수 감소; 체육 경로 주 가설)",
 "moderate_pa_rate": "- (중강도 이상 신체활동 실천율: 만성질환 악화·입원 감소; 체육 경로 주 가설)",
 "public_sport_fac_per100k": "- (공공체육시설 10만 명당: 지역 신체활동 기회 구조; 체육 경로 주 가설)",
 "village_sport_fac_per100k": "- (마을체육시설 10만 명당: 생활권 신체활동 기회, 정정판)",
 "obesity_rate": "+ (KCHS 비만율: 당뇨·심혈관·관절 질환 → 입원; 2019는 기준 데이터셋 선형 보간)",
 "smoking_rate": "+ (KCHS 현재흡연율: 호흡기·심혈관·암 입원)",
 "high_risk_drinking": "+ (KCHS 고위험음주율: 간질환·손상·알코올 관련 입원)",
 "basic_livelihood_rate": "+ (기초생활 수급률, 공식 복원판: 박탈 → 다질환·장기 입원(사회적 입원 포함); 대안 − = 접근 장벽 — 결과 전 양측 해석으로 고정)",
 "screening_uptake": "- (일반건강검진 수검률: 예방 도달 → 입원 감소; 대안 + = 의료 이용 성향 — 결과 전 양측 해석으로 고정)",
}
SETS = {"S1_narrative_full": X(HYP),
        "S2_core": X(["walking_rate", "public_sport_fac_per100k", "obesity_rate", "basic_livelihood_rate", "screening_uptake"]),
        "S3_physical_activity": X(["walking_rate", "moderate_pa_rate", "public_sport_fac_per100k", "village_sport_fac_per100k"]),
        "S4_lifestyle_risk": X(["obesity_rate", "smoking_rate", "high_risk_drinking"]),
        "S5_deprivation_care": X(["basic_livelihood_rate", "screening_uptake"]),
        "S6_risk_deprivation": X(["obesity_rate", "smoking_rate", "high_risk_drinking", "basic_livelihood_rate", "screening_uptake"])}
pov = [{"year": 2023, "fbd_region": "경상남도창원시", "population": 1009038, "source": "KOSIS DT_1YL20651E 2023 통합창원시 계(교정 인구 ~/sk_popfix)"}]
allv = ["inpt_cases", "inpt_pop"] + FORCED + NEG + SENS + HYP
LAB_M = "sex-and-age-standardised expected inpatient days, 2 sexes x 18 five-year age bands (national National Health Insurance inpatient-day rates by sex and age, KOSIS TX_35001_A047, same year; district resident population by sex and age), rescaled so that annual expected equals observed totals, 2011-2024"
WIN = {"W1_2011-15": [2011, 2015], "W2_2016-19": [2016, 2019], "W3_2020-24": [2020, 2024]}
def make():
    return {
     "tag": "INPT", "outcome_column": "ext__inpt_cases",
     "outcome_label_en": "Inpatient days of district residents insured by National Health Insurance, in-district plus out-of-district hospital admissions (NHIS district medical-use statistics, annual)",
     "disease_en": "Resident inpatient care use (inpatient days)",
     "disease_ko": "거주지 기준 입원 이용(건보 관내+관외 입원 입내원일수)",
     "pathway": "체육 경로 std3(비감염 결과, 13 §3 항목 0): 체육 세트(걷기·중강도·공공/마을체육시설) + 생활습관 위험(비만·흡연·고위험음주) + 박탈·예방 도달(수급률·검진 수검률). 신체활동 4지표 전부 null이면 v2(입원 이용 역학 세트: 병상 공급·요양병원·박탈·사망 부담).",
     "selection_basis": "데스크 리젝 재투고(33 §1): N0 선별 건너뜀, 종속변수 정의만 36 자료지도에서 가져옴",
     "offset": "ext__inpt_pop", "exclude_years": [],
     "forced": X(FORCED), "negative_controls": X(NEG),
     "offset_E_column": "E_std", "offset_E_file": "INPT_Estd3.csv", "offset_E_label": LAB_M,
     "graph": "593-edge contact graph, ST_GRAPH=04_분석/INPT/graph_593.graph",
     "population_override": pov, "sensitivity_only_extra": X(SENS),
     "version": "std3 physical-activity pathway (resident inpatient days: physical activity and sport facilities, lifestyle risk, deprivation and prevention reach), 2011-2024",
     "fixed_on": "2026-11-27 (공변량 모형 적합 전 고정; 열람한 것 = E 빌더 출력(SIR 범위·전국합 비)과 풀 기술통계뿐, 계수 미열람)",
     "sets": SETS, "transform": {XK(v): "none" for v in FORCED + NEG + SENS + HYP},
     "expected_direction": {XK(k): v for k, v in EXPECT.items()},
     "evidence": {XK(k): "근거 문헌은 S04 정독 단계에서 확정(적합 전에는 방향과 기전만 고정)" for k in HYP},
     "rules": {"likelihood": "forced-only 모형 WAIC 낮은 쪽(Poisson/NB)", "set_choice": "WAIC 최소, 4 이내면 변수 적은 세트", "credible": "95% CrI 1 제외", "vif_drop": 10, "vif_action": "drop", "multiverse_specs": 200,
               "main_model_note": "주모형 = S1_narrative_full(가설 9 + 강제통제 5). 세트를 다시 바꾸지 않음. VIF>10 후보는 80이 반복 제거, 강제통제 VIF>10이 남으면 80 뒤 게이트에서 중단.",
               "years": [2011, 2024], "windows": WIN, "covid_post_from": 2020, "windows_note": "시리즈 표준 3창(결과 전 고정)."},
     "note": "INPT = 거주지 시군구 건보 가입자의 입원 입내원일수(관내 의료기관 + 관외 의료기관). 관내 입원 결측 44칸(5개 군·시, 입원 시설 없는 해) = 관내전체−관내외래=0 확인 뒤 0. 시군구 합/전국 건보 입원 표 비 1.176–1.215(전국 표 범위 차, 율의 연령·성 구조만 사용).",
     "extra_columns": {"ext__" + v: [POOL, v] for v in allv}}
OUT = {}
VERS = {
 "std3": None,
 "s577": ("민감도: 그래프만 queen 577", {}),
 "s2": ("민감도: 80세 이상 비율 강제통제 없음", {"forced": X(FORCED[:4])}),
 "sc": ("민감도: 비표준화(조) E(인구 × 전국 연도 율)", {"offset_E_column": "E_crude", "offset_E_label": "crude E (population x national annual rate)"}),
 "s8": ("민감도: 연령만 표준화 E(성 미구분, 18구간)", {"offset_E_column": "E_std_age", "offset_E_label": "age-only standardised E"}),
 "sexE": ("민감도: 대안 E(성별만 율)", {"offset_E_column": "E_sex", "offset_E_label": "sex-only standardised E"}),
 "fix15": ("민감도: 2015 성×연령 율 고정 E(연도 재척도 없음)", {"offset_E_column": "E_fix2015", "offset_E_label": "sex-and-age standardised E with 2015 national rates fixed (no annual rescaling)"}),
 "x2021": ("민감도: 2020–2021 제외(코로나19 시기 입원 감소)", {"exclude_years": [2020, 2021]}),
 "supply": ("민감도: 공급 공변량 강제 추가(의료기관 병상 1천 명당·요양병원 65+ 10만 명당)", {"forced": X(FORCED + ["beds_per_1000", "ltc_hospitals_per100k_65"])}),
}
for ver, spec in VERS.items():
    j = make()
    if spec:
        j["sens_note"] = spec[0]; j.update(spec[1])
        if "forced" in spec[1]: j["sensitivity_only_extra"] = [v for v in j["sensitivity_only_extra"] if v not in j["forced"]]
    OUT[ver] = j
H = []
for ver, j in OUT.items():
    fp = os.path.join(HERE, f"sets_INPT_{ver}.json"); json.dump(j, open(fp, "w"), ensure_ascii=False, indent=1)
    H.append(f"sets_INPT_{ver}.json sha256={hashlib.sha256(open(fp, 'rb').read()).hexdigest()}")
rows = [dict(variable=XK(v), role="hypothesis", expected_direction=EXPECT[v], in_S1=True) for v in HYP] + [dict(variable=XK(v), role="forced_control", expected_direction="통제(방향 가설 없음)", in_S1=True) for v in FORCED] + \
       [dict(variable=XK(v), role="negative_control", expected_direction="0(연관 없음 기대)", in_S1=False) for v in NEG] + [dict(variable=XK(v), role="sensitivity_only", expected_direction="민감도 전용", in_S1=False) for v in SENS]
hp = os.path.join(HERE, "hypothesis_table_INPT.csv"); pd.DataFrame(rows).to_csv(hp, index=False, encoding="utf-8-sig")
H.append("hypothesis_table_INPT.csv sha256=" + hashlib.sha256(open(hp, "rb").read()).hexdigest())
open(os.path.join(HERE, "sets_INPT_hash_before_fit.txt"), "w").write("# sets_INPT json and hypothesis table SHA-256, before fitting\n" + "\n".join(H) + "\n")
print("\n".join(H))
