#!/usr/bin/env python3
# Part of the public analysis code for resident inpatient days in Korean districts, 2011-2024 (MIT licence). Inputs are read from environment variables; see README.md.
import os, json, hashlib
H = os.path.dirname(os.path.abspath(__file__))
j = json.load(open(os.path.join(H, "sets_INPT_std3.json")))
v = "ext__medical_aid1_per100k"
j["forced"] = j["forced"] + [v]; j["sensitivity_only_extra"] = [x for x in j["sensitivity_only_extra"] if x != v]
j["extra_columns"][v] = ["EXT_canonical_pool_INPT.csv", "medical_aid1_per100k"]; j["transform"][v] = "none"
j["sens_note"] = "사후 민감도: 의료급여 1종 진료 실인원 10만 명당 강제 추가 — 결과에 의료급여 수급자 입원일수가 포함되는 구성 효과 점검"
fp = os.path.join(H, "sets_INPT_ma.json"); json.dump(j, open(fp, "w"), ensure_ascii=False, indent=1)
print("sets_INPT_ma.json", hashlib.sha256(open(fp, "rb").read()).hexdigest())
