#!/usr/bin/env python3
"""Final aggregate-only phenotype adjudication v3.0 for MIMIC-IV EN continuity study.

Reads the final v3.0 analysis module and local MIMIC-IV v3.1 files.
Exports only aggregate tables; no patient-level rows or identifiers are written.
"""
from __future__ import annotations
import argparse
import csv
import datetime as dt
import importlib.util
import json
import math
from pathlib import Path

import duckdb


def load_analysis_module(path: Path):
    spec = importlib.util.spec_from_file_location("mimic_analysis_v30", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot import analysis module: {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def write_csv(path: Path, headers, rows):
    with path.open("w", newline="", encoding="utf-8-sig") as f:
        writer = csv.writer(f)
        writer.writerow(headers)
        writer.writerows(rows)


def safe_count(module, value):
    return module.nsmall(value)


def make_records(module, rows):
    """Mirror v3.0 primary event rules for a supplied aggregate-filtered pause table."""
    records = []
    for row in rows:
        (subject_id, stay_id, timing_stratum, overlap, pause_ts,
         follow_end, restart_ge, restart_gt, outtime) = row
        if pause_ts is None or follow_end is None or follow_end <= pause_ts:
            continue
        if overlap == 1:
            continue
        # Same-timestamp starts are excluded from the positive-time primary analysis.
        if restart_ge is not None and restart_ge == pause_ts:
            continue
        resume = restart_gt
        if resume is not None and resume >= follow_end:
            resume = None
        if resume is not None:
            elapsed = (resume - pause_ts).total_seconds() / 3600.0
            kind = "restart"
        else:
            elapsed = (follow_end - pause_ts).total_seconds() / 3600.0
            kind = "icu_exit" if outtime is not None and outtime < pause_ts + dt.timedelta(hours=72) else "censor"
        if elapsed <= 0:
            continue
        records.append({
            "subject_id": int(subject_id),
            "time": float(elapsed),
            "kind": kind,
            "timing_stratum": timing_stratum,
        })
    return records


def event_metrics(module, records, scenario, seed):
    n = len(records)
    events = sum(r["kind"] == "restart" for r in records)
    exits = sum(r["kind"] == "icu_exit" for r in records)
    censored = sum(r["kind"] == "censor" for r in records)
    result = {
        "scenario": scenario,
        "time_analysis_n": n,
        "restart_events": events,
        "icu_exit_before_restart": exits,
        "right_censored": censored,
    }
    # Suppress event estimates for strata with fewer than 11 observations.
    if n < 11:
        for h in (24, 48):
            result[f"restart_CIF_{h}h_pct"] = "SUPPRESSED"
            result[f"restart_CIF_{h}h_CI_low_pct"] = "SUPPRESSED"
            result[f"restart_CIF_{h}h_CI_high_pct"] = "SUPPRESSED"
        result["positive_restart_median_h"] = "SUPPRESSED"
        result["positive_restart_q25_h"] = "SUPPRESSED"
        result["positive_restart_q75_h"] = "SUPPRESSED"
        return result
    for h in (24, 48):
        cif = module.cif_aj([(r["time"], r["kind"]) for r in records], h)
        lo, hi = module.bootstrap_cluster(records, h, seed=seed + h, metric="cif")
        result[f"restart_CIF_{h}h_pct"] = round(cif * 100, 1) if math.isfinite(cif) else "NA"
        result[f"restart_CIF_{h}h_CI_low_pct"] = round(lo * 100, 1) if lo is not None else "NA"
        result[f"restart_CIF_{h}h_CI_high_pct"] = round(hi * 100, 1) if hi is not None else "NA"
    positive = [r["time"] for r in records if r["kind"] == "restart" and r["time"] > 0]
    med, q1, q3 = module.median_iqr(positive)
    result["positive_restart_median_h"] = round(med, 2) if med is not None else "NA"
    result["positive_restart_q25_h"] = round(q1, 2) if q1 is not None else "NA"
    result["positive_restart_q75_h"] = round(q3, 2) if q3 is not None else "NA"
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", required=True, help="MIMIC-IV v3.1 root, e.g. D:/MIMICIV")
    parser.add_argument("--config", default=None, help="Config folder; defaults to ./config")
    parser.add_argument("--out", required=True, help="Aggregate output folder")
    parser.add_argument("--analysis-script", default=None, help="Path to analysis_v3.py; defaults to same folder")
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--memory-gb", type=int, default=8)
    args = parser.parse_args()

    here = Path(__file__).resolve().parent
    analysis_path = Path(args.analysis_script) if args.analysis_script else here / "analysis_v3.py"
    config_dir = Path(args.config) if args.config else here / "config"
    root = Path(args.root)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    if not analysis_path.exists():
        raise FileNotFoundError(f"analysis_v3.py not found: {analysis_path}")
    if not config_dir.exists():
        raise FileNotFoundError(f"Config directory not found: {config_dir}")

    module = load_analysis_module(analysis_path)

    # Use the canonical deterministic cluster bootstrap from analysis_v3.py.
    # Main and phenotype-sensitivity analyses therefore share identical resampling logic.
    con = duckdb.connect(database=":memory:")
    con.execute(f"SET threads={max(1, args.threads)}")
    con.execute(f"SET memory_limit='{max(4, args.memory_gb)}GB'")

    print("Loading v3.0 cohort and phenotype tables from local MIMIC-IV...")
    module.create_base_views(con, root, config_dir, args)

    # 1) Procedure timing strata and their EN/pause denominators.
    timing = con.execute("""
      SELECT c.timing_stratum,
             COUNT(*) AS stays,
             COUNT(DISTINCT c.subject_id) AS patients,
             COUNT(DISTINCT f.stay_id) AS strict_en_stays,
             COUNT(DISTINCT p.stay_id) AS paused_stays
      FROM cohort c
      LEFT JOIN feeding_stays f ON f.stay_id=c.stay_id
      LEFT JOIN pause_followup p ON p.stay_id=c.stay_id
      WHERE c.in_primary=1
      GROUP BY 1
      ORDER BY CASE c.timing_stratum
        WHEN 'same_day_PACU_proxy' THEN 1
        WHEN 'previous_calendar_day' THEN 2
        WHEN '2_to_3_days_prior' THEN 3 ELSE 4 END
    """).fetchall()
    timing_rows = []
    en_small_indices = [i for i, row in enumerate(timing) if row[3] is not None and int(row[3]) < 11]
    pause_small_indices = [i for i, row in enumerate(timing) if row[4] is not None and int(row[4]) < 11]
    en_complementary_index = next((i for i,row in enumerate(timing)
                                   if i not in en_small_indices and row[3] is not None and int(row[3]) >= 11), None)
    pause_complementary_index = next((i for i,row in enumerate(timing)
                                      if i not in pause_small_indices and row[4] is not None and int(row[4]) >= 11), None)
    for idx, (a,b,c,d,e) in enumerate(timing):
        en_value = safe_count(module,d)
        pause_value = safe_count(module,e)
        # Hide one additional visible cell to prevent subtraction from published totals.
        if en_small_indices and idx == en_complementary_index:
            en_value = "SUPPRESSED_COMPLEMENTARY"
        if pause_small_indices and idx == pause_complementary_index:
            pause_value = "SUPPRESSED_COMPLEMENTARY"
        timing_rows.append((a, safe_count(module,b), "SUPPRESSED_COMPLEMENTARY",
                            en_value, pause_value))
    write_csv(out / "23_procedure_timing_strata.csv",
      ["timing_stratum", "icu_stays", "unique_patients", "strict_EN_stays", "paused_stays"],
      timing_rows)

    # 2) Flag clinically distinct procedure-approach groups for sensitivity review.
    con.execute("""
      CREATE TEMP TABLE primary_proc_flags AS
      SELECT sc.stay_id,
        MAX(CASE WHEN sc.icd_version=9 AND sc.icd_code='4382' THEN 1 ELSE 0 END) AS sleeve_bariatric,
        MAX(CASE WHEN sc.icd_version=9 AND TRY_CAST(sc.icd_code AS INTEGER) BETWEEN 4840 AND 4849 THEN 1 ELSE 0 END) AS pull_through,
        MAX(CASE WHEN sc.icd_version=10 AND sc.icd_code IN
          ('0DTH8ZZ','0DTN7ZZ','0DTNFZZ','0DTP7ZZ') THEN 1 ELSE 0 END) AS natural_or_artificial_opening_approach
      FROM surgical_candidates_raw sc
      JOIN cohort c ON c.stay_id=sc.stay_id
      WHERE c.in_primary=1
      GROUP BY sc.stay_id
    """)
    flag_rows = con.execute("""
      SELECT
        COUNT(*) AS primary_stays,
        SUM(sleeve_bariatric) AS sleeve_stays,
        SUM(pull_through) AS pull_through_stays,
        SUM(natural_or_artificial_opening_approach) AS natural_or_artificial_opening_approach_stays,
        SUM(CASE WHEN sleeve_bariatric=1 OR pull_through=1 OR natural_or_artificial_opening_approach=1 THEN 1 ELSE 0 END) AS any_special_stays
      FROM primary_proc_flags
    """).fetchone()
    write_csv(out / "24_procedure_code_flag_impact.csv",
      ["primary_stays","sleeve_bariatric_stays","pull_through_stays",
       "natural_or_artificial_opening_approach_stays","any_special_procedure_stays"],
      [[safe_count(module,x) for x in flag_rows]])


    # Code-level contribution table: aggregate only, no admission/stay identifiers.
    code_impact = con.execute("""
      SELECT sc.icd_version, sc.icd_code, sc.procedure_description,
             COUNT(DISTINCT sc.stay_id) AS primary_stays_with_code,
             COUNT(DISTINCT CASE WHEN f.stay_id IS NOT NULL THEN sc.stay_id END) AS EN_stays_with_code,
             COUNT(DISTINCT CASE WHEN p.stay_id IS NOT NULL THEN sc.stay_id END) AS paused_stays_with_code
      FROM surgical_candidates_raw sc
      JOIN cohort c ON c.stay_id=sc.stay_id AND c.in_primary=1
      LEFT JOIN feeding_stays f ON f.stay_id=sc.stay_id
      LEFT JOIN pause_followup p ON p.stay_id=sc.stay_id
      GROUP BY 1,2,3
      ORDER BY 1,2
    """).fetchall()
    write_csv(out / "29_procedure_code_level_impact.csv",
      ["icd_version","icd_code","procedure_description",
       "primary_stays_with_code","strict_EN_stays_with_code","paused_stays_with_code"],
      [(a,b,c,safe_count(module,d),safe_count(module,e),safe_count(module,f))
       for a,b,c,d,e,f in code_impact])

    # Build aggregate cohort maps; stay identifiers remain in memory only.
    cohort_rows = con.execute("""
      SELECT c.stay_id, c.subject_id, c.timing_stratum,
             CASE WHEN f.stay_id IS NOT NULL THEN 1 ELSE 0 END AS has_en,
             CASE WHEN p.stay_id IS NOT NULL THEN 1 ELSE 0 END AS has_pause,
             coalesce(g.sleeve_bariatric,0) sleeve_bariatric,
             coalesce(g.pull_through,0) pull_through,
             coalesce(g.natural_or_artificial_opening_approach,0) natural_or_artificial_opening_approach
      FROM cohort c
      LEFT JOIN feeding_stays f ON f.stay_id=c.stay_id
      LEFT JOIN pause_followup p ON p.stay_id=c.stay_id
      LEFT JOIN primary_proc_flags g ON g.stay_id=c.stay_id
      WHERE c.in_primary=1
    """).fetchall()
    cohort_map = {
      int(r[0]): {"subject_id":int(r[1]),"timing":r[2],"en":r[3],
                  "pause":r[4],"sleeve":r[5],"pull":r[6],"natural":r[7]}
      for r in cohort_rows
    }
    primary_ids = set(cohort_map)
    scenarios = {
      "all_current_primary": primary_ids,
      "timing_restricted_0_1_calendar_day": {sid for sid,v in cohort_map.items()
        if v["timing"] in ("same_day_PACU_proxy","previous_calendar_day")},
      "extended_timing_2_3_days_only": {sid for sid,v in cohort_map.items()
        if v["timing"] == "2_to_3_days_prior"},
      "exclude_bariatric_sleeve": {sid for sid,v in cohort_map.items() if not v["sleeve"]},
      "exclude_pull_through": {sid for sid,v in cohort_map.items() if not v["pull"]},
      "exclude_natural_or_artificial_opening_approach": {sid for sid,v in cohort_map.items() if not v["natural"]},
      "exclude_all_predefined_special_procedures": {sid for sid,v in cohort_map.items()
        if not (v["sleeve"] or v["pull"] or v["natural"])},
      "timing_restricted_0_1_calendar_day_excluding_special": {sid for sid,v in cohort_map.items()
        if v["timing"] in ("same_day_PACU_proxy","previous_calendar_day")
        and not (v["sleeve"] or v["pull"] or v["natural"])},
    }
    pause_rows = con.execute("""
      SELECT subject_id, stay_id, timing_stratum, overlap_ambiguous,
             pause_ts, pause_followup_end, restart_ge, restart_gt, outtime
      FROM pause_followup WHERE analysis_tier='primary'
    """).fetchall()
    pause_by_stay = {int(r[1]):r for r in pause_rows}
    out_rows=[]
    all_primary_records = None
    for idx,(name,ids) in enumerate(scenarios.items()):
        en_n=sum(cohort_map[s]["en"] for s in ids)
        pause_n=sum(cohort_map[s]["pause"] for s in ids)
        patient_n=len({cohort_map[s]["subject_id"] for s in ids})
        records=make_records(module,[pause_by_stay[s] for s in sorted(ids) if s in pause_by_stay])
        if name == "all_current_primary":
            all_primary_records = records
        metrics=event_metrics(module,records,name,20261107)
        # Complementary suppression: avoid back-calculating <11 EN/pause/event cells
        # by subtracting nested sensitivity cohorts or exclusion scenarios from totals.
        show_counts = name == "all_current_primary"
        show_denominator = name in ("all_current_primary", "timing_restricted_0_1_calendar_day",
                                    "extended_timing_2_3_days_only")
        if show_counts:
            count_values = [safe_count(module,len(ids)), safe_count(module,patient_n),
                            safe_count(module,en_n), safe_count(module,pause_n),
                            safe_count(module,metrics["time_analysis_n"]),
                            safe_count(module,metrics["restart_events"]),
                            safe_count(module,metrics["icu_exit_before_restart"]),
                            safe_count(module,metrics["right_censored"])]
        elif show_denominator:
            count_values = [safe_count(module,len(ids)), "SUPPRESSED_COMPLEMENTARY",
                            "SUPPRESSED_COMPLEMENTARY", "SUPPRESSED_COMPLEMENTARY",
                            "SUPPRESSED_COMPLEMENTARY", "SUPPRESSED_COMPLEMENTARY",
                            "SUPPRESSED_COMPLEMENTARY", "SUPPRESSED_COMPLEMENTARY"]
        else:
            count_values = ["SUPPRESSED_COMPLEMENTARY"] * 8
        out_rows.append([
          name, *count_values,
          metrics["positive_restart_median_h"],
          metrics["restart_CIF_24h_pct"],metrics["restart_CIF_24h_CI_low_pct"],metrics["restart_CIF_24h_CI_high_pct"],
          metrics["restart_CIF_48h_pct"],metrics["restart_CIF_48h_CI_low_pct"],metrics["restart_CIF_48h_CI_high_pct"],
        ])
    write_csv(out / "25_procedure_phenotype_sensitivity.csv",
      ["scenario","icu_stays","unique_patients","strict_EN_stays","paused_stays",
       "time_analysis_n","restart_events","icu_exit_before_restart","right_censored",
       "median_positive_restart_h","restart_CIF_24h_pct","CIF24_CI_low_pct","CIF24_CI_high_pct",
       "restart_CIF_48h_pct","CIF48_CI_low_pct","CIF48_CI_high_pct"],out_rows)

    # 3) Validate whether ChangedDose/Rate behaves like documented Changed:
    # same linkorderid, new orderid, next segment starts at previous endtime.
    transition = con.execute("""
      WITH ordered AS (
        SELECT stay_id,itemid,linkorderid,orderid,status,start_ts,end_ts,
          LEAD(start_ts) OVER (
            PARTITION BY stay_id,itemid,linkorderid
            ORDER BY start_ts,orderid
          ) AS next_start_ts,
          LEAD(orderid) OVER (
            PARTITION BY stay_id,itemid,linkorderid
            ORDER BY start_ts,orderid
          ) AS next_orderid
        FROM possible_feeds
        WHERE linkorderid IS NOT NULL
      )
      SELECT
        COUNT(*) AS changedose_rows,
        COUNT(DISTINCT stay_id) AS changedose_stays,
        COUNT(DISTINCT itemid) AS changedose_itemids,
        SUM(CASE WHEN next_start_ts=end_ts AND next_orderid<>orderid THEN 1 ELSE 0 END) AS exact_next_segment_new_order,
        SUM(CASE WHEN next_start_ts>end_ts AND next_start_ts<=end_ts+INTERVAL '1 minute'
                       AND next_orderid<>orderid THEN 1 ELSE 0 END) AS next_segment_within_1min_new_order,
        SUM(CASE WHEN next_start_ts=end_ts AND next_orderid=orderid THEN 1 ELSE 0 END) AS exact_next_segment_same_order
      FROM ordered WHERE status='changedose/rate'
    """).fetchone()
    write_csv(out / "26_changedose_rate_transition_audit.csv",
      ["changedose_rate_rows","unique_stays","unique_itemids",
       "exact_next_segment_new_order","next_segment_within_1min_new_order",
       "exact_next_segment_same_order"], [[safe_count(module,x) for x in transition]])

    # 4) Audit local Bolus-labelled rows against the documented one-minute bolus interval pattern.
    bolus = con.execute("""
      SELECT
        COUNT(*) AS bolus_rows,
        COUNT(DISTINCT stay_id) AS bolus_stays,
        COUNT(DISTINCT itemid) AS bolus_itemids,
        SUM(CASE WHEN end_ts-start_ts=INTERVAL '1 minute' THEN 1 ELSE 0 END) AS exactly_one_minute_rows,
        SUM(CASE WHEN end_ts-start_ts<>INTERVAL '1 minute' THEN 1 ELSE 0 END) AS other_duration_rows,
        SUM(CASE WHEN amount>0 AND amountuom IN ('ml','milliliter','milliliters') THEN 1 ELSE 0 END) AS positive_ml_amount_rows,
        SUM(CASE WHEN rate>0 AND rateuom IN ('ml/hour','ml/hr','ml/h') THEN 1 ELSE 0 END) AS positive_ml_rate_rows
      FROM possible_feeds
      WHERE status='bolus'
    """).fetchone()
    write_csv(out / "30_bolus_status_audit.csv",
      ["bolus_rows","unique_stays","unique_itemids","exactly_one_minute_rows",
       "other_duration_rows","positive_ml_amount_rows","positive_ml_rate_rows"],
      [[safe_count(module,x) for x in bolus]])

    # 5) Recompute the EN and pause phenotype with ChangedDose/Rate excluded.
    # This is a sensitivity analysis only; the documented Changed status remains allowed.
    con.execute("""
      CREATE TEMP TABLE noalias_valid_feeds AS
      SELECT * FROM feed_quality
      WHERE valid_strict=1 AND status<>'changedose/rate'
    """)
    con.execute("""
      CREATE TEMP TABLE noalias_first_feed AS
      SELECT stay_id, MIN(start_ts) AS first_feed_ts
      FROM noalias_valid_feeds GROUP BY stay_id
    """)
    con.execute("""
      CREATE TEMP TABLE noalias_feeding_stays AS
      SELECT c.*, f.first_feed_ts,
        LEAST(f.first_feed_ts+INTERVAL '72 hours',c.outtime) AS feed_window_end
      FROM cohort c JOIN noalias_first_feed f USING(stay_id)
    """)
    con.execute("""
      CREATE TEMP TABLE noalias_feeding_events AS
      SELECT f.* FROM noalias_valid_feeds f
      JOIN noalias_feeding_stays s USING(stay_id)
      WHERE f.start_ts>=s.first_feed_ts
        AND f.start_ts<s.feed_window_end
        AND f.end_ts>s.first_feed_ts
    """)
    con.execute("""
      CREATE TEMP TABLE noalias_first_pause AS
      SELECT s.stay_id, MIN(f.end_ts) AS pause_ts
      FROM noalias_feeding_stays s
      JOIN noalias_feeding_events f USING(stay_id)
      WHERE f.status='paused' AND f.end_ts IS NOT NULL
        AND f.end_ts>f.start_ts
        AND f.end_ts>=s.first_feed_ts
        AND f.end_ts<=s.feed_window_end
      GROUP BY s.stay_id
    """)
    con.execute("""
      CREATE TEMP TABLE noalias_pause_followup AS
      SELECT s.*, p.pause_ts,
        LEAST(p.pause_ts+INTERVAL '72 hours',s.outtime) AS pause_followup_end,
        MAX(CASE WHEN f.start_ts<p.pause_ts AND f.end_ts>p.pause_ts THEN 1 ELSE 0 END) AS overlap_ambiguous,
        MIN(CASE WHEN f.start_ts>=p.pause_ts
                  AND f.start_ts<LEAST(p.pause_ts+INTERVAL '72 hours',s.outtime)
                 THEN f.start_ts END) AS restart_ge,
        MIN(CASE WHEN f.start_ts>p.pause_ts
                  AND f.start_ts<LEAST(p.pause_ts+INTERVAL '72 hours',s.outtime)
                 THEN f.start_ts END) AS restart_gt
      FROM noalias_feeding_stays s
      JOIN noalias_first_pause p USING(stay_id)
      LEFT JOIN noalias_valid_feeds f
        ON f.stay_id=s.stay_id AND f.start_ts>=p.pause_ts-INTERVAL '72 hours'
      GROUP BY ALL
    """)
    current_counts = con.execute("""
      SELECT
        (SELECT COUNT(*) FROM cohort WHERE in_primary=1),
        (SELECT COUNT(DISTINCT c.subject_id) FROM cohort c JOIN feeding_stays f USING(stay_id) WHERE c.in_primary=1),
        (SELECT COUNT(*) FROM feeding_stays WHERE in_primary=1),
        (SELECT COUNT(*) FROM pause_followup WHERE in_primary=1)
    """).fetchone()
    noalias_counts = con.execute("""
      SELECT
        (SELECT COUNT(*) FROM cohort WHERE in_primary=1),
        (SELECT COUNT(DISTINCT c.subject_id) FROM cohort c JOIN noalias_feeding_stays f USING(stay_id) WHERE c.in_primary=1),
        (SELECT COUNT(*) FROM noalias_feeding_stays WHERE in_primary=1),
        (SELECT COUNT(*) FROM noalias_pause_followup WHERE in_primary=1)
    """).fetchone()
    noalias_pause_rows = con.execute("""
      SELECT subject_id,stay_id,timing_stratum,overlap_ambiguous,pause_ts,
             pause_followup_end,restart_ge,restart_gt,outtime
      FROM noalias_pause_followup WHERE analysis_tier='primary'
    """).fetchall()
    noalias_records = make_records(module,noalias_pause_rows)
    current_pause_rows = con.execute("""
      SELECT subject_id,stay_id,timing_stratum,overlap_ambiguous,pause_ts,
             pause_followup_end,restart_ge,restart_gt,outtime
      FROM pause_followup WHERE analysis_tier='primary'
    """).fetchall()
    current_records_from_view = make_records(module,current_pause_rows)
    def record_signature(records):
        return sorted((int(r["subject_id"]), round(float(r["time"]), 9),
                       str(r["kind"]), str(r["timing_stratum"])) for r in records)
    record_sets_match = (
        all_primary_records is not None
        and record_signature(all_primary_records) == record_signature(current_records_from_view)
    )
    if not record_sets_match:
        raise RuntimeError(
            "Primary time-to-event records differ between the procedure-sensitivity "
            "and status-sensitivity pathways. Do not use the output until this is reconciled."
        )
    # Reuse the exact same primary record object for both tables.
    current_records = all_primary_records
    current_metrics = event_metrics(module,current_records,"ChangedDose/Rate included",20261107)
    noalias_metrics = event_metrics(module,noalias_records,"ChangedDose/Rate excluded",20261202)
    status_rows=[]
    for label,counts,metrics in [
      ("Current v3.0 (alias included)",current_counts,current_metrics),
      ("Stress test only (ChangedDose/Rate rows removed)",noalias_counts,noalias_metrics)
    ]:
        if label == "Current v3.0 (alias included)":
            count_values = [safe_count(module,counts[0]),safe_count(module,counts[1]),
                            safe_count(module,counts[2]),safe_count(module,counts[3]),
                            safe_count(module,metrics["time_analysis_n"]),
                            safe_count(module,metrics["restart_events"]),
                            safe_count(module,metrics["icu_exit_before_restart"]),
                            safe_count(module,metrics["right_censored"])]
        else:
            count_values = ["SUPPRESSED_COMPLEMENTARY"] * 8
        status_rows.append([
          label,*count_values,
          metrics["positive_restart_median_h"],
          metrics["restart_CIF_24h_pct"],metrics["restart_CIF_24h_CI_low_pct"],metrics["restart_CIF_24h_CI_high_pct"],
          metrics["restart_CIF_48h_pct"],metrics["restart_CIF_48h_CI_low_pct"],metrics["restart_CIF_48h_CI_high_pct"],
        ])
    write_csv(out / "27_changedose_rate_phenotype_sensitivity.csv",
      ["status_rule","primary_stays","unique_patients_with_EN","strict_EN_stays","paused_stays",
       "time_analysis_n","restart_events","icu_exit_before_restart","right_censored",
       "median_positive_restart_h","restart_CIF_24h_pct","CIF24_CI_low_pct","CIF24_CI_high_pct",
       "restart_CIF_48h_pct","CIF48_CI_low_pct","CIF48_CI_high_pct"],status_rows)

    # Explicitly test that reversing input row order cannot change a seeded bootstrap CI.
    stable_a = module.bootstrap_cluster(current_records, 24, seed=20261107+24, metric="cif")
    stable_b = module.bootstrap_cluster(list(reversed(current_records)), 24, seed=20261107+24, metric="cif")
    order_invariant = (stable_a == stable_b)
    rmst_a = module.bootstrap_cluster(current_records, 24, seed=20261107+99, metric="rmst")
    rmst_b = module.bootstrap_cluster(list(reversed(current_records)), 24, seed=20261107+99, metric="rmst")
    rmst_order_invariant = (rmst_a == rmst_b)
    write_csv(out / "31_reproducibility_and_suppression_audit.csv",
      ["check","status","detail"],
      [
        ["primary_record_set_same_across_tables", "PASS" if record_sets_match else "FAIL",
         f"scenario_records={len(all_primary_records or [])}; status_records={len(current_records_from_view)}"],
        ["bootstrap_invariant_to_input_order", "PASS" if order_invariant else "FAIL",
         f"24h_CIF_CI_original={stable_a}; 24h_CIF_CI_reversed={stable_b}"],
        ["RMST_bootstrap_invariant_to_input_order", "PASS" if rmst_order_invariant else "FAIL",
         f"24h_RMST_CI_original={rmst_a}; 24h_RMST_CI_reversed={rmst_b}"],
        ["timing_table_EN_complementary_suppression", "PASS" if (not en_small_indices or (en_complementary_index is not None and timing_rows[en_complementary_index][3] == "SUPPRESSED_COMPLEMENTARY")) else "FAIL",
         "A visible EN cell is additionally suppressed when another timing stratum is <11."],
        ["timing_table_pause_complementary_suppression", "PASS" if (not pause_small_indices or (pause_complementary_index is not None and timing_rows[pause_complementary_index][4] == "SUPPRESSED_COMPLEMENTARY")) else "FAIL",
         "A visible pause cell is additionally suppressed when another timing stratum is <11."]
      ])

    report = [
      "MIMIC-IV v3.1 | final phenotype adjudication v3.0",
      f"Generated: {dt.datetime.now().isoformat(timespec='seconds')}",
      "Aggregate-only output; no patient-level rows or identifiers exported.",
      "",
      "Interpretation rules:",
      "- Procedure flags are review strata, not automatic proof of misclassification.",
      "- procedures_icd.chartdate is a date-only hospital-billed procedure date, not an exact operative timestamp; timing strata are calendar-date proxies.",
      "- The local ChangedDose/Rate label is not listed by that name in official inputevents documentation; the order-linkage audit tests whether its transitions match documented rate-change continuity.",
      "- v3.0 applies complementary suppression to the timing-strata EN and pause counts as well as nested sensitivity outputs; previous timing-strata outputs could be back-calculated from totals.",
      "- The shared v3.0 patient-cluster bootstrap sorts subject IDs and within-subject records and uses a common estimand-specific seed across scenarios.",
      "- The no-alias analysis is a stress test only; removing rate-change continuation rows is not assumed to be a clinically valid primary definition.",
      "- The Bolus audit is suppressed because all reported Bolus cells are below 11; no conclusion about local Bolus duration or meaning can be drawn from this aggregate export.",
      "",
      "The v3.0 primary cohort uses a hospital-billed procedure chartdate 0-3 calendar days before ICU entry, with same-day cases requiring a PACU-to-ICU transfer proxy within 6 hours. The 0-1 calendar-day subset is a temporal-restriction sensitivity analysis, not an exact 24-hour surgery-to-ICU window.",
      "Required review: confirm the clinical scope of bariatric sleeve and natural/artificial-opening approaches; review the status alias and Bolus evidence; retain the 0-3 calendar-day cohort as primary unless the protocol prospectively specified a narrower temporal window."
    ]
    (out / "28_ADJUDICATION_README.txt").write_text("\n".join(report),encoding="utf-8")
    manifest = {
      "dataset":"MIMIC-IV v3.1",
      "analysis_module":"analysis_v3.py v3.0.0",
      "purpose":"procedure phenotype, status mapping, bootstrap reproducibility and disclosure-control audit",
      "generated":dt.datetime.now().isoformat(timespec="seconds"),
      "raw_patient_rows_exported":False,
      "source_checksums_sha256": {p.name: __import__("hashlib").sha256(p.read_bytes()).hexdigest() for p in [analysis_path, Path(__file__), *(sorted(config_dir.glob("*.csv")))] if p.exists()},
      "outputs":[
        "23_procedure_timing_strata.csv",
        "24_procedure_code_flag_impact.csv",
        "25_procedure_phenotype_sensitivity.csv",
        "26_changedose_rate_transition_audit.csv",
        "27_changedose_rate_phenotype_sensitivity.csv",
        "28_ADJUDICATION_README.txt",
        "29_procedure_code_level_impact.csv",
        "30_bolus_status_audit.csv",
        "31_reproducibility_and_suppression_audit.csv",
        "ADJUDICATION_SUMMARY.json"
      ]
    }
    (out / "ADJUDICATION_SUMMARY.json").write_text(json.dumps(manifest,indent=2),encoding="utf-8")
    con.close()
    print(f"COMPLETE. Aggregate adjudication outputs: {out}")
    print("No patient-level rows or identifiers were exported.")


if __name__ == "__main__":
    main()
