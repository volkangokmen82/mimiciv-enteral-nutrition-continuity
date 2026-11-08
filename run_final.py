#!/usr/bin/env python3
"""One-shot local runner for the final MIMIC-IV EN continuity analysis release.

Runs the primary analysis, phenotype adjudication, aggregate verification,
and cross-table reproducibility checks.
All MIMIC-IV reads occur locally. Only aggregate outputs are written.
"""
from __future__ import annotations
import argparse
import csv
import importlib.util
import json
import math
import shutil
import sys
import zipfile
from datetime import datetime
from pathlib import Path


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot import {path}")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def call_main(module, argv):
    old = sys.argv[:]
    try:
        sys.argv = argv
        return module.main()
    finally:
        sys.argv = old


def read_csv(path):
    with path.open(encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f))


def num(row, key):
    try:
        return float(row[key])
    except (KeyError, TypeError, ValueError):
        return None


def close(a, b, tol=0.11):
    return a is not None and b is not None and abs(a-b) <= tol


def cross_validate(main_out: Path, adj_out: Path, report_out: Path):
    checks = []
    errors = []
    warnings = []

    def check(ok, label, detail=""):
        target = checks if ok else errors
        target.append({"check": label, "detail": detail})

    # Primary cohort flow.
    flow = {r["stage"]: int(r["n_stays"]) for r in read_csv(main_out/"00_cohort_flow.csv")}
    cohorts = {r["group"]: r for r in read_csv(main_out/"01_cohort_and_feeding_counts.csv")}
    primary = cohorts["primary"]
    primary_n = int(primary["stays"])
    en_n = int(primary["strict_EN_stays"])
    pause_n = int(primary["paused_stays"])
    check(primary_n == flow["primary_postoperative_stays"], "primary cohort reconciles to flow",
          f"primary={primary_n}; flow={flow['primary_postoperative_stays']}")
    check(0 <= pause_n <= en_n <= primary_n, "cohort nesting: pause <= EN <= primary",
          f"{pause_n} <= {en_n} <= {primary_n}")

    # Event accounting.
    event_audit = {r["scenario"]: r for r in read_csv(main_out/"21_event_timing_eligibility_audit.csv")}
    ep = event_audit["PRIMARY"]
    event_n = int(ep["analysis_n"])
    restarts = int(ep["restart_events"])
    exits = int(ep["competing_icu_exits"])
    censored = int(ep["right_censored"])
    check(event_n == restarts + exits + censored, "primary time-to-event accounting",
          "The restart, competing-event, and censoring states reconcile to the analysis denominator.")
    check(event_n == pause_n - int(ep["overlap_ambiguous_excluded_n"]) -
          int(ep["same_timestamp_pairs_n"]), "primary event denominator reconciles to pause cohort",
          "The denominator reconciles after the prespecified exclusions.")
    check(int(ep["zero_time_restarts"]) == 0, "zero-time restarts excluded from primary analysis")

    # Main-analysis and phenotype-adjudication CIF/CI must be identical for the same estimand.
    cif_main = read_csv(main_out/"08_resumption_CIF_bootstrap.csv")
    main_by = {(r["scenario"], int(r["hours_after_pause"])): r for r in cif_main}
    adj_by = {r["scenario"]: r for r in read_csv(adj_out/"25_procedure_phenotype_sensitivity.csv")}
    scenario_pairs = [
        ("PRIMARY", "all_current_primary", "primary", True),
        ("DAYS_0_1", "timing_restricted_0_1_calendar_day", "0–1 calendar-day sensitivity", False),
    ]
    event_audit_by_scenario = {r["scenario"]:r for r in read_csv(main_out/"21_event_timing_eligibility_audit.csv")}
    for main_scenario, adj_scenario, label, compare_unsuppressed_counts in scenario_pairs:
        adj_row = adj_by[adj_scenario]
        if compare_unsuppressed_counts:
            main_event = event_audit_by_scenario[main_scenario]
            check(int(adj_row["time_analysis_n"]) == int(main_event["analysis_n"]),
                  f"{label} time-analysis denominator identical across pathways",
                  "The time-analysis denominators agree.")
            check(int(adj_row["restart_events"]) == int(main_event["restart_events"]),
                  f"{label} restart-event count identical across pathways",
                  "The recorded-resumption event counts agree.")
        for h, adj_est, adj_lo, adj_hi in [
            (24, "restart_CIF_24h_pct", "CIF24_CI_low_pct", "CIF24_CI_high_pct"),
            (48, "restart_CIF_48h_pct", "CIF48_CI_low_pct", "CIF48_CI_high_pct"),
        ]:
            m = main_by[(main_scenario, h)]
            check(close(num(m,"CIF_restart_pct"), num(adj_row,adj_est)),
                  f"{label} {h}h CIF identical across analysis pathways",
                  f"main={m['CIF_restart_pct']}; adjudication={adj_row[adj_est]}")
            check(close(num(m,"bootstrap_95CI_lower_pct"), num(adj_row,adj_lo)),
                  f"{label} {h}h lower CI identical across analysis pathways",
                  f"main={m['bootstrap_95CI_lower_pct']}; adjudication={adj_row[adj_lo]}")
            check(close(num(m,"bootstrap_95CI_upper_pct"), num(adj_row,adj_hi)),
                  f"{label} {h}h upper CI identical across analysis pathways",
                  f"main={m['bootstrap_95CI_upper_pct']}; adjudication={adj_row[adj_hi]}")

    # Deterministic bootstrap and disclosure-control checks.
    reproducibility = {r["check"]: r for r in read_csv(adj_out/"31_reproducibility_and_suppression_audit.csv")}
    for name in [
        "primary_record_set_same_across_tables",
        "bootstrap_invariant_to_input_order",
        "RMST_bootstrap_invariant_to_input_order",
        "timing_table_EN_complementary_suppression",
        "timing_table_pause_complementary_suppression",
    ]:
        check(reproducibility.get(name, {}).get("status") == "PASS",
              f"adjudication check: {name}",
              "Passed." if reproducibility.get(name, {}).get("status") == "PASS" else "Missing or failed.")

    # Procedure timing strata should sum to the primary cohort.
    timing = read_csv(adj_out/"23_procedure_timing_strata.csv")
    timing_total = sum(int(r["icu_stays"]) for r in timing)
    check(timing_total == primary_n, "procedure timing strata sum to primary cohort",
          f"timing_total={timing_total}; primary={primary_n}")

    # Rate-change transition structure.
    transition = read_csv(adj_out/"26_changedose_rate_transition_audit.csv")[0]
    rate_rows = transition["changedose_rate_rows"]
    exact = transition["exact_next_segment_new_order"]
    if rate_rows not in ("<11","NA") and exact not in ("<11","NA"):
        check(int(rate_rows) == int(exact), "all observed ChangedDose/Rate rows link to an exact next segment",
              f"rows={rate_rows}; exact transitions={exact}")
    else:
        warnings.append("ChangedDose/Rate transition counts were suppressed; inspect only under authorized local access.")

    # Wilson intervals independently rechecked by the verifier.
    validation = json.loads((main_out/"VALIDATION_SUMMARY.json").read_text(encoding="utf-8"))
    if validation.get("errors"):
        errors.append({"check":"aggregate verifier errors absent",
                       "detail":"; ".join(validation["errors"])})
    else:
        checks.append({"check":"aggregate verifier has zero errors",
                       "detail":f"{len(validation.get('checks_passed',[]))} aggregate checks passed"})

    # Explicit, nonfatal author-review warnings.
    warnings.extend(validation.get("warnings", []))
    warnings.append("Clinical review of the 82 procedure codes and 36 enteral item IDs remains an author responsibility.")
    warnings.append("The locally observed ChangedDose/Rate label is structurally consistent with rate-change continuation but is not an official documented status label.")
    warnings.append("The Bolus status audit is suppressed (<11); its local meaning cannot be confirmed from the released aggregate output.")
    warnings.append("Institutional ethics determination, PhysioNet DUA compliance, author declarations, and AI disclosure must be completed by the authors.")

    payload = {
        "release": "MIMIC-IV Enteral Nutrition Continuity v3.0",
        "generated": datetime.now().isoformat(timespec="seconds"),
        "status": "PASS_ANALYTIC_CHECKS_WITH_AUTHOR_REVIEW_REQUIRED" if not errors else "FAIL",
        "checks_passed": checks,
        "errors": errors,
        "warnings": warnings,
        "primary_counts": {
            "icu_stays": primary_n,
            "strict_EN_stays": en_n,
            "paused_stays": "Masked†",
            "time_analysis_n": "Masked†",
            "restart_events": "Masked†",
            "competing_ICU_exits": "Masked†",
            "right_censored": "Masked†",
        },
        "raw_patient_level_rows_exported": False,
    }
    report_out.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")
    return payload


def zip_aggregate_outputs(main_out: Path, adj_out: Path, out: Path):
    target = out/"FINAL_AGGREGATE_OUTPUTS.zip"
    with zipfile.ZipFile(target, "w", compression=zipfile.ZIP_DEFLATED) as z:
        for folder, prefix in [(main_out,"main_analysis"),(adj_out,"phenotype_adjudication")]:
            for p in sorted(folder.iterdir()):
                if p.is_file() and p.suffix.lower() in {".csv",".json",".txt"}:
                    z.write(p, arcname=f"{prefix}/{p.name}")
        for name in ("FINAL_RELEASE_REPORT.txt","FINAL_RELEASE_VALIDATION.json"):
            p=out/name
            if p.exists():
                z.write(p, arcname=name)
    return target


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=r"D:\MIMICIV")
    ap.add_argument("--config", default=r"D:\MIMICIV\config")
    ap.add_argument("--out", default=r"D:\MIMIC_FINAL_RELEASE_V3_0")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--memory-gb", type=int, default=8)
    ap.add_argument("--skip-documents", action="store_true", help=argparse.SUPPRESS)
    args = ap.parse_args()

    repo = Path(__file__).resolve().parent
    root = Path(args.root)
    config = Path(args.config)
    out = Path(args.out)
    main_out = out/"main_analysis"
    adj_out = out/"phenotype_adjudication"
    if out.exists() and any(out.iterdir()):
        raise SystemExit(
            f"Output folder is not empty: {out}. Rename it or choose a new --out path; "
            "the final runner will not mix a new analysis with stale results."
        )
    out.mkdir(parents=True, exist_ok=True)

    required_config = ["procedure_codes.csv","formula_item_ids.csv","nursing_item_dictionary.csv"]
    missing_config = [str(config/x) for x in required_config if not (config/x).is_file()]
    if missing_config:
        raise SystemExit("Missing local config files:\n" + "\n".join(missing_config))

    print("[1/5] Running primary analysis with deterministic patient-cluster bootstrap...")
    analysis = load_module(repo/"analysis_v3.py","mimic_analysis_v30")
    call_main(analysis, [
        "analysis_v3.py","--root",str(root),"--out",str(main_out),"--config",str(config),
        "--nursing","--threads",str(args.threads),"--memory-gb",str(args.memory_gb)
    ])

    print("[2/5] Running procedure/status phenotype adjudication...")
    adjudication = load_module(repo/"phenotype_adjudication_v3.py","mimic_adjudication_v30")
    call_main(adjudication, [
        "phenotype_adjudication_v3.py","--root",str(root),"--config",str(config),
        "--analysis-script",str(repo/"analysis_v3.py"),"--out",str(adj_out),
        "--threads",str(args.threads),"--memory-gb",str(args.memory_gb)
    ])

    print("[3/5] Running aggregate validation...")
    verifier = load_module(repo/"verify_results_v3.py","mimic_verify_v30")
    audit, code = verifier.validate(main_out, config_dir=config)
    validation_path = main_out/"VALIDATION_SUMMARY.json"
    validation_path.write_text(json.dumps(audit.report(), ensure_ascii=False, indent=2), encoding="utf-8")
    if code != 0:
        raise SystemExit("Aggregate verification failed. Review main_analysis/VALIDATION_SUMMARY.json")

    print("[4/5] Cross-validating bootstrap estimates, event counts, strata, and disclosure controls...")
    payload = cross_validate(main_out, adj_out, out/"FINAL_RELEASE_VALIDATION.json")
    report_lines = [
        "MIMIC-IV ENTERAL NUTRITION CONTINUITY — FINAL RELEASE V3.0",
        f"Generated: {payload['generated']}",
        f"Status: {payload['status']}",
        "",
        "Primary cohort and event counts:",
        json.dumps(payload["primary_counts"], indent=2),
        "",
        f"Automated cross-checks passed: {len(payload['checks_passed'])}",
        f"Automated errors: {len(payload['errors'])}",
        f"Author-review warnings: {len(payload['warnings'])}",
        "",
        "This report validates internal reproducibility and aggregate consistency only.",
        "It does not establish clinical validity of all procedure or enteral item phenotypes.",
        "It does not replace local ethics approval/exemption or PhysioNet data-use compliance.",
        "No patient-level rows or identifiers are exported.",
        "The aggregate archive and raw audit CSVs are internal QA artifacts; do not submit or publish them.",
        "",
        "Warnings:",
        *[f"- {w}" for w in payload["warnings"]],
    ]
    (out/"FINAL_RELEASE_REPORT.txt").write_text("\n".join(report_lines), encoding="utf-8")

    if payload["errors"]:
        raise SystemExit("Final cross-validation failed. Review the local validation report.")

    print("[5/5] Archiving local aggregate QA outputs...")
    aggregate_zip = zip_aggregate_outputs(main_out, adj_out, out)
    print("FINAL RUN COMPLETE")
    print("Main results:", main_out)
    print("Phenotype adjudication:", adj_out)
    print("Final report:", out/"FINAL_RELEASE_REPORT.txt")
    print("Aggregate archive:", aggregate_zip)
    print("Status:", payload["status"])


if __name__ == "__main__":
    main()
