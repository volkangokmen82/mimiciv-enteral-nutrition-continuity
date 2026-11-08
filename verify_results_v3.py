#!/usr/bin/env python3
"""Aggregate-level consistency audit for a *local* MIMIC-IV analysis run.

No individual records are loaded. Results are never uploaded anywhere.
Exit status 0 = no computed inconsistencies, including Wilson interval checks; it is NOT a phenotype validation.
Exit status 1 = a hard, internally inconsistent aggregate result.
Exit status 2 = missing required aggregate artifacts.
"""
from __future__ import annotations
import argparse
import csv
import json
import math
import hashlib
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

FLOW_STAGES = (
    'adult_ICU_stays',
    'adult_ICU_stays_with_any_curated_resection_code',
    'stays_with_curated_resection_code_0_to_7_days_prior',
    'primary_postoperative_stays',
    'primary_stays_with_strict_continuous_EN',
    'primary_EN_stays_with_first_Paused_record',
)
REQUIRED = {
    '23_first_feed_timing_summary.csv',
    '00_cohort_flow.csv','01_cohort_and_feeding_counts.csv',
    '02_procedure_code_audit.csv','03_formula_item_audit.csv',
    '04_input_unit_audit.csv','05_input_status_audit.csv',
    '06_procedure_day_difference.csv','07_baseline_characteristics.csv',
    '08_resumption_CIF_bootstrap.csv','09_resumption_followup_summary.csv',
    '10_recorded_formula_volume.csv','11_ICU_entry_recording_states.csv',
    '12_descriptive_clinical_outcomes.csv',
    '13_bedside_nursing_relevant_charting.csv','14_bedside_charting_window_denominators.csv',
    '15_VALIDATION_REPORT.txt','18_order_linkage_audit.csv',
    '19_amount_rate_consistency.csv','20_cohort_unit_audit.csv',
    '21_event_timing_eligibility_audit.csv','22_status_dictionary_review.csv',
    'manifest.json',
}

def load_csv(path: Path):
    with path.open(encoding='utf-8-sig', newline='') as f:
        return list(csv.DictReader(f))

def numeric(raw):
    try:
        z = float(raw)
        return z if math.isfinite(z) else None
    except (TypeError, ValueError):
        return None

def as_count(value):
    n = numeric(value)
    if n is None or n < 0 or n != int(n):
        raise ValueError(f'Expected nonnegative integer count; got {value!r}')
    return int(n)


def wilson_interval_pct(successes, total, z=1.959963984540054):
    if total <= 0 or successes < 0 or successes > total:
        raise ValueError('Invalid numerator/denominator for Wilson interval')
    p = successes / total
    den = 1.0 + z*z/total
    center = (p + z*z/(2.0*total))/den
    half = z*math.sqrt(p*(1.0-p)/total + z*z/(4.0*total*total))/den
    return 100.0*(center-half), 100.0*(center+half)


class Audit:
    def __init__(self):
        self.errors = []
        self.warnings = []
        self.passed = []
    def check(self, predicate, message, *, warning=False):
        (self.passed if predicate else self.warnings if warning else self.errors).append(message)
    def report(self):
        return {
            'validation_type':'Aggregate consistency only; NOT clinical phenotype validation',
            'status':'FAIL' if self.errors else 'PASS_AGGREGATE_CHECKS_WITH_UNRESOLVED_CLINICAL_REVIEW',
            'errors':self.errors,'warnings':self.warnings,
            'checks_passed':self.passed,
        }

def validate(results: Path, config_dir: Path | None = None):
    audit = Audit()
    missing = sorted(x for x in REQUIRED if not (results/x).is_file())
    if missing:
        audit.errors.append('Missing required aggregate artifacts: ' + ', '.join(missing))
        return audit,2
    flow_rows = load_csv(results/'00_cohort_flow.csv')
    flow = {r['stage']:as_count(r['n_stays']) for r in flow_rows}
    for key in FLOW_STAGES:
        audit.check(key in flow, f'Flow stage present: {key}')
    if audit.errors:
        return audit,1
    stages = FLOW_STAGES
    vals = [flow[k] for k in stages]
    audit.check(all(a >= b for a,b in zip(vals,vals[1:])),
                'Cohort flow is nonincreasing through postoperative-pause selection')
    cohorts = {r['group']:r for r in load_csv(results/'01_cohort_and_feeding_counts.csv')}
    audit.check('primary' in cohorts,'Primary cohort row exists')
    if 'primary' in cohorts:
        p = cohorts['primary']
        audit.check(as_count(p['stays']) == flow['primary_postoperative_stays'], 'Primary stays match cohort flow')
        audit.check(as_count(p['strict_EN_stays']) == flow['primary_stays_with_strict_continuous_EN'], 'EN stays match cohort flow')
        audit.check(as_count(p['paused_stays']) == flow['primary_EN_stays_with_first_Paused_record'], 'Paused stays match cohort flow')
        audit.check(as_count(p['unique_patients']) <= as_count(p['stays']), 'Unique patients do not exceed stays')
        audit.check(0 <= as_count(p['paused_stays']) <= as_count(p['strict_EN_stays']) <= as_count(p['stays']),
                    'Pause ⊆ EN ⊆ primary cohort')

    unit_rows = {r['cohort']:r for r in load_csv(results/'20_cohort_unit_audit.csv')}
    audit.check('primary' in unit_rows, 'Primary ICU-stay/admission/patient unit audit exists')
    if 'primary' in unit_rows:
        u = unit_rows['primary']
        audit.check(as_count(u['icu_stays']) == flow['primary_postoperative_stays'],
                    'Unit audit ICU-stay count matches primary cohort')
        audit.check(as_count(u['unique_patients']) == as_count(cohorts['primary']['unique_patients']),
                    'Unit audit unique-patient count matches cohort summary')
        audit.check(as_count(u['unique_admissions']) <= as_count(u['icu_stays']),
                    'Unique admissions do not exceed ICU stays')
        multiple_patients = numeric(u['patients_with_multiple_icu_stays'])
        if multiple_patients is not None:
            audit.check(multiple_patients <= as_count(u['unique_patients']),
                        'Patients with multiple ICU stays do not exceed unique patients')

    outcome_rows = {r['cohort']:r for r in load_csv(results/'12_descriptive_clinical_outcomes.csv')}
    audit.check('primary' in outcome_rows, 'Primary descriptive outcome row exists')
    if 'primary' in outcome_rows:
        o = outcome_rows['primary']
        audit.check(as_count(o['icu_stays']) == flow['primary_postoperative_stays'],
                    'Outcome ICU-stay count matches primary cohort')
        audit.check(as_count(o['unique_patients']) <= as_count(o['icu_stays']),
                    'Outcome unique patients do not exceed ICU stays')
        audit.check(as_count(o['unique_admissions']) <= as_count(o['icu_stays']),
                    'Outcome unique admissions do not exceed ICU stays')
        d = numeric(o['unique_admissions_ending_in_hospital_death'])
        pdeath = numeric(o['unique_patients_with_hospital_death'])
        if d is not None:
            audit.check(d <= as_count(o['unique_admissions']),
                        'Hospital-death admissions do not exceed unique admissions')
        if pdeath is not None:
            audit.check(pdeath <= as_count(o['unique_patients']),
                        'Patients with hospital death do not exceed unique patients')

    follow = load_csv(results/'09_resumption_followup_summary.csv')
    follow_by = {r['scenario']:r for r in follow}
    expected_scenarios = {'PRIMARY','PRIMARY_same_timestamp','DAYS_0_1',
                          'DAYS_0_1_same_timestamp','PRIMARY_legacy_first_feed_window',
                          'PRIMARY_legacy_same_timestamp'}
    audit.check(expected_scenarios.issubset(set(follow_by)),
                'All six prespecified resumption scenarios are present')
    for scen,r in follow_by.items():
        n = as_count(r['analysis_n'])
        total = sum(as_count(r[c]) for c in ('restart_events','ICU_exit_competing_events','right_censored'))
        audit.check(total == n, f'{scen}: events + exits + censoring equals analysis n')
        audit.check(as_count(r['unique_patients']) <= n, f'{scen}: unique patients <= n')
        for c in ('24h_mean_time_without_recorded_restart_h','bootstrap_95CI_lower_h','bootstrap_95CI_upper_h'):
            if numeric(r[c]) is not None:
                audit.check(0 <= numeric(r[c]) <= 24.01, f'{scen}: {c} lies within 0–24 hours')
        if scen == 'PRIMARY':
            audit.check(n <= flow['primary_EN_stays_with_first_Paused_record'], 'Primary time-analysis n <= pause cohort')

    cif = load_csv(results/'08_resumption_CIF_bootstrap.csv')
    grouped = defaultdict(list)
    for r in cif: grouped[r['scenario']].append(r)
    for scen,rs in grouped.items():
        last = -1
        expected_n = as_count(follow_by[scen]['analysis_n']) if scen in follow_by else None
        for row in sorted(rs,key=lambda r:as_count(r['hours_after_pause'])):
            h=as_count(row['hours_after_pause']); est=numeric(row['CIF_restart_pct'])
            if est is not None:
                audit.check(0 <= est <= 100, f'{scen} {h}h: CIF is bounded')
                audit.check(est >= last - .11, f'{scen} {h}h: CIF is nondecreasing')
                last=est
            if expected_n is not None:
                audit.check(as_count(row['analysis_n']) == expected_n, f'{scen} {h}h: denominator matches follow-up')
                audit.check(as_count(row['at_risk_icu_stays_n']) <= expected_n, f'{scen} {h}h: at-risk count bounded')
            lo=numeric(row['bootstrap_95CI_lower_pct']); hi=numeric(row['bootstrap_95CI_upper_pct'])
            if lo is not None and hi is not None:
                audit.check(0 <= lo <= hi <= 100.01, f'{scen} {h}h: confidence interval bounded')

    volumes=load_csv(results/'10_recorded_formula_volume.csv')
    by_cohort=defaultdict(dict)
    for r in volumes:
        group=r['cohort']; h=as_count(r['complete_window_hours_from_first_EN']); n=as_count(r['n_complete_stays'])
        by_cohort[group][h]=n
        audit.check(n >= 0 and (group != 'primary' or n <= flow['primary_stays_with_strict_continuous_EN']),
                    f'{group} {h}h: complete-window count bounded')
        q=[numeric(r.get(x)) for x in ('q25_mL','median_mL','q75_mL')]
        if all(x is not None for x in q):
            audit.check(0 <= q[0] <= q[1] <= q[2], f'{group} {h}h: volume quartiles are ordered')
    for group,counts in by_cohort.items():
        if 24 in counts and 48 in counts:
            audit.check(counts[48] <= counts[24], f'{group}: 48h-complete subset <= 24h-complete subset')

    first_feed_summary = load_csv(results/'23_first_feed_timing_summary.csv')
    audit.check(len(first_feed_summary) == 1, 'First-feed timing summary contains one primary-cohort row')
    if first_feed_summary:
        ff = first_feed_summary[0]
        ff_n = as_count(ff['primary_EN_stays'])
        audit.check(ff_n == as_count(cohorts['primary']['strict_EN_stays']),
                    'First-feed timing denominator matches primary EN stays')
        ffq = [numeric(ff.get(k)) for k in ('q25_hours','median_hours_from_ICU_entry','q75_hours')]
        if all(x is not None for x in ffq):
            audit.check(0 <= ffq[0] <= ffq[1] <= ffq[2],
                        'First-feed timing quartiles are ordered and nonnegative')

    states=load_csv(results/'11_ICU_entry_recording_states.csv')
    for r in states:
        h=as_count(r['hours_after_ICU_entry'])
        denom=as_count(r['N_stays']); en=as_count(r['EN_recorded_n'])
        out=as_count(r['ICU_exit_without_EN_record_n']); staying=as_count(r['still_in_ICU_without_EN_record_n'])
        audit.check(denom==flow['primary_postoperative_stays'], f'{h}h: ICU-entry denominator matches cohort')
        audit.check(en+out+staying==denom, f'{h}h: mutually exclusive ICU-entry states sum to N')
        observed_pct = 100*en/denom
        audit.check(abs(observed_pct - float(r['EN_recorded_pct'])) < .11,
                    f'{h}h: published percentage matches numerator/denominator')
        expected_low, expected_high = wilson_interval_pct(en, denom)
        observed_low = numeric(r['Wilson_95CI_low_pct'])
        observed_high = numeric(r['Wilson_95CI_high_pct'])
        audit.check(observed_low is not None and abs(observed_low-expected_low) < .11,
                    f'{h}h: Wilson lower confidence limit matches exact calculation')
        audit.check(observed_high is not None and abs(observed_high-expected_high) < .11,
                    f'{h}h: Wilson upper confidence limit matches exact calculation')
        if observed_low is not None and observed_high is not None:
            audit.check(observed_low <= observed_pct <= observed_high,
                        f'{h}h: Wilson interval contains observed proportion')

    event_audit = {r['scenario']:r for r in load_csv(results/'21_event_timing_eligibility_audit.csv')}
    audit.check('PRIMARY' in event_audit and 'PRIMARY_same_timestamp' in event_audit,
                'Primary and same-timestamp event eligibility audits exist')
    if 'PRIMARY' in event_audit and 'PRIMARY_same_timestamp' in event_audit:
        ep = event_audit['PRIMARY']; es = event_audit['PRIMARY_same_timestamp']
        for label, row in [('PRIMARY',ep),('PRIMARY_same_timestamp',es)]:
            n = as_count(row['analysis_n'])
            parts = sum(as_count(row[c]) for c in ('restart_events','competing_icu_exits','right_censored'))
            audit.check(n == parts, f'{label}: event/exits/censoring sum to analysis n')
        audit.check(as_count(es['analysis_n']) >= as_count(ep['analysis_n']),
                    'Same-timestamp sensitivity includes at least as many stays as primary')
        audit.check(as_count(ep['zero_time_restarts']) == 0,
                    'Primary event analysis excludes zero-time restarts')
        audit.check(as_count(es['zero_time_restarts']) == as_count(es['same_timestamp_pairs_n']),
                    'Same-timestamp sensitivity explicitly counts zero-time pairs')
        pause_n = as_count(ep['primary_pause_cohort_n'])
        overlap_n = as_count(ep['overlap_ambiguous_excluded_n'])
        same_n = as_count(ep['same_timestamp_pairs_n'])
        audit.check(as_count(ep['analysis_n']) == pause_n-overlap_n-same_n,
                    'Primary analysis denominator reconciles to pause cohort minus overlap and same-time exclusions')
        audit.check(as_count(es['analysis_n']) == pause_n-overlap_n,
                    'Same-timestamp sensitivity denominator reconciles to pause cohort minus overlap exclusions')

    status_review = load_csv(results/'22_status_dictionary_review.csv')
    nonstandard = [r for r in status_review if r['documentation_review_flag'] != 'DOCUMENTED_STATUS']
    if nonstandard:
        audit.warnings.append('Observed local statusdescription labels require manual reconciliation with the exact MIMIC-IV source version: ' +
                              ', '.join(sorted({r['observed_statusdescription'] for r in nonstandard})))
    if (results/'13_bedside_nursing_relevant_charting.csv').is_file():
        chart=load_csv(results/'13_bedside_nursing_relevant_charting.csv')
        for r in chart:
            n=r['stays_with_any_record_suppressed_lt11']
            d=as_count(r['eligible_stays_complete_window'])
            if n not in ('<11','NA',''):
                v=as_count(n)
                audit.check(v==0 or v>=11,f'{r["domain"]}/{r["window"]}: small-cell suppression')
                audit.check(v<=d,f'{r["domain"]}/{r["window"]}: numerator <= denominator')
    # Confirm that the analysis source, verifier, and configuration used for this run
    # match the files currently present beside this verifier.
    manifest = json.loads((results/'manifest.json').read_text(encoding='utf-8'))
    audit.check(manifest.get('dataset') == 'MIMIC-IV v3.1',
                'Manifest identifies MIMIC-IV v3.1')
    audit.check(manifest.get('analysis_version') == '3.0.0-final-candidate',
                'Manifest identifies final deterministic analysis release v3.0')
    expected_hashes = manifest.get('source_checksums_sha256', {})
    for filename in ('analysis_v3.py','verify_results_v3.py','phenotype_adjudication_v3.py','run_final.py','requirements.txt'):
        local_path = Path(__file__).resolve().with_name(filename)
        audit.check(filename in expected_hashes,
                    f'Run manifest contains source checksum: {filename}')
        audit.check(local_path.is_file(),
                    f'Release source file exists: {filename}')
        if filename in expected_hashes and local_path.is_file():
            digest = hashlib.sha256(local_path.read_bytes()).hexdigest()
            audit.check(digest == expected_hashes[filename],
                        f'Source checksum matches run manifest: {filename}')
    local_config = config_dir if config_dir is not None else Path(__file__).resolve().parent/'config'
    for filename in ('procedure_codes.csv','formula_item_ids.csv','nursing_item_dictionary.csv'):
        local_path = local_config/filename
        audit.check(filename in expected_hashes,
                    f'Run manifest contains configuration checksum: {filename}')
        audit.check(local_path.is_file(),
                    f'Local configuration file exists: {filename}')
        if filename in expected_hashes and local_path.is_file():
            digest = hashlib.sha256(local_path.read_bytes()).hexdigest()
            audit.check(digest == expected_hashes[filename],
                        f'Configuration checksum matches run manifest: {filename}')

    audit.warnings.append('Mandatory manual review: 82 surgical-code candidates, 36 EN item candidates, units/statuses, chart-domain definitions, and study ethics.')
    return audit, (1 if audit.errors else 0)

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--results',required=True,help='Local folder with aggregate CSVs from analysis_v3.py')
    p.add_argument('--config',default=None,help='Local config folder; not included in public repository')
    args=p.parse_args()
    folder=Path(args.results).expanduser().resolve()
    if not folder.is_dir(): p.error('Results folder does not exist')
    config_dir = Path(args.config).expanduser().resolve() if args.config else None
    audit,code=validate(folder, config_dir=config_dir)
    payload=audit.report()
    report=folder/'VALIDATION_SUMMARY.json'
    report.write_text(json.dumps(payload,ensure_ascii=False,indent=2),encoding='utf-8')
    print('VALIDATION:',payload['status'])
    print('Checks passed:',len(payload['checks_passed']),'warnings:',len(payload['warnings']),'errors:',len(payload['errors']))
    for typ in ('errors','warnings'):
        for v in payload[typ]: print(f'{typ.upper()}: {v}')
    print('Aggregate audit report:',report)
    return code

if __name__=='__main__':
    raise SystemExit(main())
