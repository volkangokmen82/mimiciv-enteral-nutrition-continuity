#!/usr/bin/env python3
"""Final reproducible, aggregate-only MIMIC-IV v3.1 analysis for postoperative EN continuity.

This script is intentionally conservative. It uses explicit, versioned code/item lists,
separates the 72-h feeding-eligibility window from the 72-h post-pause outcome window,
uses only documented administration intervals, and does not export patient-level rows.

IMPORTANT: exploratory analysis code; author/clinical phenotype validation is still required.
Public repository contains NO MIMIC data and must not be considered a confirmed reproduction.
"""
from __future__ import annotations
import argparse, csv, datetime as dt, json, math, random, hashlib, sys
from pathlib import Path
from collections import defaultdict, Counter

try:
    import duckdb
except ImportError:
    print('DuckDB missing. Install with: py -3 -m pip install duckdb', file=sys.stderr)
    raise SystemExit(2)

# ChangedDose/Rate is a locally observed, undocumented label provisionally retained as a Changed alias.
ANALYSIS_VERSION = '3.0.0-final-candidate'
ALLOWED_STATUS = {'changed', 'changedose/rate', 'finishedrunning', 'paused', 'stopped'}
RATE_UNITS = {'ml/hour', 'ml/hr', 'ml/h'}
AMOUNT_UNITS = {'ml', 'milliliter', 'milliliters'}
MIN_CELL = 11
BOOTSTRAP_REPS = 3000
RNG_SEED = 20261107

FILES = {
    'icustays': 'icu/icustays.csv.gz',
    'd_items': 'icu/d_items.csv.gz',
    'inputevents': 'icu/inputevents.csv.gz',
    'patients': 'hosp/patients.csv.gz',
    'admissions': 'hosp/admissions.csv.gz',
    'procedures_icd': 'hosp/procedures_icd.csv.gz',
    'd_icd_procedures': 'hosp/d_icd_procedures.csv.gz',
    'transfers': 'hosp/transfers.csv.gz',
    'chartevents': 'icu/chartevents.csv.gz',
    'outputevents': 'icu/outputevents.csv.gz',
}


def nsmall(x):
    if x is None: return 'NA'
    try: x = int(x)
    except Exception: return x
    return '<11' if 0 < x < MIN_CELL else x


def write_csv(path: Path, headers, rows):
    with path.open('w', newline='', encoding='utf-8-sig') as f:
        w = csv.writer(f); w.writerow(headers); w.writerows(rows)


def quantile(values, p):
    a = sorted(float(x) for x in values if x is not None and math.isfinite(float(x)))
    if not a: return None
    pos = (len(a)-1)*p; lo = int(math.floor(pos)); hi = min(len(a)-1, lo+1)
    return a[lo] + (pos-lo)*(a[hi]-a[lo])


def median_iqr(values):
    return quantile(values,.5), quantile(values,.25), quantile(values,.75)


def wilson_interval_pct(successes, total, z=1.959963984540054):
    """Two-sided Wilson score interval, returned in percentage points."""
    if total <= 0 or successes < 0 or successes > total:
        return (float('nan'), float('nan'))
    p = successes / total
    den = 1.0 + (z * z / total)
    center = (p + z * z / (2.0 * total)) / den
    half = z * math.sqrt((p * (1.0 - p) / total) +
                         (z * z / (4.0 * total * total))) / den
    return (100.0 * (center - half), 100.0 * (center + half))


def cif_aj(records, horizon):
    """Aalen-Johansen CIF for restart; record=(hours, kind), kind restart/icu_exit/censor."""
    n = len(records)
    if n == 0: return float('nan')
    risk = n; surv = 1.0; cif = 0.0
    times = sorted(set(t for t,k in records if t <= horizon))
    for t in times:
        kinds = [k for tt,k in records if tt == t]
        r = kinds.count('restart'); c = kinds.count('icu_exit'); cens = kinds.count('censor')
        if risk <= 0: break
        cif += surv * r / risk
        surv *= (1.0 - (r+c)/risk)
        risk -= (r+c+cens)
    return cif


def at_risk_n(records, horizon):
    # Number of ICU-stay records still under observation immediately before horizon.
    return sum(1 for t,k in records if t >= horizon)


def rmst_without_restart(records, horizon):
    """Integral of 1-CIF_restart(t), i.e. restricted mean time without recorded ICU restart.
    If there are no competing exits this equals conventional restart-free RMST.
    """
    times = sorted(set([0.0, horizon] + [t for t,k in records if 0 < t < horizon]))
    area = 0.0
    for a,b in zip(times[:-1], times[1:]):
        area += (b-a) * (1.0-cif_aj(records, a))
    return area


def bootstrap_cluster(records, horizon, reps=BOOTSTRAP_REPS, seed=RNG_SEED, metric='cif'):
    """Order-invariant patient-cluster percentile bootstrap.

    Subject IDs and within-subject event records are sorted before sampling.
    Identical analytic records therefore produce identical intervals regardless
    of SQL row order, Python set iteration order, or pathway used to assemble them.
    All records belonging to a sampled subject are retained as one cluster.
    """
    by_subject = defaultdict(list)
    for r in records:
        by_subject[int(r['subject_id'])].append((float(r['time']), str(r['kind'])))
    subjects = sorted(by_subject)
    if len(subjects) < 2:
        return (None, None)
    for sid in subjects:
        by_subject[sid] = sorted(by_subject[sid], key=lambda x: (x[0], x[1]))
    rng = random.Random(seed)
    vals = []
    for _ in range(reps):
        sampled = [rng.choice(subjects) for _ in subjects]
        boot = []
        for sid in sampled:
            boot.extend(by_subject[sid])
        v = cif_aj(boot, horizon) if metric == 'cif' else rmst_without_restart(boot, horizon)
        if math.isfinite(v):
            vals.append(v)
    return quantile(vals, .025), quantile(vals, .975)


def normalize_status(x):
    return (x or '').strip().lower().replace(' ', '')


def procedure_group(version, code):
    c = str(code).upper()
    if int(version)==10 and len(c)>=4 and c.startswith('0DT'):
        body = c[3]
        if body in {'6','7'}: return 'Stomach/pylorus'
        if body in {'8','9','A','B','C'}: return 'Small bowel/duodenum'
        if body in {'E','F','G','H','K','L','M','N','P'}: return 'Colon/rectum'
    if int(version)==9:
        try: n=int(c)
        except: return 'Other/uncategorized'
        if 435 <= n <= 4399: return 'Stomach/pylorus'
        if n in {4561,4562,4563,4602}: return 'Small bowel/duodenum'
        if 1731 <= n <= 1739 or 4571 <= n <= 4583 or n == 458 or 4840 <= n <= 4869 or n in {485,4850,4851,4852,4859}: return 'Colon/rectum'
    return 'Other/uncategorized'


def read_csv_view(con, root, name, path=None):
    p = root / (path or FILES[name])
    if not p.exists(): raise FileNotFoundError(str(p))
    psql = p.as_posix().replace("'", "''")
    con.execute(f"CREATE OR REPLACE VIEW v_{name} AS SELECT * FROM read_csv_auto('{psql}', HEADER=TRUE, ALL_VARCHAR=TRUE, COMPRESSION='gzip', IGNORE_ERRORS=FALSE)")


def load_config(config_dir):
    def read_rows(name):
        with (config_dir/name).open(encoding='utf-8-sig', newline='') as f: return list(csv.DictReader(f))
    codes=[(int(r['icd_version']),r['icd_code'].strip().upper()) for r in read_rows('procedure_codes.csv')]
    formulas=[(int(r['itemid']),r['candidate_label']) for r in read_rows('formula_item_ids.csv')]
    nurses=[(int(r['itemid']),r['label'],r['domain'],r['source_table']) for r in read_rows('nursing_item_dictionary.csv')]
    if len(codes)!=82 or len(set(codes))!=82: raise ValueError('Procedure code list must contain 82 unique code/version pairs.')
    if len(formulas)!=36 or len(set(x[0] for x in formulas))!=36: raise ValueError('Formula candidate list must contain 36 unique item IDs.')
    return codes,formulas,nurses


def create_base_views(con, root, config_dir, args):
    for name in ['icustays','patients','admissions','procedures_icd','d_icd_procedures','transfers','d_items','inputevents']:
        read_csv_view(con,root,name)
    codes,formula_items,nursing_items=load_config(config_dir)
    con.execute('CREATE TEMP TABLE code_config(icd_version INTEGER, icd_code VARCHAR)')
    con.executemany('INSERT INTO code_config VALUES (?,?)', codes)
    con.execute('CREATE TEMP TABLE formula_config(itemid BIGINT, candidate_label VARCHAR)')
    con.executemany('INSERT INTO formula_config VALUES (?,?)', formula_items)
    con.execute('CREATE TEMP TABLE nursing_config(itemid BIGINT, label VARCHAR, domain VARCHAR, source_table VARCHAR)')
    con.executemany('INSERT INTO nursing_config VALUES (?,?,?,?)', nursing_items)

    con.execute("""
      CREATE TEMP TABLE adult_icu AS
      SELECT TRY_CAST(i.subject_id AS BIGINT) subject_id,
             TRY_CAST(i.hadm_id AS BIGINT) hadm_id,
             TRY_CAST(i.stay_id AS BIGINT) stay_id,
             TRY_CAST(i.intime AS TIMESTAMP) intime,
             TRY_CAST(i.outtime AS TIMESTAMP) outtime,
             i.first_careunit,
             CASE WHEN TRY_CAST(p.anchor_age AS INTEGER)=91 THEN 91
                  ELSE TRY_CAST(p.anchor_age AS INTEGER) + YEAR(TRY_CAST(i.intime AS TIMESTAMP))-TRY_CAST(p.anchor_year AS INTEGER) END AS age,
             p.gender,
             a.admission_type,
             TRY_CAST(a.admittime AS TIMESTAMP) admittime,
             TRY_CAST(a.dischtime AS TIMESTAMP) dischtime,
             TRY_CAST(a.deathtime AS TIMESTAMP) deathtime,
             TRY_CAST(a.hospital_expire_flag AS INTEGER) hospital_expire_flag
      FROM v_icustays i
      JOIN v_patients p ON TRY_CAST(p.subject_id AS BIGINT)=TRY_CAST(i.subject_id AS BIGINT)
      JOIN v_admissions a ON TRY_CAST(a.hadm_id AS BIGINT)=TRY_CAST(i.hadm_id AS BIGINT)
      WHERE CASE WHEN TRY_CAST(p.anchor_age AS INTEGER)=91 THEN 91
                 ELSE TRY_CAST(p.anchor_age AS INTEGER)+YEAR(TRY_CAST(i.intime AS TIMESTAMP))-TRY_CAST(p.anchor_year AS INTEGER) END >= 18
        AND i.hadm_id IS NOT NULL
    """)
    con.execute("""
      CREATE TEMP TABLE proc_events AS
      SELECT TRY_CAST(hadm_id AS BIGINT) hadm_id, upper(trim(icd_code)) icd_code,
             TRY_CAST(icd_version AS INTEGER) icd_version, TRY_CAST(chartdate AS DATE) chartdate
      FROM v_procedures_icd WHERE hadm_id IS NOT NULL AND icd_code IS NOT NULL AND chartdate IS NOT NULL
    """)
    con.execute("""
      CREATE TEMP TABLE pacu_transfers AS
      SELECT TRY_CAST(hadm_id AS BIGINT) hadm_id, TRY_CAST(outtime AS TIMESTAMP) pacu_out
      FROM v_transfers
      WHERE (lower(coalesce(careunit,'')) LIKE '%pacu%'
         OR lower(coalesce(careunit,'')) LIKE '%post%anesthe%')
        AND outtime IS NOT NULL
    """)
    con.execute("""
      CREATE TEMP TABLE surgical_candidates_raw AS
      SELECT s.*, p.icd_code, p.icd_version, p.chartdate,
             d.long_title AS procedure_description,
             date_diff('day', p.chartdate, CAST(s.intime AS DATE)) AS surgical_day_delta,
             EXISTS (SELECT 1 FROM pacu_transfers t
                WHERE t.hadm_id=s.hadm_id AND CAST(t.pacu_out AS DATE)=p.chartdate
                  AND t.pacu_out<=s.intime AND t.pacu_out>=s.intime-INTERVAL '6 hours') AS pacu_before_icu,
             CASE WHEN p.icd_version=10 AND substr(p.icd_code,1,3)='0DT' THEN 0
                  WHEN p.icd_version=9 THEN 0 ELSE 1 END AS code_review_flag,
             CASE WHEN p.icd_version=10 AND substr(p.icd_code,1,3)='0DT' AND substr(p.icd_code,4,1) IN ('6','7') THEN 'Stomach/pylorus'
                  WHEN p.icd_version=10 AND substr(p.icd_code,1,3)='0DT' AND substr(p.icd_code,4,1) IN ('8','9','A','B','C') THEN 'Small bowel/duodenum'
                  WHEN p.icd_version=10 AND substr(p.icd_code,1,3)='0DT' AND substr(p.icd_code,4,1) IN ('E','F','G','H','K','L','M','N','P') THEN 'Colon/rectum'
                  WHEN p.icd_version=9 AND TRY_CAST(p.icd_code AS INTEGER) BETWEEN 435 AND 4399 THEN 'Stomach/pylorus'
                  WHEN p.icd_version=9 AND TRY_CAST(p.icd_code AS INTEGER) IN (4561,4562,4563,4602) THEN 'Small bowel/duodenum'
                  WHEN p.icd_version=9 AND (TRY_CAST(p.icd_code AS INTEGER) BETWEEN 1731 AND 1739 OR TRY_CAST(p.icd_code AS INTEGER) BETWEEN 4571 AND 4583 OR TRY_CAST(p.icd_code AS INTEGER)=458 OR TRY_CAST(p.icd_code AS INTEGER) BETWEEN 4840 AND 4869 OR TRY_CAST(p.icd_code AS INTEGER) IN (485,4850,4851,4852,4859)) THEN 'Colon/rectum'
                  ELSE 'Other/uncategorized' END AS procedure_group
      FROM adult_icu s
      JOIN proc_events p ON p.hadm_id=s.hadm_id
      JOIN code_config c ON c.icd_code=p.icd_code AND c.icd_version=p.icd_version
      LEFT JOIN v_d_icd_procedures d ON upper(trim(d.icd_code))=p.icd_code AND TRY_CAST(d.icd_version AS INTEGER)=p.icd_version
      WHERE date_diff('day',p.chartdate,CAST(s.intime AS DATE)) BETWEEN 0 AND 7
        AND NOT (
          date_diff('day',p.chartdate,CAST(s.intime AS DATE)) BETWEEN 1 AND 7
          AND EXISTS (SELECT 1 FROM adult_icu prev
             WHERE prev.hadm_id=s.hadm_id AND prev.stay_id<>s.stay_id
               AND prev.intime<s.intime AND prev.outtime>CAST(p.chartdate AS TIMESTAMP))
        )
    """)
    con.execute("""
      CREATE TEMP TABLE target_stays AS
      WITH ranked AS (
        SELECT *, ROW_NUMBER() OVER(PARTITION BY stay_id ORDER BY
          CASE WHEN surgical_day_delta=0 AND pacu_before_icu THEN 0
               WHEN surgical_day_delta=1 THEN 1
               WHEN surgical_day_delta BETWEEN 2 AND 3 THEN 2
               WHEN surgical_day_delta BETWEEN 4 AND 7 THEN 3
               WHEN surgical_day_delta=0 THEN 4 ELSE 5 END,
          surgical_day_delta, chartdate DESC, icd_code) rnk
        FROM surgical_candidates_raw
      ), best AS (
        SELECT * EXCLUDE(rnk) FROM ranked WHERE rnk=1
      ), groups AS (
        SELECT stay_id, string_agg(DISTINCT procedure_group, '; ' ORDER BY procedure_group) AS procedure_groups
        FROM surgical_candidates_raw GROUP BY stay_id
      )
      SELECT b.* EXCLUDE(procedure_group), g.procedure_groups FROM best b JOIN groups g USING(stay_id)
    """)
    con.execute("""
      CREATE TEMP TABLE cohort AS
      SELECT *,
        CASE WHEN surgical_day_delta=0 AND pacu_before_icu THEN 'same_day_PACU_proxy'
             WHEN surgical_day_delta=1 THEN 'previous_calendar_day'
             WHEN surgical_day_delta BETWEEN 2 AND 3 THEN '2_to_3_days_prior'
             WHEN surgical_day_delta BETWEEN 4 AND 7 THEN '4_to_7_days_prior'
             ELSE 'same_day_unverified' END AS timing_stratum,
        procedure_groups AS procedure_group,
        CASE WHEN (surgical_day_delta=0 AND pacu_before_icu) OR surgical_day_delta BETWEEN 1 AND 3
             THEN 'primary' ELSE 'sensitivity_only' END AS analysis_tier,
        CASE WHEN (surgical_day_delta=0 AND pacu_before_icu) OR surgical_day_delta BETWEEN 1 AND 3
             THEN 1 ELSE 0 END AS in_primary
      FROM target_stays
    """)

    # Formula phenotype: exact candidate item IDs, explicit accepted administration statuses and units.
    con.execute("""
      CREATE TEMP TABLE formula_dictionary AS
      SELECT f.itemid, f.candidate_label, d.label AS dictionary_label,
             d.category, d.linksto,
             CASE WHEN lower(coalesce(d.category,''))='nutrition - enteral'
                       AND lower(coalesce(d.linksto,''))='inputevents'
                  THEN 1 ELSE 0 END AS dictionary_match
      FROM formula_config f LEFT JOIN v_d_items d ON TRY_CAST(d.itemid AS BIGINT)=f.itemid
    """)
    con.execute("""
      CREATE TEMP TABLE possible_feeds AS
      SELECT TRY_CAST(e.stay_id AS BIGINT) stay_id,
             TRY_CAST(e.subject_id AS BIGINT) subject_id,
             TRY_CAST(e.itemid AS BIGINT) itemid,
             TRY_CAST(e.starttime AS TIMESTAMP) start_ts,
             TRY_CAST(e.endtime AS TIMESTAMP) end_ts,
             TRY_CAST(e.amount AS DOUBLE) amount,
             lower(trim(coalesce(e.amountuom,''))) amountuom,
             TRY_CAST(e.rate AS DOUBLE) rate,
             lower(trim(coalesce(e.rateuom,''))) rateuom,
             lower(replace(trim(coalesce(e.statusdescription,'')),' ','')) status,
             TRY_CAST(e.orderid AS BIGINT) orderid,
             TRY_CAST(e.linkorderid AS BIGINT) linkorderid,
             lower(coalesce(e.ordercategoryname,'')) order_category,
             lower(coalesce(e.ordercategorydescription,'')) order_category_description,
             d.dictionary_label, d.category, d.dictionary_match,
             c.intime icu_in, c.outtime icu_out
      FROM v_inputevents e
      JOIN cohort c ON TRY_CAST(e.stay_id AS BIGINT)=c.stay_id
      JOIN formula_dictionary d ON TRY_CAST(e.itemid AS BIGINT)=d.itemid
    """)
    con.execute("""
      CREATE TEMP TABLE feed_quality AS
      SELECT *,
        CASE WHEN dictionary_match=1
          AND status IN ('changed','changedose/rate','finishedrunning','paused','stopped')
          AND amount>0 AND rate>0
          AND amountuom IN ('ml','milliliter','milliliters')
          AND rateuom IN ('ml/hour','ml/hr','ml/h')
          AND start_ts IS NOT NULL AND end_ts>start_ts
          AND start_ts>=icu_in AND start_ts<icu_out
          THEN 1 ELSE 0 END AS valid_strict
      FROM possible_feeds
    """)
    con.execute('CREATE TEMP TABLE valid_feeds AS SELECT * FROM feed_quality WHERE valid_strict=1')
    con.execute("""
      CREATE TEMP TABLE first_feed AS
      SELECT stay_id, MIN(start_ts) first_feed_ts FROM valid_feeds GROUP BY stay_id
    """)
    con.execute("""
      CREATE TEMP TABLE feeding_stays AS
      SELECT c.*, f.first_feed_ts,
             LEAST(f.first_feed_ts+INTERVAL '72 hours',c.outtime) AS feed_window_end
      FROM cohort c JOIN first_feed f USING(stay_id)
    """)
    con.execute("""
      CREATE TEMP TABLE feeding_events AS
      SELECT f.* FROM valid_feeds f JOIN feeding_stays s USING(stay_id)
      WHERE f.start_ts>=s.first_feed_ts AND f.start_ts<s.feed_window_end AND f.end_ts>s.first_feed_ts
    """)
    con.execute("""
      CREATE TEMP TABLE first_pause AS
      SELECT s.stay_id, MIN(f.end_ts) pause_ts
      FROM feeding_stays s JOIN feeding_events f USING(stay_id)
      WHERE f.status='paused' AND f.end_ts IS NOT NULL
        AND f.end_ts>f.start_ts AND f.end_ts>=s.first_feed_ts AND f.end_ts<=s.feed_window_end
      GROUP BY s.stay_id
    """)
    # Primary event follow-up is 72 h after the pause (or earlier ICU exit), not the residual part
    # of a 72-h window measured from first feeding. Legacy first-feed-anchored follow-up is retained
    # as a sensitivity analysis so the change is transparent.
    con.execute("""
      CREATE TEMP TABLE pause_followup AS
      SELECT s.*, p.pause_ts,
        LEAST(p.pause_ts+INTERVAL '72 hours',s.outtime) AS pause_followup_end,
        LEAST(s.first_feed_ts+INTERVAL '72 hours',s.outtime) AS legacy_followup_end,
        MAX(CASE WHEN f.start_ts<p.pause_ts AND f.end_ts>p.pause_ts THEN 1 ELSE 0 END) AS overlap_ambiguous,
        MIN(CASE WHEN f.start_ts>=p.pause_ts AND f.start_ts<LEAST(p.pause_ts+INTERVAL '72 hours',s.outtime)
                 THEN f.start_ts END) AS restart_ge,
        MIN(CASE WHEN f.start_ts>p.pause_ts AND f.start_ts<LEAST(p.pause_ts+INTERVAL '72 hours',s.outtime)
                 THEN f.start_ts END) AS restart_gt,
        MIN(CASE WHEN f.start_ts>=p.pause_ts AND f.start_ts<LEAST(s.first_feed_ts+INTERVAL '72 hours',s.outtime)
                 THEN f.start_ts END) AS legacy_restart_ge,
        MIN(CASE WHEN f.start_ts>p.pause_ts AND f.start_ts<LEAST(s.first_feed_ts+INTERVAL '72 hours',s.outtime)
                 THEN f.start_ts END) AS legacy_restart_gt
      FROM feeding_stays s JOIN first_pause p USING(stay_id)
      LEFT JOIN valid_feeds f ON f.stay_id=s.stay_id AND f.start_ts>=p.pause_ts-INTERVAL '72 hours'
      GROUP BY ALL
    """)


def make_output_dirs(out):
    out.mkdir(parents=True,exist_ok=True)
    return out


def run_cohort_summaries(con,out):
    # Cohort flow
    flow=[]
    def scalar(q): return int(con.execute(q).fetchone()[0] or 0)
    flow.append(('adult_ICU_stays',scalar('SELECT COUNT(*) FROM adult_icu')))
    flow.append(('adult_ICU_stays_with_any_curated_resection_code',scalar('SELECT COUNT(DISTINCT a.stay_id) FROM adult_icu a JOIN proc_events p USING(hadm_id) JOIN code_config c ON c.icd_code=p.icd_code AND c.icd_version=p.icd_version')))
    flow.append(('stays_with_curated_resection_code_0_to_7_days_prior',scalar('SELECT COUNT(*) FROM cohort')))
    flow.append(('primary_postoperative_stays',scalar("SELECT COUNT(*) FROM cohort WHERE in_primary=1")))
    flow.append(('primary_stays_with_strict_continuous_EN',scalar("SELECT COUNT(*) FROM feeding_stays WHERE in_primary=1")))
    flow.append(('primary_EN_stays_with_first_Paused_record',scalar("SELECT COUNT(*) FROM pause_followup WHERE in_primary=1")))
    flow.append(('primary_pause_stays_without_overlap_ambiguity',scalar("SELECT COUNT(*) FROM pause_followup WHERE in_primary=1 AND overlap_ambiguous=0")))
    write_csv(out/'00_cohort_flow.csv',['stage','n_stays'],flow)

    # Cohort and feed counts by stratum
    groups=['primary','sensitivity_only']
    rows=[]
    for g in groups:
        n=scalar(f"SELECT COUNT(*) FROM cohort WHERE analysis_tier='{g}'")
        en=scalar(f"SELECT COUNT(*) FROM feeding_stays WHERE analysis_tier='{g}'")
        p=scalar(f"SELECT COUNT(*) FROM pause_followup WHERE analysis_tier='{g}'")
        overlap=scalar(f"SELECT COUNT(*) FROM pause_followup WHERE analysis_tier='{g}' AND overlap_ambiguous=1")
        rows.append((g,n,scalar(f"SELECT COUNT(DISTINCT subject_id) FROM cohort WHERE analysis_tier='{g}'"),en,p,overlap))
    write_csv(out/'01_cohort_and_feeding_counts.csv',['group','stays','unique_patients','strict_EN_stays','paused_stays','overlap_ambiguous_pauses'],rows)

    # Explicitly distinguish ICU stays, hospital admissions, and unique patients.
    unit_rows = []
    for label, table, where in [
        ('primary', 'cohort', 'in_primary=1'),
        ('primary_EN', 'feeding_stays', 'in_primary=1'),
        ('primary_paused', 'pause_followup', 'in_primary=1')
    ]:
        q = con.execute(f"""
          WITH base AS (
            SELECT stay_id, subject_id, hadm_id FROM {table} WHERE {where}
          ), by_hadm AS (
            SELECT hadm_id, COUNT(*) AS n_stays FROM base GROUP BY hadm_id
          ), by_subject AS (
            SELECT subject_id, COUNT(*) AS n_stays FROM base GROUP BY subject_id
          )
          SELECT
            (SELECT COUNT(*) FROM base) AS icu_stays,
            (SELECT COUNT(DISTINCT subject_id) FROM base) AS unique_patients,
            (SELECT COUNT(DISTINCT hadm_id) FROM base) AS unique_admissions,
            (SELECT COUNT(*) FROM by_hadm WHERE n_stays>1) AS admissions_with_multiple_icu_stays,
            (SELECT COUNT(*) FROM by_subject WHERE n_stays>1) AS patients_with_multiple_icu_stays
        """).fetchone()
        unit_rows.append((label, *[nsmall(x) for x in q]))
    write_csv(out/'20_cohort_unit_audit.csv',
      ['cohort','icu_stays','unique_patients','unique_admissions',
       'admissions_with_multiple_icu_stays','patients_with_multiple_icu_stays'], unit_rows)

    # Procedure code audit
    code_rows=con.execute("""
      SELECT c.icd_version,c.icd_code,d.long_title,
        COUNT(DISTINCT p.hadm_id) admissions_n,
        COUNT(DISTINCT s.stay_id) ICU_stays_n,
        CASE
          WHEN c.icd_version=9 AND c.icd_code='4382' THEN 'AUTHOR_REVIEW_BARIATRIC_SLEEVE'
          WHEN c.icd_version=9 AND TRY_CAST(c.icd_code AS INTEGER) BETWEEN 4840 AND 4849 THEN 'AUTHOR_REVIEW_PULL_THROUGH'
          WHEN c.icd_version=10 AND c.icd_code IN ('0DTH8ZZ','0DTN7ZZ','0DTNFZZ','0DTP7ZZ') THEN 'AUTHOR_REVIEW_NATURAL_ORIFICE_APPROACH'
          WHEN lower(coalesce(d.long_title,'')) LIKE '%other%' OR lower(coalesce(d.long_title,'')) LIKE '%unspecified%' OR lower(coalesce(d.long_title,'')) LIKE '%not otherwise specified%' OR lower(coalesce(d.long_title,'')) LIKE '%partial excision%' THEN 'AUTHOR_REVIEW'
          ELSE 'No obvious text flag'
        END AS description_review_flag
      FROM code_config c
      LEFT JOIN v_d_icd_procedures d ON upper(trim(d.icd_code))=c.icd_code AND TRY_CAST(d.icd_version AS INTEGER)=c.icd_version
      LEFT JOIN proc_events p ON p.icd_code=c.icd_code AND p.icd_version=c.icd_version
      LEFT JOIN adult_icu s ON s.hadm_id=p.hadm_id
      GROUP BY 1,2,3,6 ORDER BY 1,2
    """).fetchall()
    write_csv(out/'02_procedure_code_audit.csv',['icd_version','icd_code','dictionary_description','admissions_with_code','adult_ICU_stays_with_code','description_review_flag'],
      [(a,b,c,nsmall(d),nsmall(e),f) for a,b,c,d,e,f in code_rows])

    # Feeding item audit and status/unit quality
    item_rows=con.execute("""
      SELECT d.itemid,d.candidate_label,d.dictionary_label,d.category,d.linksto,d.dictionary_match,
        COUNT(q.itemid) raw_rows,
        SUM(CASE WHEN q.valid_strict=1 THEN 1 ELSE 0 END) strict_valid_rows,
        COUNT(DISTINCT CASE WHEN q.valid_strict=1 THEN q.stay_id END) strict_valid_stays,
        SUM(CASE WHEN q.valid_strict=1 AND q.in_primary=1 THEN 1 ELSE 0 END) primary_valid_rows,
        COUNT(DISTINCT CASE WHEN q.valid_strict=1 AND q.in_primary=1 THEN q.stay_id END) primary_valid_stays
      FROM formula_dictionary d LEFT JOIN (SELECT q.*, c.in_primary FROM feed_quality q JOIN cohort c USING(stay_id)) q ON q.itemid=d.itemid
      GROUP BY 1,2,3,4,5,6 ORDER BY 1
    """).fetchall()
    write_csv(out/'03_formula_item_audit.csv',['itemid','candidate_label','d_items_label','category','linksto','dictionary_match','raw_rows','strict_valid_rows','strict_valid_stays','primary_valid_rows','primary_valid_stays'],
      [(a,b,c,d,e,f,nsmall(g),nsmall(h),nsmall(i),nsmall(j),nsmall(k)) for a,b,c,d,e,f,g,h,i,j,k in item_rows])
    unit_rows=con.execute("SELECT rateuom,amountuom,COUNT(*) n FROM possible_feeds GROUP BY 1,2 ORDER BY n DESC").fetchall()
    write_csv(out/'04_input_unit_audit.csv',['rateuom','amountuom','raw_candidate_rows'],[(a,b,nsmall(c)) for a,b,c in unit_rows])
    status_rows=con.execute("SELECT status,COUNT(*) n FROM possible_feeds GROUP BY 1 ORDER BY n DESC").fetchall()
    write_csv(out/'05_input_status_audit.csv',['statusdescription_normalized','raw_candidate_rows'],[(a,nsmall(b)) for a,b in status_rows])

    # Preserve local status labels and their order-category context for manual review.
    status_context = con.execute("""
      SELECT status, order_category_description,
             COUNT(*) AS raw_rows,
             COUNT(DISTINCT itemid) AS distinct_itemids
      FROM possible_feeds
      GROUP BY status, order_category_description
      ORDER BY status, order_category_description
    """).fetchall()
    documented_statuses = {'changed','paused','finishedrunning','stopped','flushed','rewritten'}
    status_review_rows = []
    for status, category_desc, n, n_items in status_context:
        if status in documented_statuses:
            flag = 'DOCUMENTED_STATUS'
        elif status == 'changedose/rate':
            flag = 'UNVERIFIED_ALIAS_CANDIDATE_FOR_CHANGED'
        else:
            flag = 'UNDOCUMENTED_STATUS_REVIEW'
        status_review_rows.append((status, category_desc or 'Missing', nsmall(n), nsmall(n_items), flag))
    write_csv(out/'22_status_dictionary_review.csv',
      ['observed_statusdescription','ordercategorydescription','raw_candidate_rows',
       'distinct_itemids','documentation_review_flag'], status_review_rows)

    link=con.execute("""
      SELECT COUNT(*) total_valid_rows,
        SUM(CASE WHEN linkorderid IS NULL THEN 1 ELSE 0 END) missing_linkorderid,
        SUM(CASE WHEN orderid IS NULL THEN 1 ELSE 0 END) missing_orderid,
        SUM(CASE WHEN linkorderid IS NULL AND orderid IS NULL THEN 1 ELSE 0 END) missing_both,
        COUNT(DISTINCT linkorderid) distinct_linkorderids,
        COUNT(DISTINCT orderid) distinct_orderids
      FROM valid_feeds
    """).fetchone()
    write_csv(out/'18_order_linkage_audit.csv',['valid_EN_rows','rows_missing_linkorderid','rows_missing_orderid','rows_missing_both_order_keys','distinct_linkorderids','distinct_orderids'],[tuple(nsmall(x) for x in link)])
    consistency=con.execute("""
      SELECT COUNT(*) n,
        median(amount/(rate*(date_diff('second',start_ts,end_ts)/3600.0))) median_amount_rate_ratio,
        quantile_cont(amount/(rate*(date_diff('second',start_ts,end_ts)/3600.0)),0.25) q25_ratio,
        quantile_cont(amount/(rate*(date_diff('second',start_ts,end_ts)/3600.0)),0.75) q75_ratio,
        100.0*SUM(CASE WHEN abs(amount/(rate*(date_diff('second',start_ts,end_ts)/3600.0))-1.0)<=0.10 THEN 1 ELSE 0 END)/COUNT(*) within_10pct
      FROM valid_feeds WHERE rate>0 AND end_ts>start_ts
    """).fetchone()
    write_csv(out/'19_amount_rate_consistency.csv',['valid_rows','median_amount_to_rate_duration_ratio','q25_ratio','q75_ratio','rows_within_10pct_of_rate_times_duration_pct'],[consistency])
    timing=con.execute("SELECT surgical_day_delta,COUNT(DISTINCT stay_id) n FROM surgical_candidates_raw GROUP BY 1 ORDER BY 1").fetchall()
    write_csv(out/'06_procedure_day_difference.csv',['days_from_procedure_code_date_to_ICU_entry','ICU_stays'],[(a,nsmall(b)) for a,b in timing])

    # Baseline characteristics in long form: stay-level descriptives only.
    baseline=[]
    cohorts=[('all_primary','cohort','in_primary=1'),('EN_primary','feeding_stays','in_primary=1'),('paused_primary','pause_followup','in_primary=1')]
    for label,table,where in cohorts:
        rows0=con.execute(f'SELECT * FROM {table} WHERE {where}').fetchall()
        cols=[x[0] for x in con.description]
        # rows are small; construct dicts
        data=[dict(zip(cols,r)) for r in rows0]
        n=len(data)
        if not n: continue
        ages=[r['age'] for r in data if r.get('age') is not None]
        med,q1,q3=median_iqr(ages)
        baseline.append((label,'Age, years','median (IQR)',f'{med:.1f} ({q1:.1f}-{q3:.1f})',len(ages)))
        for var,title in [('gender','Recorded sex'),('first_careunit','First ICU care unit'),('admission_type','Admission type'),('procedure_group','Procedure group'),('timing_stratum','Procedure-to-ICU timing stratum')]:
            counts=Counter(str(r.get(var) or 'Missing') for r in data)
            for cat,count in sorted(counts.items(),key=lambda x:(-x[1],x[0])):
                baseline.append((label,title,cat,(f'{count}/{n} ({100*count/n:.1f}%)' if count>=MIN_CELL or count==0 else '<11; percentage suppressed'),n))
        icu_los=[(r['outtime']-r['intime']).total_seconds()/86400 for r in data if r.get('outtime') and r.get('intime')]
        hosp_los=[(r['dischtime']-r['admittime']).total_seconds()/86400 for r in data if r.get('dischtime') and r.get('admittime')]
        for title,vals in [('ICU length of stay, days',icu_los),('Hospital length of stay, days',hosp_los)]:
            m,a,b=median_iqr(vals)
            baseline.append((label,title,'median (IQR)',f'{m:.2f} ({a:.2f}-{b:.2f})',len(vals)))
        n_admissions = len({r['hadm_id'] for r in data if r.get('hadm_id') is not None})
        death_admissions = len({r['hadm_id'] for r in data
                                if r.get('hadm_id') is not None and r.get('hospital_expire_flag') == 1})
        baseline.append((label,'Hospital admissions ending in death','n/N (%)',
                         (f'{death_admissions}/{n_admissions} ({100*death_admissions/n_admissions:.1f}%)'
                          if n_admissions and (death_admissions == 0 or death_admissions >= MIN_CELL)
                          else ('<11; percentage suppressed' if n_admissions else 'NA')),
                         n_admissions))
    write_csv(out/'07_baseline_characteristics.csv',['cohort','variable','statistic','value','denominator'],baseline)

    # Time from ICU entry to first eligible continuous EN record (primary cohort).
    first_feed_timing = con.execute("""
      SELECT COUNT(*) AS n,
             median(date_diff('second',c.intime,f.first_feed_ts)/3600.0) AS median_h,
             quantile_cont(date_diff('second',c.intime,f.first_feed_ts)/3600.0,0.25) AS q25_h,
             quantile_cont(date_diff('second',c.intime,f.first_feed_ts)/3600.0,0.75) AS q75_h
      FROM cohort c JOIN first_feed f USING(stay_id)
      WHERE c.in_primary=1 AND f.first_feed_ts>=c.intime
    """).fetchone()
    write_csv(out/'23_first_feed_timing_summary.csv',
      ['primary_EN_stays','median_hours_from_ICU_entry','q25_hours','q75_hours'],
      [(first_feed_timing[0],
        f'{first_feed_timing[1]:.2f}' if first_feed_timing[1] is not None else 'NA',
        f'{first_feed_timing[2]:.2f}' if first_feed_timing[2] is not None else 'NA',
        f'{first_feed_timing[3]:.2f}' if first_feed_timing[3] is not None else 'NA')])


def make_event_records(con, cohort='primary', same_time=False, legacy=False, max_days=None):
    rows=con.execute('SELECT subject_id,analysis_tier,timing_stratum,pause_ts,pause_followup_end,legacy_followup_end,restart_ge,restart_gt,legacy_restart_ge,legacy_restart_gt,overlap_ambiguous,outtime,first_feed_ts,intime FROM pause_followup').fetchall()
    records=[]
    for sid,tier,stratum,pause,end,legacy_end,rge,rgt,lrge,lrgt,overlap,out,first_feed_ts,icu_in in rows:
        if cohort=='primary' and tier!='primary': continue
        if cohort=='all' and tier not in ('primary','sensitivity_only'): continue
        if cohort=='days_0_1' and stratum not in ('same_day_PACU_proxy','previous_calendar_day'): continue
        if cohort=='days_4_7' and stratum!='4_to_7_days_prior': continue
        if overlap==1: continue
        if max_days is not None:
            # only strata where the procedure code date was 0-1 days before ICU entry
            if stratum not in ('same_day_PACU_proxy','previous_calendar_day'): continue
        follow_end=legacy_end if legacy else end
        same_timestamp_pair = (lrge == pause) if legacy else (rge == pause)
        if same_timestamp_pair and not same_time:
            continue
        resume=lrge if (same_time and legacy) else (lrgt if legacy else (rge if same_time else rgt))
        target_end=(first_feed_ts+dt.timedelta(hours=72)) if legacy else (pause+dt.timedelta(hours=72))
        if out is not None: target_end=min(target_end,out)
        if follow_end is None or pause is None or follow_end<=pause: continue
        if resume is not None and resume>=follow_end: resume=None
        if resume is not None:
            t=(resume-pause).total_seconds()/3600
            kind='restart'
        else:
            t=(follow_end-pause).total_seconds()/3600
            if out is not None and out < ((first_feed_ts+dt.timedelta(hours=72)) if legacy else (pause+dt.timedelta(hours=72))):
                kind='icu_exit'
            else:
                kind='censor'
        if not same_time and t<=0: continue
        records.append({'subject_id':int(sid),'time':max(0.0,float(t)),'kind':kind,'timing_stratum':stratum})
    return records


def run_time_to_restart(con,out):
    scenarios=[('PRIMARY','primary',False,False),('PRIMARY_same_timestamp','primary',True,False),
               ('DAYS_0_1','days_0_1',False,False),('DAYS_0_1_same_timestamp','days_0_1',True,False),
               ('PRIMARY_legacy_first_feed_window','primary',False,True),
               ('PRIMARY_legacy_same_timestamp','primary',True,True)]
    rows=[]; follow=[]
    horizons=[6,12,24,36,48]
    for idx,(name,cohort,same,legacy) in enumerate(scenarios):
        rec=make_event_records(con,cohort=cohort,same_time=same,legacy=legacy)
        n=len(rec); unique=len(set(r['subject_id'] for r in rec))
        for h in horizons:
            est=cif_aj([(r['time'],r['kind']) for r in rec],h)
            lo,hi=bootstrap_cluster(rec,h,seed=RNG_SEED+h,metric='cif') if n>=2 else (None,None)
            rows.append((name,h,n,unique,round(est*100,1) if math.isfinite(est) else 'NA',round(lo*100,1) if lo is not None else 'NA',round(hi*100,1) if hi is not None else 'NA',at_risk_n([(r['time'],r['kind']) for r in rec],h)))
        pos=[r['time'] for r in rec if r['kind']=='restart' and r['time']>0]
        q2,q1,q3=median_iqr(pos)
        events=sum(r['kind']=='restart' for r in rec)
        comp=sum(r['kind']=='icu_exit' for r in rec)
        cens=sum(r['kind']=='censor' for r in rec)
        rmst=rmst_without_restart([(r['time'],r['kind']) for r in rec],24) if rec else None
        lo_r,hi_r=bootstrap_cluster(rec,24,seed=RNG_SEED+99,metric='rmst') if n>=2 else (None,None)
        follow.append((name,n,unique,events,comp,cens,
                       f'{q2:.2f}' if q2 is not None else 'NA',
                       f'{q1:.2f}' if q1 is not None else 'NA',
                       f'{q3:.2f}' if q3 is not None else 'NA',
                       f'{rmst:.2f}' if rmst is not None else 'NA',
                       f'{lo_r:.2f}' if lo_r is not None else 'NA',
                       f'{hi_r:.2f}' if hi_r is not None else 'NA'))
    write_csv(out/'08_resumption_CIF_bootstrap.csv',
      ['scenario','hours_after_pause','analysis_n','unique_patients','CIF_restart_pct','bootstrap_95CI_lower_pct','bootstrap_95CI_upper_pct','at_risk_icu_stays_n'],rows)
    write_csv(out/'09_resumption_followup_summary.csv',
      ['scenario','analysis_n','unique_patients','restart_events','ICU_exit_competing_events','right_censored','positive_restart_median_h','positive_restart_q25_h','positive_restart_q75_h','24h_mean_time_without_recorded_restart_h','bootstrap_95CI_lower_h','bootstrap_95CI_upper_h'],follow)

    # Transparent accounting of overlap and same-timestamp exclusions.
    primary_pause_n, overlap_n, same_timestamp_n = con.execute("""
      SELECT COUNT(*),
             SUM(CASE WHEN overlap_ambiguous=1 THEN 1 ELSE 0 END),
             SUM(CASE WHEN overlap_ambiguous=0 AND restart_ge=pause_ts THEN 1 ELSE 0 END)
      FROM pause_followup WHERE analysis_tier='primary'
    """).fetchone()
    primary_records = make_event_records(con, cohort='primary', same_time=False, legacy=False)
    same_time_records = make_event_records(con, cohort='primary', same_time=True, legacy=False)
    event_audit = []
    for label, records in [('PRIMARY', primary_records), ('PRIMARY_same_timestamp', same_time_records)]:
        event_audit.append((
            label, primary_pause_n, overlap_n or 0, same_timestamp_n or 0, len(records),
            sum(r['kind']=='restart' for r in records),
            sum(r['kind']=='icu_exit' for r in records),
            sum(r['kind']=='censor' for r in records),
            sum(r['time']==0 and r['kind']=='restart' for r in records)
        ))
    write_csv(out/'21_event_timing_eligibility_audit.csv',
      ['scenario','primary_pause_cohort_n','overlap_ambiguous_excluded_n',
       'same_timestamp_pairs_n','analysis_n','restart_events',
       'competing_icu_exits','right_censored','zero_time_restarts'], event_audit)


def integrate_volume_by_order(intervals, window_start, window_end):
    """Integrate mL/h. De-duplicate overlapping rows within linkorderid/orderid; sum independent orders."""
    groups=defaultdict(list)
    for start,end,rate,order_key in intervals:
        if start is None or end is None or rate is None or rate<=0: continue
        lo=max(start,window_start); hi=min(end,window_end)
        if lo>=hi: continue
        groups[str(order_key)].append((lo,hi,float(rate)))
    all_times=set([window_start,window_end])
    for segs in groups.values():
        for lo,hi,rate in segs: all_times.add(lo); all_times.add(hi)
    ts=sorted(all_times); total=0.0
    for a,b in zip(ts[:-1],ts[1:]):
        if b<=a: continue
        mid=(a+(b-a)/2)
        sum_rate=0.0
        for segs in groups.values():
            active=[r for lo,hi,r in segs if lo<=mid<hi]
            if active: sum_rate += max(active)
        total += sum_rate*(b-a).total_seconds()/3600
    return total


def run_volume(con,out):
    obs=con.execute('SELECT stay_id,analysis_tier,first_feed_ts,outtime FROM feeding_stays').fetchall()
    ev=con.execute('SELECT stay_id,start_ts,end_ts,rate,COALESCE(linkorderid,orderid) order_key FROM feeding_events').fetchall()
    by=defaultdict(list)
    for stay,start,end,rate,key in ev: by[int(stay)].append((start,end,float(rate),key if key is not None else f'{stay}:{start}'))
    rows=[]
    for group in ('primary','sensitivity_only'):
        for hours in (24,48):
            vals=[]
            for stay,tier,t0,outtime in obs:
                if tier!=group or outtime is None or t0 is None or outtime<t0+dt.timedelta(hours=hours): continue
                vals.append(integrate_volume_by_order(by.get(int(stay),[]),t0,t0+dt.timedelta(hours=hours)))
            med,q1,q3=median_iqr(vals)
            rows.append((group,hours,len(vals),f'{med:.1f}' if med is not None else 'NA',f'{q1:.1f}' if q1 is not None else 'NA',f'{q3:.1f}' if q3 is not None else 'NA'))
    write_csv(out/'10_recorded_formula_volume.csv',['cohort','complete_window_hours_from_first_EN','n_complete_stays','median_mL','q25_mL','q75_mL'],rows)


def run_icu_entry_states(con,out):
    obs=con.execute('SELECT stay_id,intime,outtime,first_feed_ts FROM cohort LEFT JOIN first_feed USING(stay_id) WHERE in_primary=1').fetchall()
    rows=[]
    for hours in (24,48):
        states=Counter()
        for stay,intime,outtime,first in obs:
            cut=intime+dt.timedelta(hours=hours)
            if first is not None and first<=cut:
                states['eligible_continuous_EN_recorded']=states['eligible_continuous_EN_recorded']+1
            elif outtime is not None and outtime<=cut:
                states['ICU_exit_before_eligible_EN_record']=states['ICU_exit_before_eligible_EN_record']+1
            else:
                states['still_in_ICU_without_eligible_EN_record']=states['still_in_ICU_without_eligible_EN_record']+1
        n=len(obs); en=states['eligible_continuous_EN_recorded']
        # Wilson interval
        p=en/n
        ci_low, ci_high = wilson_interval_pct(en, n)
        rows.append((hours,n,en,f'{100*p:.1f}',f'{ci_low:.1f}',f'{ci_high:.1f}',
                     states['ICU_exit_before_eligible_EN_record'],states['still_in_ICU_without_eligible_EN_record']))
    write_csv(out/'11_ICU_entry_recording_states.csv',['hours_after_ICU_entry','N_stays','EN_recorded_n','EN_recorded_pct','Wilson_95CI_low_pct','Wilson_95CI_high_pct','ICU_exit_without_EN_record_n','still_in_ICU_without_EN_record_n'],rows)


def run_outcomes(con,out):
    rows=[]
    for label,table,where in [('primary','cohort','in_primary=1'),('EN_primary','feeding_stays','in_primary=1'),('paused_primary','pause_followup','in_primary=1')]:
        r=con.execute(f"""SELECT COUNT(*) n,
          COUNT(DISTINCT subject_id) patients,
          COUNT(DISTINCT hadm_id) admissions,
          median(date_diff('minute',intime,outtime)/1440.0) median_icu_los,
          median(date_diff('minute',admittime,dischtime)/1440.0) median_hosp_los,
          COUNT(DISTINCT CASE WHEN hospital_expire_flag=1 THEN hadm_id END) admissions_with_hospital_death,
          COUNT(DISTINCT CASE WHEN hospital_expire_flag=1 THEN subject_id END) patients_with_hospital_death
          FROM {table} WHERE {where}""").fetchone()
        rows.append((label,r[0],r[1],r[2],r[3],r[4],nsmall(r[5]),nsmall(r[6])))
    write_csv(out/'12_descriptive_clinical_outcomes.csv',
      ['cohort','icu_stays','unique_patients','unique_admissions','median_ICU_LOS_days',
       'median_hospital_LOS_days','unique_admissions_ending_in_hospital_death',
       'unique_patients_with_hospital_death'],rows)


def run_nursing(con,root,out):
    # Read only selected fields from the two event tables; charttime is the observation-time proxy.
    for name in ['chartevents','outputevents']:
        p=root/FILES[name]
        if not p.exists(): raise FileNotFoundError(str(p))
        ps=p.as_posix().replace("'","''")
        # MIMIC-IV CSV column order is not stay_id/charttime/itemid.  COLUMNS={...}
        # would map positionally and silently misread bedside observations.
        # Read actual header names, then project only the three required fields.
        con.execute(f"""CREATE OR REPLACE VIEW v_{name} AS
          SELECT TRY_CAST(stay_id AS BIGINT) AS stay_id,
                 TRY_CAST(charttime AS TIMESTAMP) AS charttime,
                 TRY_CAST(itemid AS BIGINT) AS itemid
          FROM read_csv_auto('{ps}', HEADER=TRUE, ALL_VARCHAR=TRUE, COMPRESSION='gzip')""")
    con.execute("""
      CREATE TEMP TABLE nursing_events AS
      SELECT TRY_CAST(c.stay_id AS BIGINT) stay_id, TRY_CAST(c.charttime AS TIMESTAMP) chart_ts,
             n.domain, 'chartevents' source_table
      FROM v_chartevents c JOIN nursing_config n ON TRY_CAST(c.itemid AS BIGINT)=n.itemid
      JOIN pause_followup p ON p.stay_id=TRY_CAST(c.stay_id AS BIGINT)
      WHERE p.analysis_tier='primary' AND n.source_table='chartevents'
        AND TRY_CAST(c.charttime AS TIMESTAMP)>=p.intime
        AND TRY_CAST(c.charttime AS TIMESTAMP)<p.outtime
      UNION ALL
      SELECT TRY_CAST(o.stay_id AS BIGINT), TRY_CAST(o.charttime AS TIMESTAMP),
             n.domain, 'outputevents'
      FROM v_outputevents o JOIN nursing_config n ON TRY_CAST(o.itemid AS BIGINT)=n.itemid
      JOIN pause_followup p ON p.stay_id=TRY_CAST(o.stay_id AS BIGINT)
      WHERE p.analysis_tier='primary' AND n.source_table='outputevents'
        AND TRY_CAST(o.charttime AS TIMESTAMP)>=p.intime
        AND TRY_CAST(o.charttime AS TIMESTAMP)<p.outtime
    """)
    data=con.execute("""
      SELECT stay_id,subject_id,intime,outtime,pause_ts,
        CASE WHEN overlap_ambiguous=0 AND restart_ge IS NOT NULL
                    AND restart_ge>pause_ts AND restart_gt<pause_followup_end
             THEN restart_gt ELSE NULL END AS positive_time_restart
      FROM pause_followup WHERE analysis_tier='primary'
    """).fetchall()
    stays={int(s):(sub,ii,oo,pp,rr) for s,sub,ii,oo,pp,rr in data}
    ev=con.execute('SELECT stay_id,domain,chart_ts FROM nursing_events').fetchall()
    by=defaultdict(list)
    for s,dom,t in ev: by[int(s)].append((dom,t))
    domains=con.execute('SELECT DISTINCT domain FROM nursing_config ORDER BY domain').fetchall()
    windows=['pre_pause_24h','at_pause_timestamp','post_pause_24h','pause_to_restart','post_restart_6h']
    counts=defaultdict(set); den=defaultdict(set)
    for sid,(sub,intime,outtime,pause,restart) in stays.items():
        if pause is None: continue
        if intime<=pause-dt.timedelta(hours=24): den['pre_pause_24h'].add(sid)
        den['at_pause_timestamp'].add(sid)
        if outtime>=pause+dt.timedelta(hours=24): den['post_pause_24h'].add(sid)
        if restart is not None and restart>pause: den['pause_to_restart'].add(sid)
        if restart is not None and outtime>=restart+dt.timedelta(hours=6): den['post_restart_6h'].add(sid)
        for dom,t in by.get(sid,[]):
            if pause-dt.timedelta(hours=24)<=t<pause and sid in den['pre_pause_24h']: counts[(dom,'pre_pause_24h')].add(sid)
            if t==pause: counts[(dom,'at_pause_timestamp')].add(sid)
            if pause<t<pause+dt.timedelta(hours=24) and sid in den['post_pause_24h']: counts[(dom,'post_pause_24h')].add(sid)
            if restart is not None and pause<t<restart: counts[(dom,'pause_to_restart')].add(sid)
            if restart is not None and restart<=t<restart+dt.timedelta(hours=6) and sid in den['post_restart_6h']: counts[(dom,'post_restart_6h')].add(sid)
    rows=[]
    domain_names=[r[0] for r in con.execute('SELECT DISTINCT domain FROM nursing_config ORDER BY domain').fetchall()]
    for dom in domain_names:
      for w in windows:
        ss=counts.get((dom,w),set())
        d=len(den[w]); n=len(ss)
        rows.append((dom,w,nsmall(n),d,(f'{100*n/d:.1f}' if d and (n>=MIN_CELL or n==0) else ('SUPPRESSED' if d else 'NA'))))
    write_csv(out/'13_bedside_nursing_relevant_charting.csv',['domain','window','stays_with_any_record_suppressed_lt11','eligible_stays_complete_window','percent_of_eligible'],rows)
    den_rows=[(w,len(den[w])) for w in windows]
    write_csv(out/'14_bedside_charting_window_denominators.csv',['window','eligible_stays'],den_rows)


def write_validation_report(con,out,root,config_dir):
    lines=[]
    lines.append('MIMIC-IV v3.1 | final reproducibility and phenotype validation report')
    lines.append(f'Generated: {dt.datetime.now().isoformat(timespec="seconds")}')
    lines.append('Dataset root: [LOCAL PATH WITHHELD]')
    lines.append('Status: FINAL ANALYTIC CANDIDATE — AUTHOR CLINICAL/ETHICS REVIEW STILL REQUIRED')
    lines.append('')
    lines.append('Status mapping is provisional. Current MIMIC-IV documentation describes Changed, Paused, FinishedRunning, Stopped, and Flushed; legacy documentation also describes Rewritten. Local candidate-item rows include ChangedDose/Rate and Bolus labels. ChangedDose/Rate is provisionally retained as a possible Changed alias but requires source-version confirmation; Bolus and other unrecognized labels are excluded.')
    lines.append('The first Paused row is the interruption anchor. A later eligible feed start is the operational restart marker; it is not the clinical decision time.')
    lines.append('ICU exit before restart is treated as a competing observation-ending event; this may represent transfer/discharge or death and is not synonymous with live ICU discharge.')
    lines.append('Primary follow-up is 72 hours from the pause or earlier ICU exit. A legacy sensitivity window is 72 hours from first eligible feeding or earlier ICU exit.')
    lines.append('Selected fields are structured bedside charting relevant to nutrition care. MIMIC-IV v3.1 caregiver records contain caregiver IDs but no professional-role label; records are not attributed specifically to nurses.')
    lines.append('')
    n_code=con.execute('SELECT COUNT(*) FROM code_config').fetchone()[0]
    n_item=con.execute('SELECT COUNT(*) FROM formula_config').fetchone()[0]
    invalid_items=con.execute("SELECT COUNT(*) FROM formula_dictionary WHERE dictionary_match=0").fetchone()[0]
    unknown_status=con.execute("SELECT COUNT(*) FROM possible_feeds WHERE status NOT IN ('changed','changedose/rate','finishedrunning','paused','stopped','flushed')").fetchone()[0]
    lines.append(f'Curated procedure codes: {n_code}; candidate formula item IDs: {n_item}.')
    lines.append(f'Formula item IDs not matching d_items category/linksto expectation: {invalid_items}.')
    lines.append(f'Candidate input rows with an unrecognized statusdescription: {nsmall(unknown_status)}.')
    observed_statuses = [r[0] for r in con.execute('SELECT DISTINCT status FROM possible_feeds ORDER BY status').fetchall()]
    lines.append('Observed candidate-item statusdescription labels: ' + ', '.join(observed_statuses))
    lines.append('')
    lines.append('Required author actions before submission:')
    lines.append('1. Clinically review the 82 procedure codes and confirm the intended anatomical/procedure scope; specific flags identify bariatric sleeve gastrectomy, pull-through procedures, and natural-orifice endoscopic approaches for separate adjudication.')
    lines.append('2. Confirm the 36 candidate formula item IDs; retain only item IDs with valid administered records in the final phenotype table.')
    lines.append('3. Inspect the status/unit/item audit outputs and reconcile any deviations from the prespecified filters.')
    lines.append('4. Confirm the 72-h post-pause primary estimand and compare it with the legacy first-feed-anchored sensitivity analysis.')
    lines.append('5. Review the bedside charting domains; do not attribute entries to nurses unless a validated caregiver-role source exists.')
    lines.append('6. Complete local ethics determination, PhysioNet DUA compliance, author roles, funding, conflicts, and AI disclosure.')
    lines.append('7. Confirm the observed ChangedDose/Rate and Bolus labels against the exact local MIMIC-IV v3.1 file schema and codebook before finalizing the status mapping.')
    lines.append('8. Update manuscript tables/figures from this exact run; never carry forward old numbers if counts differ.')
    (out/'15_VALIDATION_REPORT.txt').write_text('\n'.join(lines),encoding='utf-8')


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--root',default=r'D:\MIMICIV')
    ap.add_argument('--out',default=None)
    ap.add_argument('--config',default=None)
    ap.add_argument('--nursing',action='store_true')
    ap.add_argument('--threads',type=int,default=4)
    ap.add_argument('--memory-gb',type=int,default=8)
    args=ap.parse_args()
    root=Path(args.root); out=Path(args.out) if args.out else root.parent/'MIMIC_ENTERAL_V2_RESULTS'
    config=Path(args.config) if args.config else Path(__file__).parent/'config'
    missing=[str(root/p) for k,p in FILES.items() if k not in ('chartevents','outputevents') and not (root/p).exists()]
    if args.nursing:
        missing += [str(root/FILES[k]) for k in ('chartevents','outputevents') if not (root/FILES[k]).exists()]
    if missing: raise SystemExit('Missing required files:\n'+'\n'.join(missing))
    out.mkdir(parents=True,exist_ok=True)
    con=duckdb.connect(':memory:')
    con.execute(f'SET threads={max(1,args.threads)}')
    con.execute(f"SET memory_limit='{max(4,args.memory_gb)}GB'")
    tmp=out/'_PRIVATE_TEMP_DELETE_AFTER_RUN'; tmp.mkdir(exist_ok=True)
    con.execute("SET temp_directory='"+tmp.as_posix().replace("'","''")+"'")
    try:
        print('Loading MIMIC-IV v3.1 source tables...')
        create_base_views(con,root,config,args)
        print('Summarizing cohort, feeding phenotype, code and unit audits...')
        run_cohort_summaries(con,out)
        print('Estimating time to restart with 72-h post-pause follow-up and sensitivity windows...')
        run_time_to_restart(con,out)
        print('Summarizing formula volumes and ICU-entry recording states...')
        run_volume(con,out); run_icu_entry_states(con,out); run_outcomes(con,out)
        if args.nursing:
            print('Scanning chartevents/outputevents for nursing-relevant structured bedside charting...')
            run_nursing(con,root,out)
        write_validation_report(con,out,root,config)
        manifest={
          'dataset':'MIMIC-IV v3.1','analysis_version':ANALYSIS_VERSION,'generated':dt.datetime.now().isoformat(timespec='seconds'),
          'software_environment': {'python':sys.version.split()[0], 'duckdb':getattr(duckdb,'__version__','unreported')},
          'procedure_codes':82,'candidate_formula_item_ids':36,'nursing_item_dictionary_rows':15,
          'primary_pause_followup':'72 hours from first documented Paused end-time or earlier ICU exit',
          'legacy_sensitivity':'72 hours from first eligible EN administration or earlier ICU exit',
          'bootstrap_replicates':BOOTSTRAP_REPS,'bootstrap_seed':RNG_SEED,
          'raw_patient_rows_exported':False,'nursing_scan_run':bool(args.nursing),
          'outputs':sorted(p.name for p in out.iterdir() if p.is_file()),
          'source_checksums_sha256': {q.name: hashlib.sha256(q.read_bytes()).hexdigest() for q in [Path(__file__), Path(__file__).with_name('verify_results_v3.py'), Path(__file__).with_name('phenotype_adjudication_v3.py'), Path(__file__).with_name('run_final.py'), Path(__file__).with_name('build_submission_documents.py'), Path(__file__).with_name('requirements.txt'), *(sorted(config.glob('*.csv')))] if q.exists()}
        }
        (out/'manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
        print('COMPLETE. Aggregate outputs:',out)
        print('Do not submit until clinical phenotype validation and author declarations are complete.')
    finally:
        con.close()
        try:
            if tmp.exists() and not any(tmp.iterdir()): tmp.rmdir()
        except OSError: pass

if __name__=='__main__': main()
