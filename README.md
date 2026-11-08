# MIMIC-IV enteral nutrition continuity analysis

Analysis code for a retrospective cohort study of recorded continuous enteral nutrition resumption after documented pauses in intensive care stays with coded abdominal resection.

## Requirements

- Python 3.10 or newer and DuckDB from `requirements.txt`.
- Authorized local access to MIMIC-IV v3.1 from [PhysioNet](https://physionet.org/content/mimiciv/3.1/). Dataset DOI: 10.13026/kpb9-mt58.
- Three study configuration CSVs in a local directory. See [configuration requirements](config/README.md).

## Run locally

```text
python -m pip install -r requirements.txt
python -m unittest discover -s tests -v
python run_final.py --root PATH_TO_MIMICIV --config PATH_TO_LOCAL_CONFIG --out PATH_TO_NEW_LOCAL_OUTPUT --threads 4 --memory-gb 8
```

The output directory must be empty or absent. The pipeline builds the cohort, runs the primary and phenotype sensitivity analyses, and validates local aggregate outputs. The primary pause-to-resumption analysis uses the Aalen–Johansen estimator, with intensive care exit as a competing event and administrative right censoring at 72 hours. Bootstrap intervals use patient-level clusters. The code records checksums of local inputs to support reproducibility.

This public distribution excludes manuscript generation and internal comparisons against fixed historical results. It retains the legacy follow-up definition as a computed sensitivity analysis. Exact replication requires the original study configuration. Clinical author review remains necessary.

Keep MIMIC-IV tables, credentials, local configuration and all generated results outside the public repository. The included tests use constructed inputs and require no patient data.
