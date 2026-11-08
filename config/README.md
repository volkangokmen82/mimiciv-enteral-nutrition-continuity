# Local configuration files

This repository deliberately does not distribute MIMIC-derived configuration files.

For the local one-shot runner, place the following files in a local configuration directory and pass it with `--config`:
- `procedure_codes.csv`
- `formula_item_ids.csv`
- `nursing_item_dictionary.csv`

Use the exact locally curated files already held under the authorized MIMIC-IV research workspace. Do not commit these files, MIMIC-IV source tables, patient-level extracts, or derived patient-level data to a public repository. The run manifest records SHA-256 hashes of the local configuration files so the analysis can be traced to the exact local versions.
