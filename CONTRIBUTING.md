# Contributing

## Working locally

1. Create a branch from `main` using a short descriptive name.
2. Use a project-specific virtual environment and install `requirements.txt`.
3. Keep raw data and generated analysis artifacts in the ignored local directories.
4. Keep each numbered Python stage self-contained; do not introduce hidden imports
   between pipeline files.
5. Run `python3 scripts/check_repository.py` before opening a pull request.

## Pipeline changes

- Record dataset accessions, random seeds, package versions, and material parameters.
- Preserve sparse matrices throughout the full-data workflow.
- Keep generated-data paths relative to the configured workspace; do not add machine-specific absolute paths.
- Add validations at artifact boundaries so interrupted stages can resume safely.
- Keep each stage's step comments accurate when its logic, inputs, outputs, or
  checkpoint behavior changes. Those comments are the per-stage documentation; there
  is no separate file describing what the scripts do.

## Notebook and reference changes

- Treat everything under `archive/` as historical provenance rather than maintained
  pipeline code.
- Place additional superseded or diagnostic notebooks under `archive/notebooks/`.
- Keep the paper authors' own code (`figure1d_author_workflow.R` and
  `author_preprocessing/`) separate from `scripts/`. It is reference material that
  documents what the published analysis did; it is never run by this pipeline.
- Keep unused third-party source snapshots under `archive/vendor/`; active dependencies
  belong in `requirements.txt` or should be installed through R.
- Document any scientific interpretation changes in the root README or a focused reference note.

## Pull requests

Describe the scientific or engineering question, the inputs used, the outputs affected, and the validation performed. Do not commit raw participant-level data, credentials, local absolute paths containing secrets, large binary artifacts, or model checkpoints.
