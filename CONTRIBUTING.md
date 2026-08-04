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
- Update `scripts/README.md` when a stage's sequential logic, inputs, outputs, or
  checkpoint behavior changes.

## Notebook and reference changes

- Treat everything under `archive/` as historical provenance rather than maintained
  pipeline code.
- Place additional superseded or diagnostic notebooks under `archive/notebooks/`.
- Store supplied external analysis code under `archive/references/`, not among
  executable scripts.
- Keep unused third-party source snapshots under `archive/vendor/`; active dependencies
  belong in `requirements.txt` or should be installed through R.
- Document any scientific interpretation changes in the root README or a focused reference note.

## Pull requests

Describe the scientific or engineering question, the inputs used, the outputs affected, and the validation performed. Do not commit raw participant-level data, credentials, local absolute paths containing secrets, large binary artifacts, or model checkpoints.
