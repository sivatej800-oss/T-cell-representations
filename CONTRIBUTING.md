# Contributing

## Working locally

1. Create a branch from `main` using a short descriptive name.
2. Use a project-specific virtual environment and install `requirements.txt`.
3. Keep raw data and generated analysis artifacts in the ignored local directories.
4. Clear sensitive outputs from notebooks before sharing them.
5. Run `python scripts/check_repository.py` before opening a pull request.

## Notebook changes

- Keep configuration and input paths in an early, clearly labelled cell.
- Record dataset accessions, random seeds, package versions, and material parameters.
- Preserve sparse matrices throughout the full-data workflow.
- Add a short Markdown explanation when a result changes scientific interpretation.
- Place superseded exploratory work under `notebooks/experiments/` rather than mixing it into the canonical workflow.

## Pull requests

Describe the scientific or engineering question, the inputs used, the outputs affected, and the validation performed. Do not commit raw participant-level data, credentials, local absolute paths containing secrets, large binary artifacts, or model checkpoints.
