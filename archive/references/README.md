# Research references

This directory contains source material used to compare the maintained pipeline with the published analysis.

- `figure1d_author_workflow.R` is the authors' shared analysis code for integration, clustering, manual cell-type annotation, and downstream TCR work.

The reference script is not an executable project entry point: it contains placeholder
paths, depends on unpublished helper functions and metadata fields, and includes
downstream analyses that are outside the maintained pipeline. The core reproducible
parameters have been incorporated into
`../../scripts/03_reproduce_figure1d.R`.
