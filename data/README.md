# Local data directory

Data files are intentionally excluded from Git.

The current workflow downloads public GSE200164 source files and creates sparse intermediates, AnnData objects, Seurat objects, and checkpoints locally. Follow the path configuration in the first data-preparation notebook and retain the original GEO accession and checksums in analysis logs.

Do not place credentials, controlled-access data, directly identifying data, or unapproved participant-level exports in this directory.
