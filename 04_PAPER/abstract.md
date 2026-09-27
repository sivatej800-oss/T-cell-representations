# Abstract

*Draft. Bracketed text is a placeholder for the T-cell and monocyte subclustering,
which has not been run. Word count excluding placeholders: ~260.*

## Background

Single-cell foundation models such as scGPT are increasingly proposed as drop-in
replacements for principal component analysis (PCA) when representing cells, but
comparisons rarely hold the rest of the analysis fixed or separate geometric
cluster quality from biological correctness. Human cerebrospinal fluid (CSF) immune
cells, dominated by continuous T-cell and myeloid states, are a demanding test case.

## Methods

Using the public GSE200164 dataset (70,391 CSF immune cells from 59 individuals;
Piehl et al., 2022), the published Seurat clustering workflow was reproduced and
compared against a frozen, zero-shot scGPT embedding under identical downstream
clustering. Both representations were corrected for donor effects with Harmony.
Cluster structure was evaluated by donor NMI, local donor mixing and silhouette
width, with a 512-component PCA control to match scGPT's dimensionality. The two
frozen partitions were then assessed biologically using marker genes computed from
one shared RNA expression matrix, followed by manual annotation into the authors'
cell-type categories.

## Results

After Harmony, both representations carried equal donor information (NMI 0.051), but
PCA mixed donors more thoroughly at the local level (0.959 vs 0.867). scGPT produced
substantially better cluster separation: in the embedding space its partition scored
a mean silhouette of 0.207 against 0.103 for the PCA partition. The advantage
persisted at matched dimensionality (0.219 vs 0.160). PCA extended to 512 components
degraded rather than improved, losing its rare populations (9 clusters, smallest 378
cells) and scoring below zero in its own space. Despite these geometric differences,
both representations recovered the same five broad immune populations, with counts
differing by at most 83 cells. [T-cell and monocyte subclustering: RESULT PENDING.]

## Conclusions

A frozen scGPT embedding organizes CSF immune cells into more cleanly separated
clusters than PCA, and the advantage is not explained by dimensionality. At the level
of broad cell types, however, the two representations are biologically equivalent.
Better geometric separation is therefore not by itself evidence of better biology.
[Whether scGPT's advantage translates into more accurate subtype resolution: PENDING.]
