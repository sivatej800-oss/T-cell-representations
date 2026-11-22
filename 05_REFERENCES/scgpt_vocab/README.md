# scGPT gene vocabulary

`default_census_vocab.json` is the gene-token vocabulary from the upstream scGPT
repository (https://github.com/bowang-lab/scGPT, `scgpt/tokenizer/`). `LICENSE` is
that repository's MIT license. Stage 05 audits gene coverage against this file.

- sha256: `9fe73cf29929b5c66a5fffcfb2bcbdd3e38a3a086848cd12f8f918bcf458d75f`
- Tokens: 60,694.
- It is not byte-identical to the `vocab.json` in the model directory used for the
  GPU run, which prepends special tokens. It reproduces that run's 15,960 matched /
  3,729 unmatched gene split exactly, which is what makes it a valid stand-in.
