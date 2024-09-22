# Bioinformatics analysis for determining adaptive mutations under Cobalt stress

E.coli cells were subjected to gradually increasing levels of Cobalt in a [replifactory](replifactory.com) device.
Ancestral and adapted strains were sequenced using PacBio long read technology.

The bioinformatics pipeline is an effort to automate the dry lab part of experimental evolution.

**Input**: whole genome sequencing reads (fastq, before and after adaptation experiment)

**Output**: annotated mutations, sequence visualization, read mapping (single html, examples below)

When the pipeline is configured correctly and all dependencies are installed, the outputs below are generated automatically.
The [snakefile](https://github.com/catalin-rusnac/CobaltSeq/blob/main/Snakefile) includes tools for quality control and visualizing intermediate results

## Mutations: 

[Vial 1](https://catalin-rusnac.github.io/CobaltSeq/data/processed/igv/hifiasm/Mut1.part_001.vs.WT.part_001.igv.html)
[Vial 2](https://catalin-rusnac.github.io/CobaltSeq/data/processed/igv/hifiasm/Mut2.part_001.vs.WT.part_001.igv.html)
[Vial 3](https://catalin-rusnac.github.io/CobaltSeq/data/processed/igv/hifiasm/Mut3.part_001.vs.WT.part_001.igv.html)
[Vial 4](https://catalin-rusnac.github.io/CobaltSeq/data/processed/igv/hifiasm/Mut4.part_001.vs.WT.part_001.igv.html)
[Vial 5](https://catalin-rusnac.github.io/CobaltSeq/data/processed/igv/hifiasm/Mut5.part_001.vs.WT.part_001.igv.html)
