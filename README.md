# Bioinformatics Pipeline for Experimental Evolution under Cobalt Stress

### Overview

This bioinformatics pipeline automates the genomic analysis of *E. coli* cells subjected to experimental evolution under gradually increasing levels of Cobalt, using a [replifactory](https://replifactory.com) device.

## Features

### 1. **Adaptive Mutation Detection**

- **Goal**: Identify mutations in *E. coli* adapted to Cobalt stress.
- **Sequencing Technology**: PacBio long-read sequencing.
- **Comparison**: Adapted strains vs. ancestral (wild-type) strains.

### 2. **Sequence Visualization**

- **Mutation identification and analysis**: Generates a list of annotated mutations in the adapted strains.
- **Read Mapping**: Visualizes raw sequencing reads aligned to the de-novo assembled genomes for both ancestral and adapted strains.

These visualizations assist in:
- Adjusting [mutation caller parameters](https://github.com/catalin-rusnac/CobaltSeq/blob/bcb1d3675da24ca7d4e2441fabca33684bd2c45d/Snakefile#L287C104-L287C120).
- Validating that the identified mutations are supported by the raw sequencing reads.

### 3. **Quality Control and Intermediate Results Verification**

The pipeline includes quality control tools, ensuring the integrity of the sequencing data and the accuracy of resulting analysis.

---

### Output examples: 

[Vial 1](https://catalin-rusnac.github.io/CobaltSeq/data/processed/igv/hifiasm/Mut1.part_001.vs.WT.part_001.igv.html)
[Vial 2](https://catalin-rusnac.github.io/CobaltSeq/data/processed/igv/hifiasm/Mut2.part_001.vs.WT.part_001.igv.html)
[Vial 3](https://catalin-rusnac.github.io/CobaltSeq/data/processed/igv/hifiasm/Mut3.part_001.vs.WT.part_001.igv.html)
[Vial 4](https://catalin-rusnac.github.io/CobaltSeq/data/processed/igv/hifiasm/Mut4.part_001.vs.WT.part_001.igv.html)
[Vial 5](https://catalin-rusnac.github.io/CobaltSeq/data/processed/igv/hifiasm/Mut5.part_001.vs.WT.part_001.igv.html)

---
### Mutation effect visualizations:

Vial 3 mgtA point mutation:
![Vial 3 mgtA mutation](https://github.com/catalin-rusnac/CobaltSeq/blob/main/3D_protein_visualizations/mgtA815.gif)

Vial 5 mgtA point mutation:
![Vial 5 mgtA mutation](https://github.com/catalin-rusnac/CobaltSeq/blob/main/3D_protein_visualizations/mgtA779.gif)

Vial 3 corA knockout:
![Vial 3 corA knockout](https://github.com/catalin-rusnac/CobaltSeq/blob/main/3D_protein_visualizations/corA_del.gif)

Vial 5 corA knockout: 
![Vial 5 corA knockout](https://github.com/catalin-rusnac/CobaltSeq/blob/main/3D_protein_visualizations/corA_stop.gif)


