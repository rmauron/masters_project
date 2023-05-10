# Chromatin accessibility analysis of spaceflight mouse brains using ArchR
During my masters thesis project I worked with a **single-nuclei ATAC-sequencing dataset** in order to better understand the chromatin dynamics of mouse brain behavior in response to spaceflight experiment.
A multiomic dataset was available with both RNA-seq and ATAC-seq from the same nuclei. In order to better leverage what was happening, I needed a performant tool which focused on ATAC data but which could also incorporate RNA-seq.
I used [ArchR](https://github.com/GreenleafLab/ArchR) which is one of the most updated and performant R-based tool for single-cell ATAC-seq analysis. I don't intend to sell the software but however recommend to read the [documentation](https://www.archrproject.com/) before use.

## Context
The investigation focused on analyzing the behavior of the mouse brain in the context of spaceflight. Notably, intriguing changes resembling neurodegenerative alterations were observed.

## Code availability
The finalized script utilized for generating the presented results can be accessed in the **archr_multi** folder under ```doc > archr_multi > archr_multiome.qmd```

## Reproducibility
Reproducibility of the analysis necessitates a specific version of ```ggplot2``` (version 3.3.6) within the ArchR framework as well as a Python tool, ```MACS2```, which may be inconvenient for some users. To address this concern, a Singularity container is provided on the corresponding [GitHub page](https://github.com/rmauron/Singularity/tree/main/ArchR), accompanied by basic instructions on employing ```Singularity```, accessible [here](https://github.com/rmauron/Singularity/tree/main).

## Written Thesis
The written report detailing the findings of this study will soon be published by **KTH**. The link to access the thesis will be made available in due course.
