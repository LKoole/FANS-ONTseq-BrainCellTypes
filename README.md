# Unmasking cell type-specific DNA (hydroxy)methylation within the brain: Integrating nuclei sorting and low-input nanopore sequencing
Epigenetic modifications display cell type-specific patterns, which can be studied within heterogeneous tissues by isolating cell populations. This approach often yields relatively low DNA levels (≤400 ng), thereby limiting the ability to detect epigenetic changes. Here, we outline a methodological workflow for profiling cell type-specific DNA methylation and hydroxymethylation in central nervous system (CNS)-derived samples with limited DNA quantities. Fluorescence-activated nuclei sorting (FANS) was used to isolate neuronal (NeuN+), oligodendroglial (Olig2+), and microglia- and astrocyte-enriched (NeuN-/Olig2-) nuclei from middle frontal gyrus tissue of patients with Alzheimer’s disease, mild cognitive impairment, or non-demented controls. Oxford Nanopore sequencing was used to generate DNA (hydroxy)methylation data. Reduced representation methylation sequencing (RRMS) enriched for regions with high CpG density. Different experimental set-ups were assessed to maximize data throughput. Five analytical models were evaluated using a simulated dataset. FANS of post-mortem brain tissue yielded between 50-800 ng high-integrity DNA. Maximum data throughput was achieved when supplying >45 fmol library and using flow cells with >6500 initial pore count. The dispersion shrinkage for sequencing data and adapted Linear Models for Microarray and Omics Data analyses showed the highest sensitivity and specificity for assessing DNA methylation. The adapted FANS-Nanopore sequencing workflow and associated data analysis pipeline represent a novel method suitable for studying DNA modifications on a cell type-specific level in samples with low DNA content. 


## Repository
The GitHub repository contains all scripts used for this project. 
(1) Evaluation of wet-lab experimental setups 
(2) Preprocessing, quality control, and analyzing data generated through Oxford Nanopore sequencing
(3) Simulating methylation data and evaluating differential methylation methods 

[![](https://img.shields.io/badge/Repository%20license-GPL_3.0-green)](https://github.com/LKoole/ONT-seq-Epi-AD/LICENSE.md) 


## Download files

- Annotation file and add to Anno folder
- GCA_000001405.15_GRCh38_no_alt_analysis_set.fna and add to Reference-files
- GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fai and add to Reference-files


## Contact
Lisa Koole: ([LKoole](https://github.com/LKoole)) - lisa.koole@maastrichtuniversity.nl
