# Immuno-cross

Search possible cross-reactivity between the human nervous system and Epstein-Barr Virus (EBV) proteins.

## Description

*Immuno-cross* is a set of Python scripts that identifies critical contact residues of the T-cell receptor (TCR) based on sequence identity, developed by [Helmut K. A. Patrocínio](https://github.com/helkennedy) and [Tayná Fiúza](https://github.com/fiuzatayna), under the supervision of [João Paulo M.S. Lima](https://github.com/jpmslima). Our approach helped to reveal that several peptides derived from the nervous system and Epstein-Barr virus (EBV) proteins share identical residues at these critical contact points. This suggests the possibility of cross-reactivity between them. The pipeline can be used to search for nonamer pairs from other sources and other human proteins.

## Folders content

The folder "files" contains data obtained from the consulted databases.

The folder "prediction" contains data obtained from HLA-binding predictions.

> *Some files of this folder are compressed (.tar.gz) due to size limitations.*

The folder "scripts" contains scripts created for data extraction and processing, as well as to compare the nonamer identity, as follows:

* ***Haplotypes.ipynb and Haplotypes.py*** - getting alleles and haplotypes
* ***prediction_treatment.ipynb and prediction_treatment.py*** - process data from allele-binding prediction, and generate nonamers fasta files 
* ***nonamers_identity_analysis.ipynb and nonamers_identity_analysis.py***- read nonamers fasta files, compare identity, and generate .xlsx files with the nonamers identity result.
* ***upset_plot.ipynb and upset_plot.py*** - read .xlsx files with the nonamers' identity results and generate upset plots.

## Dependencies

To execute the scripts, the following packages/libraries must be installed in your system/notebook:

- biopython>=1.79.

- matplotlib>=3.5.1.

- numpy>=1.22.1.

- openpyxl==3.0.9.

- pandas>=1.4.2.

- seaborn>=0.11.2.

- UpSetPlot==0.6.1.

- XlsxWriter==3.0.3.

- matplotlib-venn==0.11.7
  bs4==0.0.1.

They are also listed in the file requirements.txt, in the scripts folder.

## The Pipeline

The scheme below describes the methodology workflow that Immuno-cross uses to identify cross-identity peptides. 

![The Immuno-cross pipeline](pipeline-new.jpg)
**The Immuno-cross pipeline:** **A**. Data retrieving steps (blue), input data, and subsequent analysis (gray). **B**. The search for cross-reactivity of the peptides from human nervous system proteins against peptides from Epstein-Barr virus. **C**. The pipeline uses four rules for selecting nonamers with relevant identities interacting with the TCR. It starts with a minimal identity of 44% between the residues of two given nonamers and identical residues at position 5 (P5) as base criteria and the cross-identity of nonamers that bind to the same HLA allele or HLA alleles contained in a common haplotype. **D**. The three additional rules used to select nonamers where TCR-relevant cross-identity may occur. Criterion 1: identical residues at P2, P3, and P8 (P2-P3-P5-P8). Criterion 2: two identical residues at P2, P3 or P8, and at least one identical residue at P4, P6 or P7 (P5 - 2 x P2/P3/P8 -  1 x P4/P6/P7). Criterion 3: identical residues at P4, P6 and P7 and at least one identical residue at P2, P3 and P8 (P5 - 1 x P2/P3/P8 -P4-P6-P7).

The pipeline is organized in the following order:  

1. Run '*haplotypes*' (which generates files for '*nonamers_identity_analysis*').

2. Then '*prediction_treatment*' (which generates files for '*nonamers_identity_analysis*').

3. Now execute '*nonamers_identity_analysis*' (which generates files for '*upset_plot*').

4.  Run '*upset_plot*'.

>  *Ps.: All files are included, so you can run 'nonamers_identity_analysis' or 'upset_plot' scripts directly.*

We used the packages matplotlib ([https://matplotlib.org/](https://matplotlib.org/)) (version 3.5) and seaborn ([https://seaborn.pydata.org/](https://seaborn.pydata.org/)) (version 0.11.2) to build the graphs. For visualizing the identity analysis of the nonamers, we used the package UpSetPlot ([UpSetPlot documentation — upsetplot 0.8.0 documentation](https://upsetplot.readthedocs.io/en/stable/)) (v. 0.6.1).

## Citation

Preprint at BioRxiv:

[Patrocinio et al. 2023]().
