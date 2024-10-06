# Interactive phylogeny of early SARS-CoV-2 in 2024 study by Crits-Christoph et al

## Overview
This repository builds a phylogeny of the sequences from [Crits-Christoph et al, Cell (2024)](https://www.cell.com/cell/fulltext/S0092-8674(24)00901-2), which reports a phylogeny of SARS-CoV-2 sequences collected prior to Jan-20-2020.

The phylogeny is rooted on lineage A like the image in Figure 1 of [Crits-Christoph et al, Cell (2024)](https://www.cell.com/cell/fulltext/S0092-8674(24)00901-2), but shows the full set of sequences that paper refers to analyzing, not just the subset shown in Figure 1A.

Note that the best way to root the SARS-CoV-2 phylogeny is still an open topic with substantial uncertainty (see [Pipes et al (2021)](https://academic.oup.com/mbe/article/38/4/1537/6028993)).
There are also open questions about the best way to subsample, de-duplicate, and quality control early SARS-CoV-2 sequences. Further, the high similarity of the sequences mean there is also just limited statistical support to reliably draw conclusions about some aspects of the phylogeny (see the [paper "Phylogenetic analysis of SARS-CoV-2 data is difficult" by Morel et al (2020)](https://academic.oup.com/mbe/article/38/5/1777/6030946)).
So recall the conclusion of [Pipes et al (2021)](https://academic.oup.com/mbe/article/38/4/1537/6028993) that _"These results suggest that phylogenetic evidence alone is unlikely to identify the origin of the SARS-CoV-2 virus and we caution against strong inferences regarding the early spread of the virus based solely on such evidence."_

## Details

The `conda` environment to run the analysis is in [environment.yml](environment.yml).

The input data are in [./data/](data):
  - [data/jointWHO_sequences.csv](data/jointWHO_sequences.csv): sequences from the Tables 6 and 7 of the [joint WHO-China report](https://www.who.int/publications/i/item/who-convened-global-study-of-origins-of-sars-cov-2-china-part), which are stated to be the only sequences from patients with onset date before Dec-31-2019. Note that I added NMDC60013002-09 and NMDC60013002-06 for consistency with the annotations of the samples in Table S1 of [Crits-Christoph et al, Cell (2024)](https://www.cell.com/cell/fulltext/S0092-8674(24)00901-2)
  - [data/Table_S1.xlsx](data/Table_S1.xlsx): Table S1 from [Crits-Christoph et al, Cell (2024)](https://www.cell.com/cell/fulltext/S0092-8674(24)00901-2)
  - [data/extract_accessions.ipynb](data/extract_accessions.ipynb): Jupyter notebook that extracts lists of accessions of the SARS-CoV-2 sequences from Table S1. It should be run manually to create the following files:
    - [data/seq_metadata.csv](data/seq_metadata.csv): all metadata extracted from the table of accessions.
    - [data/custom_accessions.csv](data/custom_accessions.csv): metadata for sequences with no accession (custom built)
    - [data/gisaid_accessions.csv](data/gisaid_accessions.csv): metadata for sequences with GISAID accessions
    - [data/ngdc_accessions.csv](data/ngdc_accessions.csv): metadata for sequences with NGDC accessions
    - [data/genbank_accessions.csv](data/genbank_accessions.csv): metadata for sequences with Genbank accessions
  - The following files are downloaded manually on Oct-5-2024 to match the accessions listed immediately above:
    - [data/custom_sequences.fa](data/custom_sequences.fa): downloaded from [here](https://github.com/sars-cov-2-origins/huanan-market-environment/tree/main/sars2_phylogenetics/HSM_sequences) in the GitHub repo of [Crits-Christoph et al, Cell (2024)](https://www.cell.com/cell/fulltext/S0092-8674(24)00901-2). Note that the sequence they provide for sample B5 is all `N` nucleotides, so they must have some error as they still show this sequence on the tree in their paper.
    - [data/gisaid_sequences.fa](data/gisaid_sequences.fa): note my GISAID query turned up 782 sequences not the 783 GISAID accessions given by [Crits-Christoph et al, Cell (2024)](https://www.cell.com/cell/fulltext/S0092-8674(24)00901-2). The accession `EPI_ISL_1143994` for strain *hCoV-19/Germany/BY-RKI-I-013858/2020* is not found on GISAID. **Due to GISAID data sharing restrictions, this file is not tracked in the GitHub repo.**
    - [data/ngdc_sequences.fa](data/ngdc_ssequences.fa): sequences from NGDC
    - [data/genbank_sequences.fa](data/genbank_sequences.fa): sequences from Genbank.


The pipeline to build the tree mostly uses the [augur](https://docs.nextstrain.org/projects/augur) pipeline, and is run by the `snakemake` file in [Snakefile](Snakefile).

All the files created by the pipeline are placed in `./results/` subdirectory which is not tracked in this GitHub repo due to GISAID data-sharing restrictions.
The final tree is in [auspice/build_CritsChristoph_SARS2_seqs_tree.json](auspice/build_CritsChristoph_SARS2_seqs_tree.json), and can be viewed with Nextstrain using the [instructions here](https://docs.nextstrain.org/en/latest/guides/share/community-builds.html) at the link [https://nextstrain.org/community/jbloom/build_CritsChristoph_SARS2_seqs_tree](https://nextstrain.org/community/jbloom/build_CritsChristoph_SARS2_seqs_tree).
