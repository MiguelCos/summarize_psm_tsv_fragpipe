# Summarize `psm.tsv` file from Fragpipe output after TMT search

Small script to summarize peptides + modifications from the psm.tsv file (Fragpipe output) after TMT experiment.

The `peptide.tsv` file from Fragpipe is summarizing different modified versions of the same peptide into a single peptide sequence. This is not optimal when the modified identity of the peptide is important (i.e. when doing a search using N-terminal acetylation as variable modification).

This simple script is currently only working for one-mixture TMT experiments.

## How to use the script? 

1. Download and unzip this repo into a new folder in your PC.
2. Ideally, initiallize an RStudio project based on the downloaded repo.
3. Add your `psm.tsv` and `annotation.txt` file in the `data/` folder.
  Note: the `annotation.txt` file is the one generated and/or required by TMT integrator for summarization.
4. Open the `summarize_psm_tsv_fragpipe.R` script.
5. Click on `Source`, on the top right corned.

## What does it do? 

It will:

1. Load the `psm.tsv` file.
2. Load the `annotation.txt` file
3. Map the sample names of the annotation file to the column names of the TMT reporter ion intensities in the PSM file.
4. Sum the total reporter ion intensities.
5. Filter the PSM file to keep only the best PSM per modified peptide under the following criteria:
  - Only PSMs with PeptideProbability > 0.9.
  - Only PSMs with summed reporter ion intensities > 0.
  - If one modified peptide has several redundant PSMs, select the PSM with best Purity.
  - If one modified peptide has severak reduntant PSMs with the same purity, then select the PSM with the highest summed reporter ion intensity.
6. Generate a new filtered `psm.tsv` file, located in the `output/` folder and called `best_modified_peptides_psm.tsv`.


