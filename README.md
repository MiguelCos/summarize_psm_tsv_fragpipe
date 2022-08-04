# Summarize `psm.tsv` file from Fragpipe output after TMT search and normalize intensities

This repo contains two scripts:

## `summarize_psm_tsv_fragpipe.R`:

The `peptide.tsv` file from Fragpipe is summarizing different modified versions of the same peptide into a single peptide sequence. This is not optimal when the modified identity of the peptide is important (i.e. when doing a search using N-terminal acetylation as variable modification).

This simple script is currently only working for one-mixture TMT experiments.

### How to use the script? 

1. Download and unzip this repo into a new folder in your PC.
2. Ideally, initiallize an RStudio project based on the downloaded repo.
3. Add your `psm.tsv` and `annotation.txt` file in the `data/` folder.
  Note: the `annotation.txt` file is the one generated and/or required by TMT integrator for summarization.
4. Open the `summarize_psm_tsv_fragpipe.R` script.
5. Click on `Source`, on the top right corner.

### What does it do? 

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

## `normalize_peptide_abund_by_prot_abund.R`:

After summarizing the peptides and quantitative features from the `psm.tsv` file, this script would help you process the intensity values from the peptide list.

It will normalize the peptide intensities against the protein intensities and apply log2 transformation and median centering to generate a normalized abundance matrix that can be used for downstream analysis.

### How to use the script? 

1. You probably had already downloaded this repo and initialized an RStudio project to execute the previous script... Then:
2. Add the `protein.tsv` file from your Fragpipe search into the `data/` folder.
3. Open the `normalize_peptide_abund_by_prot_abund.R` script.
4. Click on `Source`, on the top right corner.

### What does it do? 

It will:

1. Load the `best_modified_peptides_psm.tsv` file (this is the PSM-summarized peptide list generated with the previous script).
2. Load the `protein.tsv` file.
3. Map the sample names of the annotation file to the column names of the TMT reporter ion intensities in the PSM file and also for the Protein file.
4. Select quantitative columns from the peptide and protein files.
5. Merge the psm file and protein file, to get a peptide-to-protein mapping table with quantitative information for both.
6. Calculate the ratio of `peptide_quant`/`protein_quant` per sample (ratio).
7. Calculate fraction of the peptide intensity normalized by protein intensity (fraction).
8. Apply log2 trans formation and median centering on:
  - Quant table of peptide intensity ratios
  - Quant table of peptide intensity fractions
  - Quant table of non-protein-normalized peptide/PSM intensities.
9. Generate an abundance table, with 1 row per peptide and 1 column per sample for each of these normalized values.
  - `log2_n_MD_peptide_fraction_of_protein_intensity.tsv`
  - `log2_n_MD_peptide_from_psm_wo_prot_norm.tsv`
  - `log2_n_MD_peptide2protein_ratio.tsv`




