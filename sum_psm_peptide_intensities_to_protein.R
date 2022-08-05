# get protein level intensities from summarized PSM file
# Miguel Cosenza
# v 0.1

# required packages ---------------

library(tidyverse)
library(here)

# load data -----------------------

# summarized psm file
psms <- read_tsv(here("output/best_modified_peptides_psm.tsv"))

# load data annotation
annot <- read_delim(here("data/annotation.txt"),
                    col_names = FALSE) %>%
  # set column names to channel and sample respectively
  # this is necessary to select the TMT channel columns from the psm.tsv
  dplyr::rename(channel = X1, sample = X2)

# pre-process annotation file --------------------------

# we need to modify the sample names in the annotation file to match the 
# names of the of the TMT intensities after loading into R.

psm_cols_trim <- colnames(psms) %>% # get column names of psm file
  str_remove("\\..*") # eliminate sufix 

psm_cols <- colnames(psms) # get column names of psm file wo eliminate sufix

# create df mapping trimmed sample names and the ones with sufix
trimed2original_sample <- tibble(sample = psm_cols_trim,
                                 sample_trim = psm_cols)

# correct annotation file to match column names for TMT channel intensities
# in the psm file
annot_2 <- annot %>%
  left_join(., trimed2original_sample) %>% # merge annotation file w untrimmed
  slice(1:nrow(annot)) # keep only sample names present in columns of psm file

# processing ----------------------

## select quant columns from psm.tsv and add suffix

psms_q <- psms %>% 
  filter(`Is Unique` == TRUE) %>% # keep only unique peptides/PSMs
  dplyr::select(`Protein ID`, `Modified Peptide`,
                annot_2$sample_trim) %>%
  group_by(`Protein ID`) %>% 
  summarise_if(is.numeric, sum, na.rm = TRUE)

# generate the file of summed protein intensities from PSM/peptides

write_tsv(x = psms_q,
          file = here("output/abundance_mat_summed_best_psm_intensities_per_protein.tsv"))


