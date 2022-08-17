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

# function to calculate fraction of usable intensity by purity -----

refine_int <- function(col){
  # extract columns for peptide id, purity and quant
  # col should be a vector with names corresponding to the names of 
  # the quant columns from the summarized psm.tsv file.
  # this was taken care of during summarization of annotation file
  df_q <- psms %>%
    dplyr::select(`Modified Peptide`, Purity, starts_with(col)) 
  
  # get column names of the previous data frame
  df_q_names <- colnames(df_q)
  
  # exection of purification
  df_q_rat <- df_q %>%
    #multiply the intensity times the purity, for each quant column
    mutate({{col}} := .data[[{{col}}]] * .data[["Purity"]]) %>%
    # then select only the quant column used for the purification
    dplyr::select(all_of(col))
  
  return(df_q_rat)
  
}

# purification step -----

# apply the purification function for each quantitative column
# each list element in the output would represent a purified column of the 
# psm summarized intensity matrix
purified_ <- purrr::map(.x = annot_2$sample_trim,
                        .f = refine_int)

# then we merge each of the elements of the list above
purified_quant <- bind_cols(purified_) %>% 
  # add information about peptide and proteins to the quant data
  mutate(`Modified Peptide` = psms$`Modified Peptide`,
         `Protein ID` = psms$`Protein ID`) %>%
  mutate(protein_id_modif_pep = paste(`Protein ID`, `Modified Peptide`,
                                      sep = "_")) %>%
  dplyr::relocate(protein_id_modif_pep) %>%
  dplyr::select(-c(`Modified Peptide`, `Protein ID`))


# median centering and log2 transformation -----

# median centering and log2 transformation of purified quants ----
mat_ratios_pure <- purified_quant %>% 
  column_to_rownames("protein_id_modif_pep") %>% 
  as.matrix()

mat_ratios_pure[mat_ratios_pure == 0] <- NA

log2_mat_pure <- mutate_all(as.data.frame(mat_ratios_pure), 
                            log2)

scaled_mat_pure <- scale(log2_mat_pure,
                         scale = F,
                         center = apply(log2_mat_pure, 2, median, 
                                        na.rm = TRUE) - median(as.matrix(log2_mat_pure), 
                                                               na.rm = TRUE)) %>%
  abs(.) # log2 values of fractions are negative; take absolute values

scaled_mat_pure_df <- scaled_mat_pure %>% 
  as.data.frame() %>%
  rownames_to_column("protein_peptidemod")

# median centering and log2 transformation of non purified quants ---- 

psms_q2 <- psms %>% 
  dplyr::select(`Protein ID`, `Modified Peptide`,
                annot_2$sample_trim) %>% 
  mutate(protein_id_modif_pep = paste(`Protein ID`, `Modified Peptide`,
                                      sep = "_")) %>%
  dplyr::relocate(protein_id_modif_pep) %>%
  dplyr::select(-c(`Modified Peptide`, `Protein ID`)) 

mat_nopure <- psms_q2 %>% 
  column_to_rownames("protein_id_modif_pep") %>% 
  as.matrix()

mat_nopure[mat_nopure == 0] <- NA

log2_mat_nopure <- mutate_all(as.data.frame(mat_nopure), 
                              log2)

scaled_mat_nopure <- scale(log2_mat_nopure,
                           scale = F,
                           center = apply(log2_mat_nopure, 2, median, 
                                          na.rm = TRUE) - median(as.matrix(log2_mat_nopure), 
                                                                 na.rm = TRUE)) %>%
  abs(.) # log2 values of fractions are negative; take absolute values

scaled_mat_nopure_df <- scaled_mat_nopure %>% 
  as.data.frame() %>%
  rownames_to_column("protein_peptidemod")

# export tabular files ----

# standardized and purified abundance matrix

write_tsv(x = scaled_mat_pure_df,
          file = here("output/abundance_mat_standardized_n_purified.tsv"))

# standardized and non-purified abundance matrix

write_tsv(x = scaled_mat_nopure_df,
          file = here("output/abundance_mat_standardized_n_non_purified.tsv"))


