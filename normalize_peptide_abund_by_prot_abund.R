# Normalize peptide intensities from summarized psm.tsv to protein intensities
# Miguel Cosenza
# v 0.1

# required packages ----

library(tidyverse)
library(here)

# required data ----------------------------------------

# summarized psm file
psms <- read_tsv(here("output/best_modified_peptides_psm.tsv"))

# load data annotation
annot <- read_delim(here("data/annotation.txt"),
                    col_names = FALSE) %>%
  # set column names to channel and sample respectively
  # this is necessary to select the TMT channel columns from the psm.tsv
  dplyr::rename(channel = X1, sample = X2)

# load protein-level data
prots <- read_tsv(here("data/protein.tsv"))

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

# get the names of the quant columns from protein.tsv 
# 

prot_cols_trim <- colnames(prots) %>% # get column names of psm file
  str_remove("\\..*") # eliminate sufix 

prot_cols <- colnames(prots) # get column names of psm file wo eliminate sufix

# create df mapping trimmed sample names and the ones with sufix
# for protein table
trimed2original_prot <- tibble(sample = prot_cols_trim,
                               sample_prot_trim = prot_cols)

# add info about column names for protein.tsv to the annotation file

annot_3 <- annot %>%
  left_join(., trimed2original_prot) %>% # merge annotation file w untrimmed
  slice(1:nrow(annot)) # keep only sample names present in columns of psm file

annot_4 <- bind_cols(annot_2, 
                     sample_prot_trim = annot_3$sample_prot_trim) %>% 
  dplyr::rename(sample_psm_trim = sample_trim)


# processing files -------------------------------------

## select quantitative columns from  protein.tsv 

prots_q <- prots %>% 
  dplyr::select(`Protein ID`, annot_4$sample_prot_trim) %>%
  rename_at(.vars = vars(annot_4$sample_prot_trim), 
            .funs = function(x) paste0(x, "_prot")) # add prot suffix to columns

## select quant columns from psm.tsv and add suffix

psms_q <- psms %>% 
  dplyr::select(`Protein ID`, `Modified Peptide`,
                annot_4$sample_psm_trim) %>%
  rename_at(.vars = vars(annot_4$sample_psm_trim), 
            .funs = function(x) paste0(x, "_psm")) %>% # add prot suffix to columns
  mutate(protein_id_modif_pep = paste0(`Protein ID`,"_",`Modified Peptide`))


## merge psm file with protein file 

psm_n_prots <- left_join(psms_q, prots_q) %>%
  dplyr::relocate(protein_id_modif_pep)

## function to get the ratio of intensity of PSM/Protein  -----
#col

pept2prot_ratios <- function(col){
  
  df_q <- psm_n_prots %>%
    dplyr::select(protein_id_modif_pep, starts_with(col)) %>%
    dplyr::relocate(protein_id_modif_pep,ends_with("psm"), ends_with("prot"))
  
  df_q_names <- colnames(df_q)
  
  df_q_rat <- df_q %>%
    mutate({{col}} := .data[[df_q_names[2]]] / .data[[df_q_names[3]]]) %>% 
    mutate("fraction_int_psm2prot_{col}" := .data[[{{col}}]]*.data[[df_q_names[2]]])
  
  df_q_rat <- df_q_rat %>%
    dplyr::select(-c(protein_id_modif_pep, ends_with("psm"), ends_with("prot")))
  
  return(df_q_rat)
  
}

## PSM intensity/protein ration calculation per sample ----

list_pept2prot_rat <- purrr::map(.x = annot_4$sample,
                                 .f = pept2prot_ratios)

pept2prot_norm1 <- bind_cols(list_pept2prot_rat) %>% 
  mutate(protein_id_modif_pep = psm_n_prots$protein_id_modif_pep) %>%
  dplyr::relocate(protein_id_modif_pep)

pept2prot_ratios <- pept2prot_norm1 %>%
  dplyr::select(-starts_with("fraction_int_"))

pep2norm1_names <- colnames(pept2prot_ratios)[-1]

pept2prot_ratios <- pept2prot_ratios %>%
  rename_at(.vars = pep2norm1_names, 
            .funs = function(x) paste0("psm2prot_ratio_",x)) # add prot suffix to columns

pept2prot_norm_ratio <- pept2prot_norm1 %>%
  dplyr::select(protein_id_modif_pep, starts_with("fraction_int_"))

## normalizations ---------------------------------------------------

### log2 and median centering of PSM / protein intensity ratios -----

mat_ratios <- pept2prot_ratios %>% 
  column_to_rownames("protein_id_modif_pep") %>% 
  as.matrix()

log2_mat_rat <- mutate_all(as.data.frame(mat_ratios), 
                           log2)

scaled_mat_rat <- scale(log2_mat_rat,
                        scale = F,
                        center = apply(log2_mat_rat, 2, median, 
                                       na.rm = TRUE) - median(as.matrix(log2_mat_rat), 
                                                              na.rm = TRUE)) %>%
  abs(.) # log2 values of fractions are negative; take absolute values

boxplot(mat_ratios,
        main = "Boxplot of Peptide to Protein Intensity ratios",
        sub  = "Before log2 and median centering")

boxplot(scaled_mat_rat,
        main = "Boxplot of Peptide to Protein Intensity ratios",
        sub  = "After log2 and median centering")

### log2 and median centering of fraction of PSM / protein intensity ratios -----

mat_ratios2 <- pept2prot_norm_ratio %>% 
  column_to_rownames("protein_id_modif_pep") %>% 
  as.matrix()

log2_mat_rat2 <- mutate_all(as.data.frame(mat_ratios2), 
                           log2)

scaled_mat_rat2 <- scale(log2_mat_rat2,
                        scale = F,
                        center = apply(log2_mat_rat2, 2, median, 
                                       na.rm = TRUE) - median(as.matrix(log2_mat_rat2), 
                                                              na.rm = TRUE)) %>%
  abs(.) # log2 values of fractions are negative; take absolute values

boxplot(mat_ratios2,
        main = "Boxplot of Fraction of Peptide intensity over protein",
        sub  = "Before log2 and median centering")

boxplot(scaled_mat_rat2,
        main = "Boxplot of Fraction of Peptide intensity over protein",
        sub  = "After log2 and median centering")


### log2 and median centering of PSM intensities without protein normalization ------

psms_q2 <- psms_q %>%
  dplyr::select(-c(`Protein ID`, `Modified Peptide`))

mat_ratios3 <- psms_q2 %>% 
  column_to_rownames("protein_id_modif_pep") %>% 
  as.matrix()

log2_mat_rat3 <- mutate_all(as.data.frame(mat_ratios3), 
                            log2)

scaled_mat_rat3 <- scale(log2_mat_rat3,
                         scale = F,
                         center = apply(log2_mat_rat3, 2, median, 
                                        na.rm = TRUE) - median(as.matrix(log2_mat_rat3), 
                                                               na.rm = TRUE)) %>%
  abs(.) # log2 values of fractions are negative; take absolute values

boxplot(mat_ratios3,
        main = "Boxplot of peptide/PSM intensities - No protein-level correction",
        sub  = "Before log2 and median centering")

boxplot(scaled_mat_rat3,
        main = "Boxplot of peptide/PSM intensities - No protein-level correction",
        sub  = "After log2 and median centering")


### generate tabular outputs of the different types of normalization approaches ----

# normalization of the psm/protein ratio

df_scaled_mat_rat <- as.data.frame(scaled_mat_rat) %>% 
  rownames_to_column("protein_id_modif_pep") %>%
  separate(col = protein_id_modif_pep,
           into = c("Protein ID", "Modified Peptide"),
           sep = "_")

write_tsv(df_scaled_mat_rat,
          file = here("output/log2_n_MD_peptide2protein_ratio.tsv"))

# normalization of the psm intensity fraction of protein intensity

df_scaled_mat_fraction <- as.data.frame(scaled_mat_rat2) %>% 
  rownames_to_column("protein_id_modif_pep") %>%
  separate(col = protein_id_modif_pep,
           into = c("Protein ID", "Modified Peptide"),
           sep = "_")

write_tsv(df_scaled_mat_fraction,
          file = here("output/log2_n_MD_peptide_fraction_of_protein_intensity.tsv"))


# normalization of the psm intensity without protein-level correction

df_scaled_mat_psm <- as.data.frame(scaled_mat_rat3) %>% 
  rownames_to_column("protein_id_modif_pep") %>%
  separate(col = protein_id_modif_pep,
           into = c("Protein ID", "Modified Peptide"),
           sep = "_")

write_tsv(df_scaled_mat_psm,
          file = here("output/log2_n_MD_peptide_from_psm_wo_prot_norm.tsv"))





