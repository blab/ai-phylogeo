library(tidyverse)

args <- commandArgs(trailingOnly = T)
file_path_metadata_auspice <- as.character(args[1])
file_path_traits <- as.character(args[2])

file_path_metadata_auspice <- '../results_larger/metadata_for_dta_calibration.tsv'
file_path_traits <- '../results_larger/traits.tsv'

# Load metadata
metadata <- read_tsv(file_path_metadata_auspice)
inferred_traits <- read_tsv(file_path_traits)

#
df_inference_accuracy <- metadata %>%
  filter(type == 'Test') %>% 
  left_join(inferred_traits) %>% 
  pivot_longer(cols = 6:8, names_to = 'proba_group', values_to = 'proba', names_prefix = 'proba_') %>% 
  group_by(strain) %>% 
  filter(proba == max(proba)) %>% 
  mutate(accurate_inference = (proba_group == true_subgroup)) %>% 
  ungroup()

df_inference_accuracy %>% 
  ungroup() %>% 
  summarise(accuracy = sum(accurate_inference)/n())
