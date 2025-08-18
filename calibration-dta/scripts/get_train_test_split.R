library(tidyverse)

args <- commandArgs(trailingOnly = T)
file_path_input_metadata <- as.character(args[1])
file_path_output_metadata <- as.character(args[2])
prop_train <- as.numeric(args[3])
random_seed <- as.numeric(args[4])

# Read input metadata file
metadata_input <- read_csv(file_path_input_metadata)

# Draw train / test status
set.seed(random_seed)
is_train <- rbinom(n = nrow(metadata_input), size = 1, prob = prop_train)

# Create output metadata file
metadata_for_auspice <- metadata_input %>% 
  mutate(is_train = is_train) %>% 
  mutate(type = c('Test', 'Train')[is_train + 1]) %>% 
  select(- is_train) %>% 
  rename(true_subgroup = subgroup) %>% 
  mutate(subgroup = ifelse(type == 'Test', '?', true_subgroup))


# Save output metadata file
write_tsv(metadata_for_auspice, file_path_output_metadata)