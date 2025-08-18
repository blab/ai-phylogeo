library(jsonlite)
library(tidyverse)
# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
input_traits_json <- as.character(args[1])
output_traits_tsv <- as.character(args[2])

# Load the traits.json file
traits_data <- fromJSON(input_traits_json)

# Extract subgroup confidence data from the nodes
df_for_tsv <- Reduce('bind_rows', lapply(names(traits_data$nodes), FUN = function(strain_name){
  node_data <- traits_data$nodes[[strain_name]]
  confidence <- node_data$subgroup_confidence
  c('strain' = strain_name, 
    'proba_A' = ifelse(is.null(confidence$A), 0.0, confidence$A),
    'proba_B' = ifelse(is.null(confidence$B), 0.0, confidence$B),
    'proba_C' = ifelse(is.null(confidence$C), 0.0, confidence$C)
    )
}))

# Write to TSV file
write.table(df_for_tsv, output_traits_tsv,
            sep = "\t", row.names = FALSE, quote = FALSE)
