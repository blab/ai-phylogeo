source('scripts/utils_remaster.R')

args <- commandArgs(trailingOnly = T)
xml_file_path <- as.character(args[1])
remaster_trees <- as.character(args[2])
remaster_nexus <- as.character(args[3])
remaster_traj <- as.character(args[4])

## Define parameters used for simulation
proba_trans_before_mut <- 0.6 # Probability that transmission occurs before mutation
p_seq <- 0.03 # Probability of sampling an individual for sequencing

## Define generation time parameters
rate_out_of_E <- 0.33 # Rate out which individuals exit the E compartment
rate_out_of_I <- rate_out_of_E # Rate out which individuals exit the I compartment

## Corresponding parameters for a gamma distribution
alpha_GT <- 3.
beta_GT <- 1. / rate_out_of_E

get_mu_from_proba_trans_before_mut <- function(alpha_GT, beta_GT, proba_trans_before_mut){
  return(
    beta_GT * (exp(- 1/ alpha_GT * log(proba_trans_before_mut)) - 1)
  )
}

## Compute mutation rate corresponding to probability that transmission occurs before mutation
## This is a per genome mutation rate
mu_per_genome <- get_mu_from_proba_trans_before_mut(alpha_GT, beta_GT, proba_trans_before_mut) 
sequence_length <- 5000 # Length of simulated genome
mu_per_site <- mu_per_genome / sequence_length

## Demographic parameters
n_demes <- 3
pop_per_deme <- 500000
vec_S_init_per_deme <- rep(pop_per_deme, n_demes)

## Transmission parameters
R0 <- 1.4 # Basic reproduction number
beta_rate_SEIR <- R0 * rate_out_of_I
beta_rate = R0 * rate_out_of_I / (pop_per_deme) # Remaster rate for the force of infection (needs to be scaled by population size)

## Load mixing matrix
mat_proba_migration <- read.csv('input/mixing_matrix.csv', header = F)
mat_proba_migration <- as.matrix(mat_proba_migration)

## Write XML
write_xml(output_file = xml_file_path,
          remaster_nexus_path = remaster_nexus, remaster_trees_path = remaster_trees,
          remaster_traj_path = remaster_traj,
          mat_proba_migration = mat_proba_migration,
          sequence_length = sequence_length, mutation_rate = mu_per_site, 
          beta_rate = beta_rate, 
          rate_out_of_E = rate_out_of_E, 
          rate_out_of_I = rate_out_of_I, 
          p_sample = p_seq,
          n_demes = n_demes, min_sample_per_deme = -1 , 
          vec_S_init_per_deme = vec_S_init_per_deme)
