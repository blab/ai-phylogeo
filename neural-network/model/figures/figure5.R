## Code to create boxplots of the ratio between predicted probabilities for uncalibrated data (20k dataset)

#Importing necessary libraries and data
library(tidyverse)
library(ggplot2)
library(scales)
uncalibrated <- read_csv("../results/testUncalibrated.csv")

comparisonData <- function(df, title){
  df <- df |>
    mutate(pab = NA,
           pac = NA, 
           pba = NA,
           pbc = NA,
           pca = NA,
           pcb = NA)
  for(i in 1:nrow(df)){
    if(df$actual_label[i] == "A"){
      df$pab[i] <- df$prob_class_A[i] / df$prob_class_B[i]
      df$pac[i] <- df$prob_class_A[i] / df$prob_class_C[i]
    } else if(df$actual_label[i] == "B"){
      df$pba[i] <- df$prob_class_B[i] / df$prob_class_A[i]
      df$pbc[i] <- df$prob_class_B[i] / df$prob_class_C[i]
    } else if(df$actual_label[i] == "C"){
      df$pca[i] <- df$prob_class_C[i] / df$prob_class_A[i]
      df$pcb[i] <- df$prob_class_C[i] / df$prob_class_B[i]
    }
  }
  
  df <- df |>
    pivot_longer(pab:pcb, 
                 names_to = "ratioGroup",
                 values_to = "ratio") |>
    mutate(ratioGroup = recode(ratioGroup,
                               "pab" = "P(A) / P(B)",
                               "pac" = "P(A) / P(C)",
                               "pba" = "P(B) / P(A)",
                               "pbc" = "P(B) / P(C)",
                               "pca" = "P(C) / P(A)",
                               "pcb" = "P(C) / P(B)"
                               ))
  
  df |>
    filter(ratio != 1) |>
    ggplot(aes(x = ratioGroup, y = ratio)) +
    geom_boxplot() +
    facet_wrap(~actual_label, scales = 'free_x', labeller = labeller(actual_label = 
                                                                       c("A" = "True label A",
                                                                         "B" = "True label B",
                                                                         "C" = "True label C"))) +
    scale_y_log10(labels = label_number(accuracy = 0.01)) +
    geom_hline(yintercept = 1, color = "red") +
    theme_bw() +
    labs(x = "",
         y = "Ratio",
         title = title)
  
}

plot <- comparisonData(uncalibrated, "Ratio Between Predicted Probabilities, Uncalibrated Model")

ggsave("figure5.png", plot, width = 9, height = 6)
