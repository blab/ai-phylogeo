## Code to create reliability diagrams for uncalibrated data (20k dataset)

#Importing necessary libraries and data
library(tidyverse)
library(ggplot2)
uncalibrated <- read_csv("../results/testUncalibrated.csv")

#Setting colors for chart
group.colors = c("A" = "#3288BDFF", "B" = "#5E4FA2FF", "C" = "#FDAE61FF")

#Set ideal reliability for comparison
reliabilityIdeal <- data_frame(
  conf = c(0, .1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1),
  acc = c(0, .1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)
)

#Wrangle data for reliability diagrams
BREAKS <- seq(0, 1, by = .1)

reliability <- function(df){
  df <- df |>
    rename(A = prob_class_A, 
           B = prob_class_B, 
           C = prob_class_C) |>
    pivot_longer(A:C, names_to = "group",
                 values_to = "conf") |>
    add_column(acc = 0)
  
  df$acc <- as.integer(df$group == df$actual_label)
  
  df$bin <- cut(df$conf, breaks = BREAKS, labels = FALSE, include.lowest = TRUE)
  
  df <- df |>
    group_by(bin, group) |>
    summarize(conf = mean(conf),
              acc = mean(acc),
              diff = acc-conf) |>
    arrange(bin)
  
  df$bin <- as.numeric(df$bin)
  
  return(df)
}

#Calculate ECE
ECE <- function(df, i){
  dfLong <- df |>
    rename(A = prob_class_A,
           B = prob_class_B,
           C = prob_class_C) |>
    pivot_longer(A:C,
                 names_to = "group",
                 values_to = "conf") |>
    mutate(bin = cut(conf, breaks = BREAKS, labels = FALSE, include.lowest = TRUE),
           acc = as.integer(group == actual_label))
  
  agg <- dfLong |>
    group_by(bin, group) |>
    summarize(conf = mean(conf),
              acc = mean(acc),
              n = n(),
              .groups = 'drop') |>
    mutate(w = n/sum(n),
           ece_term = abs(acc-conf)*w)
  
  
  agg |>
    filter(group == i) |>
    summarize(ECE = sum(ece_term)) |>
    pull(ECE)
}

#Create reliability chart
reliabilityChart <- function(testdf, title){
  test <- reliability(testdf)
  
  ggplot() +
    geom_col(data = test,
             aes(x = (bin/10-0.05), y = acc, fill = group)) +
    geom_line(data = reliabilityIdeal, 
              aes(x = conf, y = acc), linetype = 'dashed', 
              color = 'black') +
    facet_wrap(~group, nrow = 1)+ 
    labs(x = "Confidence", y = "Accuracy",
         colour = NULL,
         title = title,
         subtitle = glue::glue(
           "ECE A: {sprintf('%.4f', ECE(testdf, 'A'))} | ",
           "B: {sprintf('%.4f', ECE(testdf, 'B'))} | ",
           "C: {sprintf('%.4f', ECE(testdf, 'C'))}"
         )) +
    scale_fill_manual(values = group.colors) +
    theme_bw() +
    theme(legend.position = 'none')
}
  
plot <- reliabilityChart(uncalibrated, "Reliability Plots for Uncalibrated Data")

ggsave("figure9.png", plot, width = 11, height = 4)
