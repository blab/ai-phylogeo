## Code to create density plots of predicted probabilities for calibrated data (20k dataset)

#Importing necessary libraries and data
library(tidyverse)
library(ggplot2)
uncalibrated <- read_csv("../results/testUncalibrated.csv")

#Setting colors for chart
group.colors = c("A" = "#3288BDFF", "B" = "#5E4FA2FF", "C" = "#FDAE61FF")

densityGrid <- function(df, title){
  df <- df |>
    pivot_longer(cols = prob_class_A:prob_class_C,
                 names_to = "group",
                 values_to = "prob")
  
  ggplot(data = df, aes(x = prob, color = actual_label)) +
    geom_density(linewidth = 1) +
    facet_grid(rows = vars(group),
               cols = vars(actual_label),
               scales = 'free_y',
               labeller = labeller(group = as_labeller(c(
                 "prob_class_A" = "P(A)",
                 "prob_class_B" = "P(B)",
                 "prob_class_C" = "P(C)"
               )),
               actual_label = as_labeller(c(
                 "A" = "True label A",
                 "B" = "True label B",
                 "C" = "True label C"
               ))
               )) +
    scale_color_manual(values = group.colors) +
    theme_bw(base_family = "Helvetica") +
    labs(x = "Predicted Probability",
         y = "Density",
         title = title) +
    theme(legend.position = "none")
}

plot <- densityGrid(uncalibrated, "Density of Predicted Probabilities, Uncalibrated Data")

ggsave("figure7.png", plot, width = 9, height = 7)
