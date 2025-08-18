## Code to create a boxplot of predicted probabilities for calibrated test data (20k dataset)

#Importing necessary libraries and data
library(tidyverse)
library(ggplot2)
calibrated <- read_csv("../results/testVectorScaled.csv")

#Setting colors for the chart
correctColors <- c("FALSE" = "#D53E4FFF", "TRUE" = "#84cd79")

boxChartDots <- function(df, title, subtitle){
  df_long <- df |>
    rename("A" = prob_class_A,
           "B" = prob_class_B,
           "C" = prob_class_C) |>
    pivot_longer(cols = A:C,
                 names_to = "group",
                 values_to = "prob")
  
  ggplot(data = df_long, aes(x = group, y = prob)) +
    geom_boxplot(outliers = FALSE, 
                 aes(fill = as.factor(group == actual_label))) +
    facet_wrap(~actual_label, labeller = labeller(actual_label = 
                                                    c("A" = "True label A",
                                                      "B" = "True label B",
                                                      "C" = "True label C"))) +
    theme_bw() +
    labs(x = "Subgroups",
         y = "Predicted Probability",
         title = title,
         subtitle = subtitle) +
    scale_fill_manual(values = correctColors) +
    theme(legend.position = "none") +
    geom_hline(yintercept = .33, linetype = 'dashed')
}

plot <- boxChartDots(calibrated, "Distribution of Predicted Probabilities", "Calibrated Test Data")

ggsave("figure4.png", plot, width = 9, height = 6)
