## Code to create a chart of correct predictions per subgroup for calibrated test data (20k dataset)

#Importing necessary libraries and data
library(tidyverse)
library(ggplot2)
calibrated <- read_csv("../results/testVectorScaled.csv")

#Setting colors for the chart
correctColors <- c("0" = "#D53E4FFF", "1" = "#84cd79")

probChart <- function(df, title, subtitle){
  df <- df |>
    group_by(actual_label, predicted_label) |>
    summarize(value = n()) |>
    group_by(actual_label) |>
    mutate(prop = value/sum(value)) |>
    ungroup() |>
    mutate(equal = as.integer(actual_label == predicted_label))
  
  
  ggplot(data = df) +
    geom_col(aes(x = predicted_label, y = prop, fill = as.factor(equal))) +
    facet_wrap(~actual_label, labeller = labeller(actual_label = 
                                                    c("A" = "True label A",
                                                      "B" = "True label B",
                                                      "C" = "True label C"))) +
    scale_y_continuous(labels = scales::percent,
                       limits = c(0, 1)) +
    theme_bw() +
    theme(legend.position = "none") +
    labs(x = "Predicted Label",
         y = "Proportion",
         title = title,
         subtitle = subtitle) +
    scale_fill_manual(values = correctColors)
}

plot <- probChart(vectorScaled, "Proportion of Correct Predictions Across Subgroups", "Calibrated Test Data")

ggsave("figure2.png", plot, width = 9, height = 6)
