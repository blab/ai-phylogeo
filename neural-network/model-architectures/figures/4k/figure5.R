## Code to create a viz comparing the # of neurons per layer that minimizes loss for training and test groups
## For the model that trained over ~4,000 sequences

#Importing necessary libraries and data
library(tidyverse)
library(ggplot2)
fourk <- read_csv("../../results/results-4k.csv")

#Extracting minimum loss across both seeds
minTrain <- fourk |>
  filter(layers != 1) |>
  group_by(layers) |>
  slice_min(trainLoss)

minTest <- fourk |>
  filter(layers != 1) |>
  group_by(layers) |>
  slice_min(testLoss)

#Establishing colors for chart
colors <- c("Train" = "blue", "Test" = "red")

#Generating plot
plot <- ggplot() +
  geom_point(aes(x = layers, y= neuronsPerLayer, color = "Train"), data = minTrain) +
  geom_line(aes(x = layers, y= neuronsPerLayer, color = 'Train'), data = minTrain) +
  geom_point(aes(x = layers, y= neuronsPerLayer, color = 'Test'), data = minTest) +
  geom_line(aes(x = layers, y= neuronsPerLayer, color = 'Test'), data = minTest) +
  theme_bw() +
  labs(x = "Layers",
       y = "Neurons per layer",
       title = "Number of neurons that minimizes train/test loss, per layer",
       subtitle = "Model trained on 4,000 sequences", 
       color = "") +
  scale_color_manual(values = colors) +
  theme_bw()

ggsave('figure5.png',
       plot,
       width = 6, height = 4)
