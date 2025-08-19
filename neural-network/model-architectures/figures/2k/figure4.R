## Code to create a viz of train loss as a function of complexity
## For the model that trained over ~2,000 sequences

#Importing necessary libraries and data
library(tidyverse)
library(ggplot2)
twok <- read_csv("../../results/results-2k.csv")

#Averaging loss across seeds
twokMean <- twok |>
  group_by(layers, neuronsPerLayer) |>
  summarize(testLoss = mean(testLoss),
         trainLoss = mean(trainLoss))

#Generating plot
plot <- ggplot(data = twokMean) +
  geom_point(aes(x = layers*neuronsPerLayer, y = trainLoss), color = 'blue') +
  geom_smooth(se = FALSE, aes(x = layers*neuronsPerLayer, y = trainLoss), color = 'blue') +
  scale_x_log10() +
  labs(title = "Train Loss and Model Complexity",
       x = "Neurons in model (log scale)",
       y = "Loss",
       subtitle = "Model trained over 2,000 sequences") +
  theme_bw()

ggsave('figure4.png',
       plot,
       width = 6, height = 4)
