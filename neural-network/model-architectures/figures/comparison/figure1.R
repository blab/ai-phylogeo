## Code to create a viz comparing train loss as a function of complexity
## For the models that trained over ~2,000 and ~4,000 sequences

#Importing necessary libraries and data
library(tidyverse)
library(ggplot2)
twok <- read_csv("../../results/results-2k.csv")
fourk <- read_csv("../../results/results-4k.csv")

#Averaging loss across seeds
twokMean <- twok |>
  group_by(layers, neuronsPerLayer) |>
  summarize(testLoss = mean(testLoss),
            trainLoss = mean(trainLoss))

fourkMean <- fourk |>
  group_by(layers, neuronsPerLayer) |>
  summarize(testLoss = mean(testLoss),
         trainLoss = mean(trainLoss))

#Generating plot
plot <- ggplot() +
  geom_point(data = twokMean, aes(x = layers*neuronsPerLayer, y = trainLoss, color = '2,000')) +
  geom_smooth(data = twokMean, se = FALSE, 
              aes(x = layers*neuronsPerLayer, y = trainLoss, color = '2,000')) +
  geom_point(data = fourkMean, aes(x = layers*neuronsPerLayer, y = trainLoss, color = '4,000')) +
  geom_smooth(data = fourkMean, se = FALSE, aes(x = layers*neuronsPerLayer, y = trainLoss, color = '4,000')) +
  scale_x_log10() +
  labs(title = "Train Loss and Model Complexity",
       x = "Neurons in model (log scale)",
       y = "Loss",
       color = 'Pool') +
  theme_bw() + 
  scale_color_manual(values = c(
    '4,000' = '#08306B',
    '2,000' = '#6BAED6'
  )) +
  ylim(0, .65)

ggsave('figure1.png', plot, width = 6, height = 4)
