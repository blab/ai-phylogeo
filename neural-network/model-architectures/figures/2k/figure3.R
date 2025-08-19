## Code to create a viz of test loss as a function of complexity
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
  geom_point(aes(x = layers*neuronsPerLayer, y = testLoss), color = 'red') +
  geom_smooth(method = "lm", formula = y ~ x + I(x^2), se = FALSE, 
              aes(x = layers*neuronsPerLayer, y = testLoss), color = 'red') +
  scale_x_log10() +
  labs(title = "Test Loss and Model Complexity",
       x = "Neurons in model (log scale)",
       y = "Loss",
       subtitle = "Model trained over 2,000 sequences") +
  theme_bw()

ggsave('figure3.png',
       plot,
       width = 6, height = 4)
