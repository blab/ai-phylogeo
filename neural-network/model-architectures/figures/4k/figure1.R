## Code to create a viz of test loss as a function of complexity over various layers
## For the model that trained over ~4,000 sequences

#Importing necessary libraries and data
library(tidyverse)
library(ggplot2)
library(RColorBrewer)
fourk <- read_csv("../../results/results-4k.csv")

#Averaging loss across seeds
fourkMean <- fourk |>
  group_by(layers, neuronsPerLayer) |>
  summarize(testLoss = mean(testLoss),
         trainLoss = mean(trainLoss))

#Generating plot
plot <- ggplot(data = filter(fourkMean, layers > 1)) +
  geom_line(aes(x = neuronsPerLayer, 
                y = testLoss, 
                color = as.factor(layers))) +
  geom_point(aes(x = neuronsPerLayer, 
                 y = testLoss, 
                 color = as.factor(layers)))+
  scale_x_log10() +
  labs(color = "Layers",
       x = "Neurons per layer (log scale)",
       y = "Loss",
       title = "Test loss at various architectures",
       subtitle = "Model trained on 4,000 sequences") +
  scale_color_manual(values = brewer.pal(name = "Reds", n = 6)[2:9]) +
  theme_bw()

ggsave('figure1.png', plot, width = 6, height = 4)

