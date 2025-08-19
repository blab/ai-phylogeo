## Code to create a viz of train loss as a function of complexity over various layers
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

#Setting colors for plot 
#Taken from RColorBrewer Blues
blues.pal = c('#6BAED6', "#4292C6", "#2171B5", "#08519C", "#08306B")

#Generating plot
plot <- ggplot(data = filter(twokMean, layers > 1)) +
  geom_line(aes(x = neuronsPerLayer, 
                y = trainLoss, 
                color = as.factor(layers))) +
  geom_point(aes(x = neuronsPerLayer, 
                 y = trainLoss, 
                 color = as.factor(layers)))+
  scale_x_log10() +
  labs(color = "Layers",
       x = "Neurons per layer (log scale)",
       y = "Loss (log scale)",
       title = "Train loss at various architectures",
       subtitle = "Model trained on 2,000 sequences") +
  scale_color_manual(values = blues.pal) +
  scale_y_log10(labels = trans_format("log10", math_format(10^.x))) +
  theme_bw()

ggsave('figure2.png', plot, width = 6, height = 4)
