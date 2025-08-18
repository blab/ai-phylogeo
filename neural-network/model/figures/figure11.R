## Code to create a visualization of the decrease in ECE with calibration

#Importing necessary libraries and data
library(tidyverse)
library(ggplot2)

eceData <- data.frame(group = c("A", "B", "C", "A", "B", "C", "Overall", "Overall"),
                      ECE = c(0.03536857, 0.01686004, 0.02833649, 0.01437688, 0.01825119, 0.01908578, 0.05171385, 0.0805651),
                      calibrated = c(0, 0, 0, 1, 1, 1, 1, 0))

plot <- ggplot(data = eceData, aes(x = as.factor(calibrated), y = ECE, fill = as.factor(calibrated))) +
  geom_col() + 
  facet_wrap(~group, nrow = 1, 
             labeller = labeller(group = 
                                   c("A" = "True label A",
                                     "B" = "True label B",
                                     "C" = "True label C"))) +
  theme_bw() +
  theme(legend.position = "none") +
  scale_fill_brewer(palette = 'Accent', direction = -1) +
  labs(x = "",
       title = "Change in ECE with Calibration") +
  scale_x_discrete(labels = c("Uncalibrated", "Calibrated"))

ggsave('figure11.png', plot, width = 8, height = 4.5)
