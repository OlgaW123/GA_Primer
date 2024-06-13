library(ggplot2)
library(tidyverse)

pe_long <- pm %>%
  pivot_longer(cols = everything(), names_to = "Variable", values_to = "Value") %>%
  mutate(Index = rep(1:nrow(pm), each = ncol(pm)))

ggplot(pe_long, aes(x = Index, y = Value, color = Variable, group = Variable)) +
  geom_line(size = 1, lineend = "round") +  # adjust line
  labs(title = "Pm (mutation probability)", x = "Generation", y = "Best fitness score") +
  theme_minimal() +
  scale_color_discrete(name = "Pm") + 
  theme(
    plot.title = element_text(size = 15,hjust = 0.5),  
    axis.title.x = element_text(size = 14),             
    axis.title.y = element_text(size = 14),             
    axis.text.x = element_text(size = 10, hjust = 1),     
    axis.text.y = element_text(size = 10),                            
    legend.title = element_text(size = 14),             
    legend.text = element_text(size = 13)                              
  )
