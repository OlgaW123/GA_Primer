install.packages("ggplot2")

library(ggplot2)

data <- data.frame(
  max_gen = c(10, 20, 50, 100, 200),
  best_fitness_score = c(0.013422819, 0.034482759, 0.042105263, 0.042105263, 0.0625)
)

ggplot(data, aes(x = max_gen, y = best_fitness_score)) +
  geom_point(size = 2, color = 'blue') +
  labs(title = "Fitness Score vs. Number of Generations", x = "Number of Generations", y = "Best fitness score") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 15,hjust = 0.5),  
    axis.title.x = element_text(size = 14),            
    axis.title.y = element_text(size = 14),             
    axis.text.x = element_text(size = 10, hjust = 1),      
    axis.text.y = element_text(size = 10),                                                       
  )
