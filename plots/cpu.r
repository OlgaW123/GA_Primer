library(ggplot2)
library(gridExtra)

data <- read.csv("C:/Users/wiero/OneDrive/Pulpit/BLAT/wykresy/cpu/resource_usage_2.txt", skip = 1, header = FALSE)

colnames(data) <- c("Time_s", "CPU", "Memory")

head(data)

# Create a plot for CPU usage
cpu_plot <- ggplot(data, aes(x = Time_s, y = CPU)) +
  geom_line(color = "red") +
  labs(title = "CPU Usage Over Time",
       x = "Time (seconds)",
       y = "CPU Usage (%)") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 15,hjust = 0.5),  # Title size, style, and centering
    axis.title.x = element_text(size = 14),             # X-axis title size and style
    axis.title.y = element_text(size = 14),             # Y-axis title size and style
    axis.text.x = element_text(size = 10, hjust = 1),      # X-axis text size and angle
    axis.text.y = element_text(size = 10),                             # Y-axis text size                           
  )

# Create a plot for Memory usage
memory_plot <- ggplot(data, aes(x = Time_s, y = Memory)) +
  geom_line(color = "blue") +
  labs(title = "Memory Usage Over Time",
       x = "Time (seconds)",
       y = "Memory Usage (MB)") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 15,hjust = 0.5),  
    axis.title.x = element_text(size = 14),             
    axis.title.y = element_text(size = 14),             
    axis.text.x = element_text(size = 10, hjust = 1),      
    axis.text.y = element_text(size = 10),                                                      
  )

combined_plot <- grid.arrange(cpu_plot, memory_plot, ncol = 1)
