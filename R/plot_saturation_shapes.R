# Script to plot the shapes of tau(X) against x1 for saturation scenarios
# Based on the generate_data_medical function in export_import.R

library(ggplot2)
library(dplyr)
library(tidyr)

# 1. Generate data
set.seed(1848)
n <- 10000  # Large n to get a smooth density and line
x1 <- rnorm(n, mean=0, sd=1)

# 2. Define the saturation functions
# The "medical_saturation" scenario from generate_data_medical:
# Patients with low x1 get a small baseline benefit (0.5).
# As x1 increases, benefit grows rapidly but caps at 4.5 (saturation).
tau_medical <- 0.5 + 4 / (1 + exp(-2.5 * x1))

# Defining a hypothetical "mild_saturation" scenario
# Assuming "mild" means either a less steep growth or a lower maximum cap.
# Here we use a less steep slope (-1.0 instead of -2.5) as an example.
# Feel free to adjust the parameters for the "mild" scenario!
tau_mild <- 0.5 + 4 / (1 + exp(-1.0 * x1))

# Create a data frame for plotting
df <- data.frame(
  x1 = x1,
  Medical = tau_medical,
  Mild = tau_mild
)

# 3. Plot the shapes of tau(X) against x1
df_long <- df %>%
  pivot_longer(cols = c(Medical, Mild), names_to = "Scenario", values_to = "tau")

p_shapes <- ggplot(df_long, aes(x = x1, y = tau, color = Scenario)) +
  geom_line(linewidth = 1.2) +
  theme_bw(base_size = 14) +
  labs(
    title = "Shape of Treatment Effect tau(X) vs x1",
    subtitle = "Comparing Medical and Mild Saturation Scenarios",
    x = "x1 (Biomarker, Standard Normal)",
    y = "tau(X)",
    color = "Scenario"
  ) +
  scale_color_manual(values = c("Medical" = "#D55E00", "Mild" = "#0072B2")) +
  theme(legend.position = "bottom")

# 4. Plot the density support for x1
p_density <- ggplot(df, aes(x = x1)) +
  geom_density(fill = "gray50", alpha = 0.5, color = "black", linewidth = 1) +
  theme_bw(base_size = 14) +
  labs(
    title = "Density Support for x1",
    subtitle = "x1 ~ N(0, 1)",
    x = "x1",
    y = "Density"
  )

# Print the plots
print(p_shapes)
print(p_density)

# Optional: You can save them using ggsave
# ggsave("tau_shapes.png", p_shapes, width = 8, height = 6)
# ggsave("x1_density.png", p_density, width = 8, height = 6)
