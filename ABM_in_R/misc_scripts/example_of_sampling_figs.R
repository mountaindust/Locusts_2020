library(tidyverse)
library(grid)
# test data for sampling plot ------------------------------------------------------

theme_set(MCMsBasics::minimal_ggplot_theme())

# this is just making a plot to show the "sampling" procedure we use, basically showing how we're getting a 1m cross section of the band
data <- data.frame(x = 10-rexp(10000), y = runif(10000, min = 0, max = 10))
data %>% 
  filter(x > 1) %>% 
  filter(!(y < 0.5 & x < 4)) %>% 
  ggplot(aes(x = x, y = y)) + 
  geom_jitter(width = 0.8, height = 0.8, alpha = 0.3) +
  geom_rect(xmin = 0.5, xmax = 11, ymin = 3, ymax = 4, color = "red", fill = NA) +
  xlim(0,11) +
  MCMsBasics::minimal_ggplot_theme() +
  geom_segment(data = data.frame(x = 1, y = 0), xend = 4, yend = 0, color = "red", arrow = arrow(type = "closed", length = unit(0.25, "cm"))) +
  annotate("text", label = "direction of movement", x = 2.4, y = -0.25, size = 3, colour = "red") +
  xlab("band width (m)") +
  ylab("band length (m)") +
  annotate("text", label = "1 meter", angle = 90, x = 0.2, y = 3.5, size = 3, colour = "red")
ggsave("ABM_in_R/figures/example_of_sampling.jpg", width = 4, height = 8)

# trying a little fancier
g <- rasterGrob(t(c("tan", "tan", "forestgreen")), width=unit(1,"npc"), height = unit(1,"npc"), 
                interpolate = TRUE) 
data %>% 
  filter(x > 1) %>% 
  filter(!(y < 0.5 & x < 4)) %>% 
  ggplot(aes(x = x, y = y)) + 
  annotation_custom(g, xmin=0, xmax=13, ymin=-Inf, ymax=Inf) +
  geom_jitter(width = 0.8, height = 0.8, alpha = 0.3) +
  #annotate("rect", xmin = 0, xmax = 13, ymin = 3, ymax = 4, color = NA, fill = "grey80", alpha = 0.4) +
  annotate("segment", x = 0, y = 3, xend = 13, yend = 3, color = "red") +
  annotate("segment", x = 0, y = 4, xend = 13, yend = 4, color = "red") +
  xlim(0,13) +
  MCMsBasics::minimal_ggplot_theme() +
  geom_segment(data = data.frame(x = 1, y = 0), xend = 4, yend = 0, color = "black", arrow = arrow(type = "closed", length = unit(0.25, "cm"))) +
  annotate("text", label = "direction of movement", x = 2.4, y = -0.25, size = 3, colour = "black") +
  xlab("band width (m)") +
  ylab("band length (m)") +
  annotate("text", label = "1 meter", angle = 90, x = 0.2, y = 3.5, size = 3, colour = "black")
ggsave("ABM_in_R/figures/example_of_sampling_gradient.jpg", width = 4, height = 8)

# one more version
data %>% 
  filter(x > 1) %>% 
  filter(!(y < 0.5 & x < 4)) %>% 
  ggplot(aes(x = x, y = y)) + 
  annotation_custom(g, xmin=0, xmax=13, ymin=-Inf, ymax=Inf) +
  geom_jitter(width = 0.8, height = 0.8, alpha = 0.3) +
  annotate("segment", x = 0, y = 3, xend = 13, yend = 3, color = "red") +
  annotate("segment", x = 0, y = 4, xend = 13, yend = 4, color = "red") +
  xlim(0,13) +
  MCMsBasics::minimal_ggplot_theme() +
  geom_segment(data = data.frame(x = 1, y = 0), xend = 4, yend = 0, color = "red", arrow = arrow(type = "closed", length = unit(0.25, "cm"))) +
  annotate("text", label = "direction of movement", x = 2.4, y = -0.25, size = 3, colour = "red") +
  xlab("band width (m)") +
  ylab("band length (m)") +
  annotate("text", label = "1 meter", angle = 90, x = 0.2, y = 3.5, size = 3, colour = "red")
ggsave("ABM_in_R/figures/example_of_sampling_gradient_2.jpg", width = 4, height = 8)

# more zoomed out

g <- rasterGrob(t(c("tan", "tan", "tan", "forestgreen", "forestgreen")), width=unit(1,"npc"), height = unit(1,"npc"), interpolate = TRUE) 

data %>% 
  filter(x > 1) %>% 
  filter(!(y < 0.5 & x < 4)) %>% 
  ggplot(aes(x = x, y = y)) + 
  annotation_custom(g, xmin=0, xmax=15.5, ymin=-0.8, ymax=10.8) +
  geom_jitter(width = 0.8, height = 0.8, alpha = 0.3) +
  annotate("rect", xmin = 0, xmax = 15.5, ymin = 3, ymax = 4, color = NA, fill = "grey80", alpha = 0.4) +
  xlim(0,15.5) +
  geom_segment(data = data.frame(x = 1, y = 0), xend = 4, yend = 0, color = "black", arrow = arrow(type = "closed", length = unit(0.25, "cm"))) +
  annotate("text", label = "direction of\nmovement", x = 2.45, y = -0.35, size = 3.5, colour = "black",
           family= theme_get()$text[["family"]]) +
  xlab("band width") +
  ylab("band length") +
  annotate("text", label = "1 meter", angle = 90, x = 0.3, y = 3.5, size = 3.5, colour = "black",
           family= theme_get()$text[["family"]]) +
  theme(axis.text = element_blank())
ggsave("ABM_in_R/figures/example_of_sampling_gradient_3.jpg", width = 4, height = 8)

# same thing but a wider version instead of taller

text_ann <- textGrob("1 meter", gp=gpar(fontsize=8, family = theme_get()$text[["family"]]), rot = 270)

data %>% 
  filter(x > 1) %>% 
  filter(!(y < 0.5 & x < 4)) %>% 
  ggplot(aes(x = x, y = y)) + 
  annotation_custom(g, xmin=0, xmax=15.5, ymin=-0.8, ymax=10.8) +
  geom_jitter(width = 0.8, height = 0.8, alpha = 0.3) +
  annotate("rect", xmin = 0, xmax = 15.5, ymin = 3, ymax = 4, color = NA, fill = "grey80", alpha = 0.4) +
  xlim(0,15.5) +
  geom_segment(data = data.frame(x = 12, y = 1), xend = 14.5, yend = 1, color = "black", arrow = arrow(type = "closed", length = unit(0.25, "cm"))) +
  annotate("text", label = "direction of\nmovement", x = 13.2, y = 0.25, size = 3.5, colour = "black", family= theme_get()$text[["family"]]) +
  xlab("band width") +
  ylab("band length") +
  annotation_custom(text_ann, xmin = 15.7, xmax = 15.7, ymin = 3.5, ymax = 3.5) +
  theme(axis.text = element_blank()) +
  coord_cartesian(clip = 'off')
ggsave("ABM_in_R/figures/example_of_sampling_gradient_wide.jpg", width = 8, height = 4)


