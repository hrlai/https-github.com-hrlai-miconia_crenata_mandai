library(modelr)
library(tidybayes)
library(scales)
library(ggpubr)

logit <- function(x) log(x / (1 - x))
ilogit <- function(x) 1 / (1 + exp(-x))

n_new <- 100
n_draws <- 1000

pred_direct <- 
  clidemia_s %>%
  data_grid(dist_veg = seq_range(dist_veg, n = n_new),
            Litter = mean(Litter), 
            Densio = mean(Densio, na.rm = TRUE),
            Easting = mean(Easting),
            Northing = mean(Northing)) %>%
  add_linpred_draws(fit_veg, 
                    resp = c("Clid"),
                    newdata = .,
                    transform = TRUE,
                    ndraws = n_draws) %>% 
  group_by(dist_veg) %>% 
  median_qi(.linpred, .width = 0.89)

pred_total <- 
  clidemia_s %>%
  data_grid(dist_veg = seq_range(dist_veg, n = n_new),
            Easting = mean(Easting),
            Northing = mean(Northing)) %>%
  add_predicted_draws(fit_veg, 
                      resp = c("Densio", "Litter"),
                      newdata = .,
                      ndraws = n_draws) %>% 
  ungroup %>% 
  pivot_wider(names_from = .category,
              values_from = .prediction) %>% 
  select(dist_veg, Densio, Litter, Easting, Northing) %>% 
  add_linpred_draws(fit_veg, 
                  resp = c("Clid"),
                  newdata = .,
                  transform = TRUE,
                  ndraws = n_draws) %>% 
  group_by(dist_veg) %>% 
  median_qi(.linpred, .width = 0.89)

pred <- 
  bind_rows(
    pred_direct %>% mutate(Effect = "Direct"),
    pred_total  %>% mutate(Effect = "Total")
  ) %>% 
  mutate(dist_veg = dist_veg * sd(clidemia$dist_veg) + mean(clidemia$dist_veg))

p_counterfactual <- 
  ggplot(pred) +
  geom_line(aes(dist_veg, .linpred, colour = Effect),
              linewidth = 1) +
  geom_line(aes(dist_veg, .lower, colour = Effect),
              linewidth = 0.5, 
              lty = 2) +
  geom_line(aes(dist_veg, .upper, colour = Effect),
              linewidth = 0.5, 
              lty = 2) +
  scale_colour_grey(start = 0.75, end = 0) +
  scale_y_continuous(breaks = seq(0, 1, 0.2)) +
  labs(x = "Distance from edge [m]",
       y = expression(paste("Pr(", italic(Miconia), " occurrence = 1)"))) +
  coord_cartesian(expand = FALSE, ylim = c(0, 1)) +
  theme_classic() +
  theme(axis.ticks.length = unit(-0.1, "cm"),
        axis.text = element_text(colour = "black"),
        legend.position = c(0.85, 0.8),
        plot.margin = margin(20, 1, 1, 1))
  

ggarrange(p_effects, p_counterfactual,
          nrow = 1,
          widths = c(1, 1),
          labels = LETTERS[1:2])
