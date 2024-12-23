# Calculate standardised effect sizes 
# see also https://jslefche.github.io/piecewiseSEM/reference/coefs.html
# and https://doi.org/10.1002/ecs2.2283
# we will be using the latent-theoretic approach for all non-Gaussian responses

library(tidybayes)

source("code/model.R")

# Latent-theoretic SDs
linpred_post <- posterior_linpred(fit_veg)
dpar_sigma_post <- posterior_linpred(fit_veg, resp = "Litter", dpar = "sigma")
dpar_phi_post <- posterior_linpred(fit_veg, resp = "Densio", dpar = "phi")
LT_sd <- colMeans(apply(linpred_post, c(1, 3), sd))
LT_sd["Clid"] <- sqrt(LT_sd["Clid"]^2 + pi^2/3)
LT_sd <- c(LT_sd, 
           sigma = sd(colMeans(dpar_sigma_post)),
           phi = sd(colMeans(dpar_phi_post)))

# coefficients and standardisations
effects_mcmc <-
  gather_draws(
    fit_veg,
    `b_.*|bsp_.*`, 
    regex = TRUE) %>% 
  filter(!str_detect(.variable, "_Intercept")) %>% 
  mutate(.variable = str_remove(.variable, "b_|bsp_"),
         .variable = str_replace(.variable, "phi_Densio", "phi"),
         .variable = str_replace(.variable, "sigma_Litter", "sigma"),
         .variable = str_replace(.variable, "miDensio", "Densio"),
         .variable = str_replace(.variable, "dist_veg", "dist")) %>% 
  separate(.variable, c("Response", "Predictor"), sep = "_") %>% 
  # scale by predictors
  mutate_at(vars(.value),
            ~ifelse(Predictor == "Litter",
                   . * LT_sd["Litter"],
                   ifelse(Predictor == "Densio",
                          . * LT_sd["Densio"], 
                          .))) %>% 
  # then scale by responses
  mutate_at(vars(.value),
            ~ifelse(Response %in% c("Litter", "sigma"),
                    . / LT_sd["Litter"],
                    ifelse(Response %in% c("Densio", "phi"),
                           . / LT_sd["Densio"], 
                           ifelse(Response == "Clid",
                                  . / LT_sd["Clid"], 
                                  .)))) 
effects <- 
  effects_mcmc %>% 
  # summarise
  group_by(Response, Predictor) %>% 
  median_qi(.value, 
            .width = 0.89) %>% 
  mutate(Response = fct_relevel(Response, 
                                "Clid", "Densio", "Litter"),
         Predictor = fct_relevel(Predictor, 
                                "Densio", "Litter", "dist"))


# total effect of distance on clidemia
total_effect_Clid <- 
  effects_mcmc %>% 
  pivot_wider(names_from = c(Response, Predictor),
              values_from = .value) %>% 
  mutate(
    Direct_Dist = Clid_dist,
    Indirect_Dist.Canopy = Clid_Densio * Densio_dist,
    Indirect_Dist.Litter = Clid_Litter * Litter_dist,
    Total_Dist = Direct_Dist + Indirect_Dist.Canopy + Indirect_Dist.Litter,
  ) %>% 
  select(.chain, .iteration, .draw, 
         starts_with("Direct_"),
         starts_with("Indirect_"),
         starts_with("Total_")) %>% 
  pivot_longer(cols = c(starts_with("Direct_"),
                        starts_with("Indirect_"),
                        starts_with("Total_")),
               names_to = ".variable",
               values_to = ".value") %>% 
  # summarise
  group_by(.variable) %>% 
  median_qi(.value,
            .width = 0.89) %>% 
  separate(.variable, c("Type", "Path"), sep = "_") %>% 
  mutate(Path = str_replace(Path, "\\.", paste0(" ", sprintf('\u2799'), " ")),
         Path = paste(Path, sprintf('\u2799'), "Y"),
         Path = fct_rev(fct_relevel(as.factor(Path), "Dist ➙ Y")))


ggplot(total_effect_Clid) +
  facet_grid(rows = vars(Type), 
             scales = "free_y",
             space = "free_y",
             switch = "y") +
  geom_vline(xintercept = 0, 
             linewidth = 0.3,
             colour = "grey") +
  geom_pointinterval(aes(.value, Path, xmin = .lower, xmax = .upper,
                         fill = Type),
                     colour = "black",
                     pch = 21,
                     linewidth = 1,
                     point_size = 3,
                     position = ggstance::position_dodgev(height=0.5),
                     height = 0) +
  scale_fill_manual(values = c("darkgrey", "white", "black"),
                    guide = "none") +
  labs(x = "Standardised effect size") +
  theme_classic() +
  theme(strip.background = element_rect(fill = "grey", colour = "grey"), 
        strip.placement = "outside",
        strip.text = element_text(size = 11),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 10, colour = "black"),
        axis.text.x = element_text(colour = "black"),
        axis.ticks.length = unit(-0.1, "cm"),
        plot.margin = margin(20, 10, 0, 0))




# R2 ----------------------------------------------------------------------

bayes_R2(fit_veg, resp = "Clid")
