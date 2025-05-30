
# Paper: Real-time coordination of parental provisioning revealed by high-resolution biologging in the wild #
## Script relative to the Result sub-section: 
### "Within-night adjustments of parental effort" only results about piece‐wise exponential additive mixed model (PAMM),
### "Contributions of specific parental behaviours to nestling survival and growth"
## some graphs and model summaries produced here are present in Supplementary Materials


# Packages ####
packages<-function(x){
  x<-as.character(match.call()[[2]])
  if (!require(x,character.only=TRUE)){
    install.packages(pkgs=x,repos="http://cran.r-project.org")
    require(x,character.only=TRUE)
  }
}


packages(tidyverse)
packages(tidylog)
packages(data.table)
packages(brms)
packages(parallel)
packages(coda)
packages(bayesplot)
packages(bayestestR)
packages(BayesFactor)
packages(loo)
packages(emmeans)
packages(sjPlot)
packages(pammtools)
packages(mgcv)
packages(performance)
packages(viridis)

# Load data ####

mypath <- "/Users/pbecciu/Desktop/Personal/CH/Lausanne/Work/ms/Partnership/code-tables_shared/"
subset.night.z <- read.csv(paste0(mypath, "nightly_params.csv")) # movement/foraging parameters relative to the Male and Female of the pair averaged by night
subset.broodID.z <- read.csv(paste0(mypath, "Pair_params.csv"))  # movement/foraging parameters relative to the Male and Female of the pair averaged by individual
cox.foragingMF <- read.csv(paste0(mypath, "pamm_table.csv"))
chicks_table <- read.csv(paste0(mypath, "chicks_table.csv"))


#_____________________________________________####
# Time-to-event model using PAMMs (pammtools) ####
# (Fig. 2 E, F)  ####
#_____________________________________________####

## subset by sex ####
db_gamF <- cox.foragingMF |> filter(Sex == "female")
db_gamM <- cox.foragingMF |> filter(Sex == "male")


## females ####
pedF <- db_gamF |>
  dplyr::select(BroodID_night, prey.nest.Fcoop01, time, status) |>
  mutate(BroodID_night = as.factor(BroodID_night)) |>
  as_ped(Surv(time, status)~ prey.nest.Fcoop01 + BroodID_night, zero = -0.1) |>
  na.omit()


pamm_F <- pammtools::pamm(
  ped_status ~ 
    s(tend) + 
    s(prey.nest.Fcoop01) + 
    ti(tend, prey.nest.Fcoop01) +  
    s(BroodID_night, bs = "re"),
  data = pedF,
  engine = "bam",
  method = "fREML",
  discretize = T)

summary(pamm_F)


pedF_df <- pedF |> 
  make_newdata(tend = unique(tend), prey.nest.Fcoop01 = seq_range(prey.nest.Fcoop01, by = 0.1), BroodID_night = first(BroodID_night)) |>
  filter(prey.nest.Fcoop01 < 0.6)

# plot "foraging" probability

pedF_df_probs <- pedF_df |> group_by(BroodID_night, prey.nest.Fcoop01) |> add_surv_prob(pamm_F)

plot_pammF <- ggplot(pedF_df_probs, 
                     aes(x = tend, y = surv_prob, ymax = surv_lower, ymin = surv_upper, group = prey.nest.Fcoop01)) +
  geom_line(aes(col = prey.nest.Fcoop01)) + 
  geom_ribbon(aes(fill = prey.nest.Fcoop01), alpha = 0.2) +
  scale_color_viridis_c() +
  scale_fill_viridis_c() +
  labs(title = "Females",
       y = expression(hat(F)(t)),
       x = paste0(expression(t), " ", "(hours from dusk)"),
       col = "FprovSHARE",
       fill = "FprovSHARE") 

# plot hazard

pedF_df_haz <- pedF_df |> add_hazard(pamm_F)

ggplot(pedF_df_haz, aes(x = tend, group = prey.nest.Fcoop01)) +
  geom_stephazard(aes(y = hazard, col = prey.nest.Fcoop01)) +
  geom_stepribbon(aes(ymin = ci_lower, ymax = ci_upper, fill = prey.nest.Fcoop01), alpha = 0.2) +
  scale_color_viridis_c() +
  scale_fill_viridis_c() +
  labs(title = "Females",
       y = expression(hat(Lambda)(t)),
       x = paste0(expression(t), " ", "(hours from dusk)"),
       col = "prop.bip",
       fill = "prop.bip")

# plot cumulative hazard

pedF_df_cumhaz <- pedF_df |> group_by(BroodID_night, prey.nest.Fcoop01) |> add_cumu_hazard(pamm_F)

ggplot(pedF_df_cumhaz, 
       aes(x = tend, y = cumu_hazard, ymin = cumu_lower, ymax = cumu_upper, group = prey.nest.Fcoop01)) +
  geom_hazard(aes(col = prey.nest.Fcoop01)) + 
  geom_ribbon(aes(fill = prey.nest.Fcoop01), alpha = 0.2) +
  scale_color_viridis_c() +
  scale_fill_viridis_c() +
  labs(title = "Females",
       y = expression(hat(Lambda)(t)),
       x = paste0(expression(t), " ", "(hours from dusk)"),
       col = "prop.bip",
       fill = "prop.bip")

## males ####
pedM <- db_gamM |>
  dplyr::select(BroodID_night, prey.nest.Fcoop01, time, status) |>
  mutate(BroodID_night = as.factor(BroodID_night)) |>
  as_ped(Surv(time, status)~ prey.nest.Fcoop01 + BroodID_night, zero = -0.1) |>
  na.omit()


pamm_M <- pamm(
  ped_status ~ 
    s(tend) + 
    s(prey.nest.Fcoop01) + 
    ti(tend, prey.nest.Fcoop01) +  
  s(BroodID_night, bs = "re"),
  data = pedM,
  engine = "bam",
  method = "fREML",
  discretize = T)


summary(pamm_M)

pedM_df <- pedM |> 
  make_newdata(tend = unique(tend), prey.nest.Fcoop01 = seq_range(prey.nest.Fcoop01, by = 0.1), BroodID_night = first(BroodID_night)) |>
  filter(prey.nest.Fcoop01 < 0.6)

# plot "foraging" probability

pedM_df_probs <- pedM_df |> group_by(BroodID_night, prey.nest.Fcoop01) |> add_surv_prob(pamm_M)

plot_pammM <- ggplot(pedM_df_probs, 
                     aes(x = tend, y = surv_prob, ymax = surv_lower, ymin = surv_upper, group = prey.nest.Fcoop01)) +
  geom_line(aes(col = prey.nest.Fcoop01)) + 
  geom_ribbon(aes(fill = prey.nest.Fcoop01), alpha = 0.2) +
  scale_color_viridis_c() +
  scale_fill_viridis_c() +
  labs(title = "Males",
       y = expression(hat(F)(t)),
       x = paste0(expression(t), " ", "(hours from dusk)"),
       col = "prop.bip",
       fill = "prop.bip") 
plot_pammM

# plot hazard
#Hazard Function, h(t): the instantaneous potential of experiencing an event at time t, conditional on having "survived" to that time

pedM_df_haz <- pedM_df |> add_hazard(pamm_M)

ggplot(pedM_df_haz, aes(x = tend, group = prey.nest.Fcoop01)) +
  geom_stephazard(aes(y = hazard, col = prey.nest.Fcoop01)) +
  geom_stepribbon(aes(ymin = ci_lower, ymax = ci_upper, fill = prey.nest.Fcoop01), alpha = 0.2) +
  scale_color_cvi_c("cvi_purples") +
  scale_fill_cvi_c("cvi_purples") +
  labs(title = "Males",
       y = expression(hat(Lambda)(t)),
       x = paste0(expression(t), " ", "(hours from dusk)"),
       col = "prop.bip",
       fill = "prop.bip") +
  THEME

# plot cumulative hazard
pedM_df_cumhaz <- pedM_df |> group_by(BroodID_night, prey.nest.Fcoop01) |> add_cumu_hazard(pamm_M)

ggplot(pedM_df_cumhaz, 
       aes(x = tend, y = cumu_hazard, ymin = cumu_lower, ymax = cumu_upper, group = prey.nest.Fcoop01)) +
  geom_hazard(aes(col = prey.nest.Fcoop01)) + 
  geom_ribbon(aes(fill = prey.nest.Fcoop01), alpha = 0.2) +
  scale_color_cvi_c("cvi_purples") +
  scale_fill_cvi_c("cvi_purples") +
  labs(title = "Males",
       y = expression(hat(Lambda)(t)),
       x = paste0(expression(t), " ", "(hours from dusk)"),
       col = "prop.bip",
       fill = "prop.bip") +
  THEME


panel.pammMF <- gridExtra::grid.arrange(plot_pammF, plot_pammM, nrow = 2)


# _________________________________####
# Nestling survival ####
# _________________________________####



## model loop ####

list_varcombination <- c(
  "prey.nest.Fcoop01_z + prey.per.chick_z + area_wildflower.zsqrt",
  "prey.nest.Fcoop01_z + area_wildflower.zsqrt",
  "prey.per.chick_z + area_wildflower.zsqrt",
  "prey.per.chickM_z + prey.per.chickF_z + area_wildflower.zsqrt",
  "hunt.att.sum.M_z + hunt.att.sum.F_z + area_wildflower.zsqrt",
  "prop.succ.dives.F_z + prop.succ.dives.M_z + area_wildflower.zsqrt",
  "vedba.avg.F_z + vedba.avg.M_z + area_wildflower.zsqrt"
)

Nsucc.model_list <- list()
for (i in 1:length(list_varcombination)) {
  my.cores <- detectCores()-1
  
  Nsucc.m <- brms::brm(
    paste("Nb_Nestlings2 | trials(Nb_Eggs) ~", list_varcombination[i]," + yearF + (1|BroodID)"),
    chains = 3, cores = my.cores, warmup = 5000, iter = 20000, 
    family = binomial("logit"), 
    control = list(adapt_delta = 0.9),
    save_pars = save_pars(all = TRUE),
    data = subset.broodID.z)
  
  print(summary(Nsucc.m))
  
  
  Nsucc.model_list[[i]] <- Nsucc.m
}

## check models ####
pp_check(Nsucc.model_list[[1]], ndraws = 200)
pp_check(Nsucc.model_list[[2]], ndraws = 200)
pp_check(Nsucc.model_list[[3]], ndraws = 200)
pp_check(Nsucc.model_list[[4]], ndraws = 200)
pp_check(Nsucc.model_list[[5]], ndraws = 200)
pp_check(Nsucc.model_list[[6]], ndraws = 200)
pp_check(Nsucc.model_list[[7]], ndraws = 200)

plot(Nsucc.model_list[[1]])


##  plots ####
conditional_effects(Nsucc.model_list[[4]], ndraws = 1000)
conditional_effects(Nsucc.model_list[[1]], ndraws = 1000)
conditional_effects(Nsucc.model_list[[5]], ndraws = 1000)
conditional_effects(Nsucc.model_list[[6]], ndraws = 1000)
conditional_effects(Nsucc.model_list[[7]], ndraws = 1000)


plot_models(
  Nsucc.model_list[[1]], 
  Nsucc.model_list[[2]],
  Nsucc.model_list[[3]],
  Nsucc.model_list[[4]], 
  Nsucc.model_list[[5]], 
  Nsucc.model_list[[6]], 
  Nsucc.model_list[[7]], 
  rm.terms = "b_Intercept") +
  labs(subtitle = "Nestling survival") +
  scale_y_log10(limits = c(0.5, 3)) +
  geom_hline(yintercept = 1, linetype = "dashed", alpha = 0.5) +
  THEME0 + 
  theme(legend.position = "none")

## summary tables ####

for (i in 1:length(Nsucc.model_list)) {
  
  table <- describe_posterior(
    Nsucc.model_list[[i]],
    effects = "all",
    component = "all",
    centrality = "all",
    test = c("p_direction", "p_significance"),
    verbose = FALSE
  )
  
  print(table)
}

# _________________________________####
# Chick body condition models (Fig. 3)####


## _____________####
## weight ####
## _____________####

## subset ####

chicks_table1 <- chicks_table |>
  filter(!is.na(time.group)) |>
  droplevels() |>
  ungroup()


### models ####

my.cores <- detectCores()-2

bform <- bf(CalculatedMass ~ 0 + Intercept +
              prey.nest.Fcoop01_z + # biparentality
              time.group + 
              Rank +
              time.group:prey.nest.Fcoop01_z +
              Rank:prey.nest.Fcoop01_z +
              time.group:prey.nest.Fcoop01_z:Rank +
              (1|p|BroodID) + (1|q|RingId),
            sigma ~ prey.nest.Fcoop01_z + time.group + Rank)


Mchick.weight.ens <- brms::brm(bform,
                               chains = 3, cores = my.cores, warmup = 5000, iter = 20000,
                               data = chicks_table1)



### model check ####

print(Mchick.weight.ens, digits = 3)
pp_check(Mchick.weight.ens, ndraws = 100)
mcmc_plot(Mchick.weight.ens, type = "acf") #correlation plots
plot(Mchick.weight.ens)

### summary ####
describe_posterior(
  Mchick.weight.ens,
  effects = "all",
  component = "all",
  centrality = "all",
  test = c("p_direction", "p_significance"),
  verbose = FALSE
)

### contrasts ####

emmip(Mchick.weight.ens, time.group~prey.nest.Fcoop01_z|Rank, cov.reduce = range, CIs = TRUE)

emt <- emtrends(Mchick.weight.ens, c("time.group", "Rank"), var = "prey.nest.Fcoop01_z")
emt
pairs(emt)


### plots ####

grid.weight.ens = chicks_table1 |>
  group_by(time.group, Rank) |>
  data_grid(prey.nest.Fcoop01_z, BroodID=first(BroodID), RingId=first(RingId))

gc() ; means.weight.ens = grid.weight.ens |>
  add_epred_draws(Mchick.weight.ens, ndraws = 1000)

plot_weightrankTime <- ggplot(chicks_table1, aes(x = prey.nest.Fcoop01_z, y = CalculatedMass, color = ordered(Rank))) +
  facet_wrap(.~time.group) +
  stat_lineribbon(data = means.weight.ens, aes(y = .epred), .width = c(.95, .80, .50), alpha = 1/4) +
  scale_fill_brewer(palette = "Greys") +
  scale_color_brewer(palette = "Dark2") +
  THEME0
plot_weightrankTime



## _____________####
## wing growth ####
## _____________####


### model ####

my.cores <- detectCores()-2

bform <- bf(LeftWing ~ 0 + Intercept +
              prey.nest.Fcoop01_z +
              time.group + 
              Rank +
              time.group:prey.nest.Fcoop01_z +
              Rank:prey.nest.Fcoop01_z +
              time.group:prey.nest.Fcoop01_z:Rank +
              (1|p|BroodID) + (1|q|RingId),
            sigma ~ prey.nest.Fcoop01_z + time.group + Rank)


Mchick.wing.ens <- brms::brm(bform,
                             chains = 3, cores = my.cores, warmup = 5000, iter = 20000, 
                             data = chicks_table1)

### model check ####

print(Mchick.wing.ens, digits = 3)
pp_check(Mchick.wing.ens, ndraws = 100)
mcmc_plot(Mchick.wing.ens, type = "acf") #correlation plots
plot(Mchick.wing.ens)

### summary ####
describe_posterior(
  Mchick.wing.ens,
  effects = "all",
  component = "all",
  centrality = "all",
  test = c("p_direction", "p_significance"),
  verbose = FALSE
)

### contrasts ####

emmip(Mchick.wing.ens, time.group~prey.nest.Fcoop01_z|Rank, cov.reduce = range, CIs = TRUE)

emt <- emtrends(Mchick.wing.ens, c("time.group", "Rank"), var = "prey.nest.Fcoop01_z")
emt
pairs(emt)

### plots ####

grid.wing.ens = chicks_table1 |>
  group_by(time.group, Rank) |>
  data_grid(prey.nest.Fcoop01_z, BroodID=first(BroodID), RingId=first(RingId))

gc() ; means.wing.ens = grid.wing.ens |>
  add_epred_draws(Mchick.wing.ens, ndraws = 1000)

plot_wingrankTime <- ggplot(chicks_table1, aes(x = prey.nest.Fcoop01_z, y = LeftWing, color = ordered(Rank))) +
  facet_wrap(.~time.group) +
  stat_lineribbon(data = means.wing.ens, aes(y = .epred), .width = c(.95, .80, .50), alpha = 1/4) +
  # geom_point(data = test.ens) +
  scale_fill_brewer(palette = "Greys") +
  scale_color_brewer(palette = "Dark2") +
  THEME0
plot_wingrankTime


## _____________####
## daily weight change ####
## _____________####


## subset ####

chicks_table2 <- chicks_table |>
  filter(!is.na(time.group)) |>
  filter(!is.na(daily.weight.change)) |>
  filter(time.group == "2" | time.group == "3") |>
  droplevels() |>
  ungroup()



### model ####

my.cores <- detectCores()-2

bform <- bf(daily.weight.change ~ 0 + Intercept +
              prey.nest.Fcoop01_z +
              time.group + 
              Rank +
              time.group:prey.nest.Fcoop01_z +
              Rank:prey.nest.Fcoop01_z +
              time.group:prey.nest.Fcoop01_z:Rank +
              (1|p|BroodID) + (1|q|RingId),
            sigma ~ prey.nest.Fcoop01_z + time.group + Rank)


Mchick.weightC.ens <- brms::brm(bform,
                                chains = 3, cores = my.cores, warmup = 5000, iter = 20000, 
                                data = chicks_table2)

### model check ####

print(Mchick.weightC.ens, digits = 3)
pp_check(Mchick.weightC.ens, ndraws = 100)
mcmc_plot(Mchick.weightC.ens, type = "acf")
plot(Mchick.weightC.ens)

### summary ####
describe_posterior(
  Mchick.weightC.ens,
  effects = "all",
  component = "all",
  centrality = "all",
  test = c("p_direction", "p_significance"),
  verbose = FALSE
)

### contrasts ####

emmip(Mchick.weightC.ens, time.group~prey.nest.Fcoop01_z|Rank, cov.reduce = range, CIs = TRUE)

emt <- emtrends(Mchick.weightC.ens, c("time.group", "Rank"), var = "prey.nest.Fcoop01_z")
emt
pairs(emt)

### plots ####

grid.weightC.ens = chicks_table2 |>
  group_by(time.group, Rank) |>
  data_grid(prey.nest.Fcoop01_z, BroodID=first(BroodID), RingId=first(RingId))

gc() ; means.weightC.ens = grid.weightC.ens |>
  add_epred_draws(Mchick.weightC.ens, ndraws = 1000)

plot_weightCrankTime <- ggplot(chicks_table2, aes(x = prey.nest.Fcoop01_z, y = daily.weight.change, color = ordered(Rank))) +
  facet_wrap(.~time.group) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.7) +
  stat_lineribbon(data = means.weightC.ens, aes(y = .epred), .width = c(.95, .80, .50), alpha = 1/4) +
  scale_fill_brewer(palette = "Greys") +
  scale_color_brewer(palette = "Dark2") +
  THEME0
plot_weightCrankTime


## _____________####
## daily wing growth change ####
## _____________####

### models ####

my.cores <- detectCores()-2

bform <- bf(daily.wing.growth ~ 0 + Intercept +
              prey.nest.Fcoop01_z +
              time.group + 
              Rank +
              time.group:prey.nest.Fcoop01_z +
              Rank:prey.nest.Fcoop01_z +
              time.group:prey.nest.Fcoop01_z:Rank +
              (1|p|BroodID) + (1|q|RingId),
            sigma ~ prey.nest.Fcoop01_z + time.group + Rank)


Mchick.wingG.ens <- brms::brm(bform,
                              chains = 3, cores = my.cores, warmup = 5000, iter = 20000, 
                              data = chicks_table2)


### model check ####

print(Mchick.wingG.ens, digits = 3)
pp_check(Mchick.wingG.ens, ndraws = 100)
mcmc_plot(Mchick.wingG.ens, type = "acf") #correlation plots
plot(Mchick.wingG.ens)

### summary ####
describe_posterior(
  Mchick.wingG.ens,
  effects = "all",
  component = "all",
  centrality = "all",
  test = c("p_direction", "p_significance"),
  verbose = FALSE
)


### contrasts ####

emmip(Mchick.wingG.ens, time.group~prey.nest.Fcoop01_z|Rank, cov.reduce = range, CIs = TRUE)

emt <- emtrends(Mchick.wingG.ens, c("time.group", "Rank"), var = "prey.nest.Fcoop01_z")
emt
pairs(emt)

### plots ####

grid.wingG.ens = chicks_table2 |>
  group_by(time.group, Rank) |>
  data_grid(prey.nest.Fcoop01_z, BroodID=first(BroodID), RingId=first(RingId))

gc() ; means.wingG.ens = grid.wingG.ens |>
  add_epred_draws(Mchick.wingG.ens, ndraws = 1000)

plot_wingGrankTime <- ggplot(chicks_table2, aes(x = prey.nest.Fcoop01_z, y = daily.wing.growth, color = ordered(Rank))) +
  facet_wrap(.~time.group) +
  stat_lineribbon(data = means.wingG.ens, aes(y = .epred), .width = c(.95, .80, .50), alpha = 1/4) +
  scale_fill_brewer(palette = "Greys") +
  scale_color_brewer(palette = "Dark2") +
  THEME0
plot_wingGrankTime


