
# Paper: DReal-time coordination of parental provisioning revealed by high-resolution biologging in the wild #
## Script relative to the Result sub-sections: 
### "Variation in biparental provisioning", 
### "Within-night adjustments of parental effort" and 
### "Between-night carry-over effects"
## some graphs and model outputs produced here are present in Supplementary Materials


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
packages(viridis)
packages(rptR)

# Load data ####

mypath <- "/Users/pbecciu/Desktop/Personal/CH/Lausanne/Work/ms/Partnership/code-tables_shared/"
# mypath <- ""
subset.night.z <- read.csv(paste0(mypath, "nightly_params.csv")) # movement/foraging parameters relative to the Male and Female of the pair averaged by night
subset.broodID.z <- read.csv(paste0(mypath, "Pair_params.csv"))  # movement/foraging parameters relative to the Male and Female of the pair averaged by individual

## subsetting data for some models####

dat.mod <- subset.night.z |>
  dplyr::select(BroodID, night_nb_real, prey.to.nest.F10, preys.nest10, prey.nest.Fcoop01_z, prey.per.chick, prey.per.chickF, prey.per.chickM,
                vedba.avg.F_z,
                hunt.att.sum.F_z, 
                prop.succ.dives.F_z, 
                prop.prey.eaten.prey.capturedF1_z,
                vedba.avg.M_z,
                hunt.att.sum.M_z,
                prop.succ.dives.M_z,
                prop.prey.eaten.prey.capturedM_z,
                area_wildflower.zsqrt,
                yearF,
                n_nest_encounters_z,
                n_out_encounters_z,
                time_at_nest_M_z, 
                time_at_nest_F_z,
                SMI_M_z, 
                SMI_F_z,
                SMI_diff_F_z,
                SMI_diff_M_z,
                brood.size.change01F,
                youngest.chick.age_z) |>
  mutate(across(
    where(function(x) (is.matrix(x) || is.array(x)) && is.numeric(as.matrix(x))),
    ~ as.numeric(as.matrix(.))
  )) |>
  as.data.frame()

## prep lagged data (carry-over effect models)####

dat.lag <- dat.mod |>
  arrange(BroodID, night_nb_real) |>
  group_by(BroodID) |>
  mutate(
    prey.nest.Fcoop01_z_prev = lag(prey.nest.Fcoop01_z, order_by = night_nb_real),
    prey.per.chick_z = as.single(scale(prey.per.chick)),
    prey.per.chick_z_prev = lag(prey.per.chick_z, order_by = night_nb_real),
    prey.per.chickF_z = as.single(scale(prey.per.chickF)),
    prey.per.chickM_z = as.single(scale(prey.per.chickM)),
    prey.per.chickF_z_prev = lag(prey.per.chickF_z, order_by = night_nb_real),
    prey.per.chickM_z_prev = lag(prey.per.chickM_z, order_by = night_nb_real),
    vedba.avg.F_z_prev = lag(vedba.avg.F_z, order_by = night_nb_real),
    hunt.att.sum.F_z_prev = lag(hunt.att.sum.F_z, order_by = night_nb_real),
    prop.succ.dives.F_z_prev = lag(prop.succ.dives.F_z, order_by = night_nb_real),
    prop.prey.eaten.prey.capturedF1_z_prev = lag(prop.prey.eaten.prey.capturedF1_z, order_by = night_nb_real),
    vedba.avg.M_z_prev = lag(vedba.avg.M_z, order_by = night_nb_real),
    hunt.att.sum.M_z_prev = lag(hunt.att.sum.M_z, order_by = night_nb_real),
    prop.succ.dives.M_z_prev = lag(prop.succ.dives.M_z, order_by = night_nb_real),
    prop.prey.eaten.prey.capturedM_z_prev = lag(prop.prey.eaten.prey.capturedM_z, order_by = night_nb_real),
    n_nest_encounters_z_prev = lag(n_nest_encounters_z, order_by = night_nb_real),
    n_out_encounters_z_prev = lag(n_out_encounters_z, order_by = night_nb_real),
    time_at_nest_M_z_prev = lag(time_at_nest_M_z, order_by = night_nb_real),
    time_at_nest_F_z_prev = lag(time_at_nest_F_z, order_by = night_nb_real)
  ) |>
  ungroup() |>
  na.omit() 

# ______________####
# Repeatability ####
# ______________####

rep.data <- subset.night.z |>
  dplyr::select(prey.to.nest.F10, prey.to.nest.M10, brood.size.fin, yearF, BroodID) |>
  mutate(brood.size.fin = as.factor(brood.size.fin),
         yearF = as.factor(yearF))

rep.biparentality <- rpt(cbind(prey.to.nest.F10, prey.to.nest.M10) ~ brood.size.fin + (1| BroodID), 
                 grname = c("BroodID"), data = subset.night.z, 
                 datatype = "Proportion", link = "logit", nboot = 1000, npermut = 0, ratio = T)
print(rep.biparentality)
plot(rep.biparentality)

summary(rep.biparentality)

# __________________________________________________####
# Relationship between BROOD SIZE and BIPARENTALITY ####
# __________________________________________________####

brood.mod.data <- subset.night.z |>
  dplyr::select(prey.to.nest.F10, preys.nest10, brood.size.fin, BroodID) |>
  mutate(brood.size.fin = as.factor(brood.size.fin))

## model ####

my.cores <- detectCores()-1

brood_sizeM <- brms::brm(prey.to.nest.F10 | trials(preys.nest10) ~
                                 0 + brood.size.fin +
                                 (1|BroodID),
                               chains = 3, cores = my.cores, warmup = 5000, iter = 20000,
                               family = binomial("logit"), 
                               save_pars = save_pars(all = TRUE),
                               data = brood.mod.data)

## model check ####

print(brood_sizeM, digits = 3)

pp_check(brood_sizeM, ndraws = 500)
plot(brood_sizeM)

## summary ####
describe_posterior(
  brood_sizeM,
  effects = "all",
  component = "all",
  centrality = "all",
  test = c("p_direction", "p_significance"),
  verbose = FALSE
)

# basic plot
plot_model(
  brood_sizeM,
  ci.lvl = 0.9,
  dot.size = 2,
  line.size = 1,
  rm.terms = "b_Intercept") +
  THEME0 +
  geom_hline(yintercept = 1, linetype = "dashed", alpha = 0.5)

## contrasts ####

rg_brm.Bsize <- ref_grid(brood_sizeM)
em_nnet.Bsize <- emmeans(rg_brm.Bsize,
                         specs = ~brood.size.fin)
plot(em_nnet.Bsize)

c_.Bsize <- contrast(em_nnet.Bsize, interaction = "pairwise",
                     # "consec",
                     simple = "each", combine = TRUE, adjust = "mvt")
c_.Bsize
plot(c_.Bsize, plot = T)

c_.Bsize.df <- plot(c_.Bsize, plot = F)

c_.Bsize_plot <- ggplot(c_.Bsize.df, aes(y = pri.fac, x = the.emmean)) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
  geom_point(size = 3) +
  geom_linerange(aes(xmin=lower.HPD, xmax=upper.HPD), linewidth=1,show.legend = F) +
  labs(y = "Contrasts",
       x = "Estimate") +
  THEME0
c_.Bsize_plot


#_____________________________####
# Predictors of biparental provisioning (Fig. 2A) ####
#_____________________________####

## model ####
my.cores <- detectCores()-1

FlexBipM <- brms::brm(prey.to.nest.F10 | trials(preys.nest10) ~
                        
                        vedba.avg.F_z +
                        hunt.att.sum.F_z + 
                        prop.prey.eaten.prey.capturedF1_z +
                        prop.succ.dives.F_z + 
                        vedba.avg.M_z +
                        hunt.att.sum.M_z +
                        prop.succ.dives.M_z +
                        prop.prey.eaten.prey.capturedM_z +
                        n_nest_encounters_z +
                        n_out_encounters_z +
                        time_at_nest_M_z + 
                        time_at_nest_F_z +
                        SMI_M_z +
                        SMI_F_z +
                        area_wildflower.zsqrt +
                        yearF +
                        brood.size.change01F +
                        youngest.chick.age_z +
                        (1|BroodID), 
                      
                      chains = 3, cores = my.cores, warmup = 5000, iter = 20000, 
                      # control = list(adapt_delta = 0.9),
                      family = binomial("logit"),
                      save_pars = save_pars(all = TRUE),
                      data = dat.mod)

## model check ####
print(FlexBipM, digits = 3)
pp_check(FlexBipM, ndraws = 100)
plot(FlexBipM)


## summary ####

describe_posterior(
  FlexBipM,
  effects = "all",
  component = "all",
  centrality = "all",
  test = c("p_direction", "p_significance"),
  verbose = FALSE
)

## basic plots ####

plot_model(FlexBipM, bpe = "median", bpe.style ="dot", ci.style = "whisker",
                                   prob.inner = 0.5, prob.outer = 0.95) + # forest plot
  THEME1 +
  scale_y_continuous(limits = c(0.2, 2)) +
  geom_hline(yintercept = 1, linetype = "dashed", alpha = 0.4) +
  scale_color_okabe_ito() +
  scale_fill_okabe_ito()



#_____________________________####
# Predictors of PREY PER CHICK (Fig. 2B) ####
#_____________________________####

## model ####

my.cores <- detectCores()-1

PpC_M <- brms::brm(prey.per.chick ~
                     prey.nest.Fcoop01_z + # female provisioning share
                     vedba.avg.F_z +
                     hunt.att.sum.F_z + 
                     prop.prey.eaten.prey.capturedF1_z +
                     prop.succ.dives.F_z + 
                     vedba.avg.M_z +
                     hunt.att.sum.M_z +
                     prop.succ.dives.M_z +
                     prop.prey.eaten.prey.capturedM_z +
                     n_nest_encounters_z +
                     n_out_encounters_z +
                     time_at_nest_M_z + 
                     time_at_nest_F_z +
                     SMI_M_z +
                     SMI_F_z +
                     area_wildflower.zsqrt +
                     yearF +
                     brood.size.change01F +
                     youngest.chick.age_z +
                     (1|BroodID), 
                   chains = 3, cores = my.cores, warmup = 5000, iter = 20000, 
                   # control = list(adapt_delta = 0.9),
                   family = gaussian,
                   save_pars = save_pars(all = TRUE),
                   data = dat.mod)

## model check ####

print(PpC_M, digits = 3)
pp_check(PpC_M, ndraws = 100)
plot(PpC_M)

## summary ####
describe_posterior(
  PpC_M,
  effects = "all",
  component = "all",
  centrality = "all",
  test = c("p_direction", "p_significance"),
  verbose = FALSE
)

## basic plots ####

plot_model(PpC_M, bpe = "median", bpe.style ="dot", ci.style = "whisker",
           prob.inner = 0.5, prob.outer = 0.95) + # forest plot
  THEME1 +
  scale_y_continuous(limits = c(0.2, 2)) +
  geom_hline(yintercept = 1, linetype = "dashed", alpha = 0.4) +
  scale_color_okabe_ito() +
  scale_fill_okabe_ito()


#_______________________________ ####
# FEMALE prey per chick predicted by MALE foraging behaviours (Fig. 2C)####
#_______________________________ ####


## model ####

my.cores <- detectCores()-1

FppcMb <- brms::brm(prey.per.chickF ~
                      vedba.avg.M_z +
                      hunt.att.sum.M_z +
                      prop.succ.dives.M_z +
                      prop.prey.eaten.prey.capturedM_z +
                      time_at_nest_M_z + 
                      n_nest_encounters_z +
                      n_out_encounters_z +
                      area_wildflower.zsqrt +
                      yearF +
                      brood.size.change01F +
                      youngest.chick.age_z +
                      (1|BroodID),
                    chains = 3, cores = my.cores, warmup = 5000, iter = 20000, 
                    # control = list(adapt_delta = 0.9),
                    family = gaussian,
                    save_pars = save_pars(all = TRUE),
                    data = dat.mod)


## model check ####

print(FppcMb, digits = 3)
pp_check(FppcMb, ndraws = 100)
plot(FppcMb)

## summary ####
table.FppcMb <- describe_posterior(
  FppcMb,
  effects = "all",
  component = "all",
  centrality = "all",
  test = c("p_direction", "p_significance"),
  verbose = FALSE
)

## basic plots ####

plot_model(FppcMb, bpe = "median", bpe.style ="dot", ci.style = "whisker",
           prob.inner = 0.5, prob.outer = 0.95) + # forest plot
  THEME1 +
  # scale_y_continuous(limits = c(0.2, 2)) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.4) +
  scale_color_okabe_ito() +
  scale_fill_okabe_ito() +
  labs(title = "")


#_______________________________ ####
# MALE prey per chick predicted by FEMALE foraging behaviours (Fig. 2D)####
#_______________________________ ####

## model ####

my.cores <- detectCores()-1

MppcFb <- brms::brm(prey.per.chickM ~
                      vedba.avg.M_z +
                      hunt.att.sum.M_z +
                      prop.succ.dives.M_z +
                      prop.prey.eaten.prey.capturedM_z +
                      time_at_nest_M_z + 
                      n_nest_encounters_z +
                      n_out_encounters_z +
                      area_wildflower.zsqrt +
                      yearF +
                      brood.size.change01F +
                      youngest.chick.age_z +
                      (1|BroodID),
                    chains = 3, cores = my.cores, warmup = 5000, iter = 20000, 
                    # control = list(adapt_delta = 0.9),
                    family = gaussian,
                    save_pars = save_pars(all = TRUE),
                    data = dat.mod)


## model check ####

print(MppcFb, digits = 3)
pp_check(MppcFb, ndraws = 100)
mcmc_plot(MppcFb, type = "acf") #correlation plots
plot(MppcFb)

## summary ####
describe_posterior(
  MppcFb,
  effects = "all",
  component = "all",
  centrality = "all",
  test = c("p_direction", "p_significance"),
  verbose = FALSE
)

## basic plots ####

plot_model(MppcFb, bpe = "median", bpe.style ="dot", ci.style = "whisker",
           prob.inner = 0.5, prob.outer = 0.95) + # forest plot
  THEME1 +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.4) +
  scale_color_okabe_ito() +
  scale_fill_okabe_ito() +
  labs(title = "")

# _____________________ ####
# Night-to-night carry-over effects ####
# _____________________ ####

# Model 1 - Response: female provisioning share ####
# _____________________ ####

## model ####
my.cores <- detectCores()-1

bip.lag.r1 <- brms::brm(prey.to.nest.F10 | trials(preys.nest10) ~
                          prey.nest.Fcoop01_z_prev +
                          vedba.avg.F_z_prev +
                          hunt.att.sum.F_z_prev + 
                          prop.prey.eaten.prey.capturedF1_z_prev +
                          prop.succ.dives.F_z_prev + 
                          vedba.avg.M_z_prev +
                          hunt.att.sum.M_z_prev +
                          prop.succ.dives.M_z_prev +
                          prop.prey.eaten.prey.capturedM_z_prev +
                          n_nest_encounters_z_prev +
                          n_out_encounters_z_prev +
                          time_at_nest_M_z_prev + 
                          time_at_nest_F_z_prev +
                          (1|BroodID), 
                        chains = 3, cores = my.cores, warmup = 5000, iter = 20000, 
                        # control = list(adapt_delta = 0.9),
                        family = binomial("logit"),
                        save_pars = save_pars(all = TRUE),
                        data = dat.lag)

## model check ####
print(bip.lag.r1, digits = 3)
pp_check(bip.lag.r1, ndraws = 100)
plot(bip.lag.r1)


## summary ####

describe_posterior(
  bip.lag.r1,
  effects = "all",
  component = "all",
  centrality = "all",
  test = c("p_direction", "p_significance"),
  verbose = FALSE
)

## basic plots ####

plot_model(bip.lag.r1, bpe = "median", bpe.style ="dot", ci.style = "whisker",
           prob.inner = 0.5, prob.outer = 0.95) + # forest plot
  THEME1 +
  scale_y_continuous(limits = c(0.2, 2)) +
  geom_hline(yintercept = 1, linetype = "dashed", alpha = 0.4) +
  scale_color_okabe_ito() +
  scale_fill_okabe_ito()

# _____________________ ####
# Model 2 - Response: prey per chick by both parents ####
# _____________________ ####

## model ####
my.cores <- detectCores()-1

ppc.lag.r1 <- brms::brm(prey.per.chick ~
                          prey.per.chick_z_prev +
                          prey.nest.Fcoop01_z_prev +
                          vedba.avg.F_z_prev +
                          hunt.att.sum.F_z_prev + 
                          prop.prey.eaten.prey.capturedF1_z_prev +
                          prop.succ.dives.F_z_prev + 
                          vedba.avg.M_z_prev +
                          hunt.att.sum.M_z_prev +
                          prop.succ.dives.M_z_prev +
                          prop.prey.eaten.prey.capturedM_z_prev +
                          n_nest_encounters_z_prev +
                          n_out_encounters_z_prev +
                          time_at_nest_M_z_prev + 
                          time_at_nest_F_z_prev +
                          (1|BroodID), 
                        chains = 3, cores = my.cores, warmup = 5000, iter = 20000, 
                        # control = list(adapt_delta = 0.9),
                        family = gaussian,
                        save_pars = save_pars(all = TRUE),
                        data = dat.lag)

## model check ####
print(ppc.lag.r1, digits = 3)
pp_check(ppc.lag.r1, ndraws = 100)
plot(ppc.lag.r1)


## summary ####

describe_posterior(
  ppc.lag.r1,
  effects = "all",
  component = "all",
  centrality = "all",
  test = c("p_direction", "p_significance"),
  verbose = FALSE
)

## basic plots ####

plot_model(ppc.lag.r1, bpe = "median", bpe.style ="dot", ci.style = "whisker",
           prob.inner = 0.5, prob.outer = 0.95) + # forest plot
  THEME1 +
  scale_y_continuous(limits = c(0.2, 2)) +
  geom_hline(yintercept = 1, linetype = "dashed", alpha = 0.4) +
  scale_color_okabe_ito() +
  scale_fill_okabe_ito()

# _____________________ ####
# Model 3 - Response: prey per chick by female parent ####
# _____________________ ####

## model ####
my.cores <- detectCores()-1

ppcF.r1 <- brms::brm(prey.per.chickF ~
                       vedba.avg.M_z +
                       hunt.att.sum.M_z +
                       prop.succ.dives.M_z +
                       prop.prey.eaten.prey.capturedM_z +
                       time_at_nest_M_z + 
                       n_nest_encounters_z +
                       n_out_encounters_z +
                       area_wildflower.zsqrt +
                       yearF +
                       brood.size.change01F +
                       youngest.chick.age_z +
                       (1|BroodID),
                     chains = 3, cores = my.cores, warmup = 5000, iter = 20000, 
                     # control = list(adapt_delta = 0.9),
                     family = gaussian,
                     save_pars = save_pars(all = TRUE),
                     data = dat.mod)

## model check ####
print(ppcF.r1, digits = 3)
pp_check(ppcF.r1, ndraws = 100)
plot(ppcF.r1)


## summary ####

describe_posterior(
  ppcF.r1,
  effects = "all",
  component = "all",
  centrality = "all",
  test = c("p_direction", "p_significance"),
  verbose = FALSE
)

## basic plots ####

plot_model(ppcF.r1, bpe = "median", bpe.style ="dot", ci.style = "whisker",
           prob.inner = 0.5, prob.outer = 0.95) + # forest plot
  THEME1 +
  scale_y_continuous(limits = c(0.2, 2)) +
  geom_hline(yintercept = 1, linetype = "dashed", alpha = 0.4) +
  scale_color_okabe_ito() +
  scale_fill_okabe_ito()

# _____________________ ####
# Model 4 - Response: prey per chick by male parent ####
# _____________________ ####

## model ####
my.cores <- detectCores()-1

ppcM.lag.r1 <- brms::brm(prey.per.chickM ~
                           prey.per.chickM_z_prev +
                           vedba.avg.F_z_prev +
                           hunt.att.sum.F_z_prev +
                           prop.succ.dives.F_z_prev +
                           prop.prey.eaten.prey.capturedF1_z_prev +
                           time_at_nest_F_z_prev + 
                           n_nest_encounters_z_prev +
                           n_out_encounters_z_prev +
                           (1|BroodID), 
                         chains = 3, cores = my.cores, warmup = 5000, iter = 20000, 
                         control = list(adapt_delta = 0.9999),
                         family = gaussian,
                         save_pars = save_pars(all = TRUE),
                         data = dat.lag)

## model check ####
print(ppcM.lag.r1, digits = 3)
pp_check(ppcM.lag.r1, ndraws = 100)
plot(ppcM.lag.r1)


## summary ####

describe_posterior(
  ppcM.lag.r1,
  effects = "all",
  component = "all",
  centrality = "all",
  test = c("p_direction", "p_significance"),
  verbose = FALSE
)

## basic plots ####

plot_model(ppcM.lag.r1, bpe = "median", bpe.style ="dot", ci.style = "whisker",
           prob.inner = 0.5, prob.outer = 0.95) + # forest plot
  THEME1 +
  scale_y_continuous(limits = c(0.2, 2)) +
  geom_hline(yintercept = 1, linetype = "dashed", alpha = 0.4) +
  scale_color_okabe_ito() +
  scale_fill_okabe_ito()
