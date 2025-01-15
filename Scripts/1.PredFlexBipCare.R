
# Paper: Dynamic parental roles revealed by fine-scale hunting behaviour with concurrent pair tracking in the wild #
## Script relative to the Result sub-sections: 
### "Variation in biparental care", 
### "Predictors of flexible biparental care" and 
### "Relative and combined parental chick provisioning"
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

# Load data ####

mypath <- "your/path/"
subset.night.z <- read.csv(paste0(mypath, "nightly_params.csv")) # movement/foraging parameters relative to the Male and Female of the pair averaged by night
subset.broodID.z <- read.csv(paste0(mypath, "Pair_params.csv"))  # movement/foraging parameters relative to the Male and Female of the pair averaged by individual


# ______________####
# Repeatability ####
# ______________####

rep.biparentality <- rpt(cbind(prey.to.nest.F10, prey.to.nest.M10) ~ brood.size.fin + yearF + (1| BroodID), 
                 grname = c("BroodID"), data = subset.night.z, 
                 datatype = "Proportion", link = "logit", nboot = 1000, npermut = 0, ratio = T)
print(rep.biparentality)
plot(rep.biparentality)

summary(rep.biparentality)

# __________________________________________________####
# Relationship between BROOD SIZE and BIPARENTALITY ####
# __________________________________________________####

## model ####
my.cores <- detectCores()-1

brood_sizeM <- brms::brm(prey.to.nest.F10 | trials(preys.nest10) ~
                                 0 + brood.size.fin +
                                 (1|BroodID),
                               chains = 3, cores = my.cores, warmup = 5000, iter = 20000,
                               family = binomial("logit"), 
                               save_pars = save_pars(all = TRUE),
                               data = subset.broodID.z)

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
# Predictors of BIPARENTALITY (Fig. 2A) ####
#_____________________________####


## model ####
my.cores <- detectCores()-1

FlexBipM <- brms::brm(prey.to.nest.F10 | trials(preys.nest10) ~
                              
                              vedba.avg.F_z +
                              hunt.att.sum.F_z + 
                              prop.prey.eaten.prey.capturedF_z +
                              prop.succ.dives.F_z + 
                              vedba.avg.M_z +
                              hunt.att.sum.M_z +
                              prop.succ.dives.M_z +
                              prop.prey.eaten.prey.capturedM_z +
                              area_wildflower.zsqrt +
                              yearF +
                              brood.size.change01F +
                              (1|BroodID), 
                            
                            chains = 3, cores = my.cores, warmup = 5000, iter = 20000, 
                            control = list(adapt_delta = 0.9),
                            family = binomial("logit"),
                            save_pars = save_pars(all = TRUE),
                            data = subset.night.z)

## model check ####
print(FlexBipM, digits = 3)
pp_check(FlexBipM, ndraws = 100)
mcmc_plot(FlexBipM, type = "acf") #correlation plots
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
  scale_fill_okabe_ito() +
  labs(title = "Response: Biparentality")

conditional_effects(FcoopMFb.slopeM)

#_____________________________####
# Predictors of PREY PER CHICK (Fig. 2B) ####
#_____________________________####

## model ####

my.cores <- detectCores()-1

PpC_M <- brms::brm(prey.per.chick ~
                               prey.nest.Fcoop01_z + # biparentality
                               vedba.avg.F_z +
                               hunt.att.sum.F_z +
                               prop.prey.eaten.prey.capturedF_z +
                               prop.succ.dives.F_z +
                               vedba.avg.M_z +
                               hunt.att.sum.M_z +
                               prop.succ.dives.M_z +
                               prop.prey.eaten.prey.capturedM_z +
                               area_wildflower.zsqrt +
                               yearF +
                               brood.size.change01F +
                               (1|BroodID) +
                               ar(time = night_nb_real, gr = BroodID, cov = T),
                             chains = 3, cores = my.cores, warmup = 5000, iter = 20000, 
                             family = gaussian, 
                             save_pars = save_pars(all = TRUE),
                             data = subset_coop.ppc)

## model check ####

print(PpC_M, digits = 3)
pp_check(PpC_M, ndraws = 100)
mcmc_plot(PpC_M, type = "acf") #correlation plots
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
                      area_wildflower.zsqrt +
                      yearF +
                      brood.size.change01F +
                      (1|BroodID) +
                      ar(time = night_nb_real, gr = BroodID, cov = T),
                    chains = 3, cores = my.cores, warmup = 5000, iter = 20000, 
                    family = gaussian,
                    save_pars = save_pars(all = TRUE),
                    data = subset_coop.ppcF)


## model check ####

print(FppcMb, digits = 3)
pp_check(FppcMb, ndraws = 100)
mcmc_plot(FppcMb, type = "acf") 
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
                      vedba.avg.F_z +
                      hunt.att.sum.F_z +
                      prop.prey.eaten.prey.capturedF_z +
                      prop.succ.dives.F_z +
                      area_wildflower.zsqrt +
                      yearF +
                      brood.size.change01F +
                      (1|BroodID) +
                      ar(time = night_nb_real, gr = BroodID, cov = T),
                    chains = 3, cores = my.cores, warmup = 5000, iter = 20000, 
                    family = gaussian, 
                    save_pars = save_pars(all = TRUE),
                    data = subset_coop.ppcF)


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

