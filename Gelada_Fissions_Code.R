## CODE FOR GELADA FISSIONS MANUSCRIPT
## LOAD REQUIRED PACKAGES
library(ggplot2)
library(brms)
library(survival)
library(tidyverse)
library(tidybayes)
library(viridis)
library(patchwork)
theme_set(theme_classic(base_size = 20))

## USE DEVELOPMENT VERSION WITH SURVIVAL FUNCTIONS
## How to install:
## install.packages("rstanarm", repos = c('https://stan-dev.r-universe.dev', getOption("repos")))
library(rstanarm)

## SET WORKING DIRECTORY
setwd("~/Desktop/SBU Postdoc/Manuscripts/Fission - Proc B/GitHub")

## LOAD REQUIRED DATASETS
## FEMALE-FEMALE DYADS
ff.dyads<-read.csv("Fission_Dyadic_FF.csv")

## FEMALE-MALE DYADS
fm.dyads<-read.csv("Fission_Dyadic_FM.csv")

## AGGRESSION
agg.data<-read.csv("Fission_Agonism.csv")

## ACTIVITY BUDGETS
state.data<-read.csv("Fission_States.csv")

## UNIT SIZES
unit.sizes<-read.csv("Fission_Units.csv")

## INFANT SURVIVAL
inf.surv<-read.csv("Fission_Infant_Survival.csv")

## INFANT MORTALITY CAUSES
inf.mortalities<-read.csv("Fission_Infant_Mortality_Causes.csv")

## PART 1: DYADIC REGRESSIONS
## BUILD DYADIC ANALYSES

## CONSTRUCT BISON MODELS FROM FEMALE-FEMALE DYADIC DATA
library(bisonR)
priors <- get_default_priors("binary")
priors
priors$edge <- "normal(-5, 2)"
prior_check(priors,"binary")

fit_edge_D <- bison_model(
  (Groom.Total|Focal.Total) ~ dyad(ID1, ID2), 
  data=ff.dyads[ff.dyads$Unit=="D",], 
  model_type="binary",
  priors=priors, directed = F)
plot_predictions(fit_edge_D, num_draws=20, type="density")

fit_edge_C <- bison_model(
  (Groom.Total|Focal.Total) ~ dyad(ID1, ID2), 
  data=ff.dyads[ff.dyads$Unit=="C",], 
  model_type="binary",
  priors=priors, directed = F)
plot_predictions(fit_edge_C, num_draws=20, type="density")

fit_edge_M <- bison_model(
  (Groom.Total|Focal.Total) ~ dyad(ID1, ID2), 
  data=ff.dyads[ff.dyads$Unit=="M",], 
  model_type="binary",
  priors=priors, directed = F)
plot_predictions(fit_edge_M, num_draws=20, type="density")

fit_edge_Z <- bison_model(
  (Groom.Total|Focal.Total) ~ dyad(ID1, ID2), 
  data=ff.dyads[ff.dyads$Unit=="Z",], 
  model_type="binary",
  priors=priors, directed = F)

fit_edge_T <- bison_model(
  (Groom.Total|Focal.Total) ~ dyad(ID1, ID2), 
  data=ff.dyads[ff.dyads$Unit=="T",], 
  model_type="binary",
  priors=priors, directed = F)

fit_edge_B <- bison_model(
  (Groom.Total|Focal.Total) ~ dyad(ID1, ID2), 
  data=ff.dyads[ff.dyads$Unit=="B",], 
  model_type="binary",
  priors=priors, directed = F)
plot_predictions(fit_edge_B, num_draws=20, type="density")

## MAKE FUNCTION FOR Z-SCORING WITHIN UNITS
group.center <- function(var,grp) {
  return((var-tapply(var,grp,mean,na.rm=T)[grp])/
           tapply(var,grp,sd,na.rm=T)[grp])
  
}

## EXTRACT DYADIC MEASUREMENTS
## DATA LIST
bison.list<-list()
for (i in 1:20) {
  edges<-c(extract_metric(fit_edge_B,"edge_weight")[i,],
           extract_metric(fit_edge_C,"edge_weight")[i,],
           extract_metric(fit_edge_D,"edge_weight")[i,],
           extract_metric(fit_edge_M,"edge_weight")[i,],
           extract_metric(fit_edge_T,"edge_weight")[i,],
           extract_metric(fit_edge_Z,"edge_weight")[i,])
  new.data<-ff.dyads
  j<-1
  for (j in 1:nrow(new.data)) {
    new.data$groom.value[j]<-edges[j]
  }
  new.data$groom.z<-group.center(new.data$groom.value, new.data$Unit)
  bison.list[[i]]<-new.data
}

## CONSTRUCT STRUCTURAL EQUATION MODELS
f1<-bf(groom.z ~ scale(Age.Diff) + scale(Pedigree.R) + scale(RankDiff) + (1|mm(ID1,ID2)), family="student") 
f2<-bf(Fissioned~ scale(Age.Diff) + scale(Pedigree.R) + scale(RankDiff) + scale(groom.z) + (1|mm(ID1,ID2)), family="bernoulli")

get_prior(f1+f2 + set_rescor(FALSE), data=bison.list)

m2.females <- brm_multiple(f1 + f2 + set_rescor(FALSE), 
                           prior = c(prior(normal(0, 2), class="b", resp="groomz"),
                                     prior(normal(0, 2), class="b", resp="Fissioned")), data = bison.list[1:20], cores = 4, 
                           control=list(adapt_delta=0.99))
summary(m2.females, prob=0.89)

## MAKE CONDITIONAL EFFECTS PLOTS
plot1<-plot(conditional_effects(m2.females, "Pedigree.R", resp="Fissioned", spaghetti = T, ndraws=250), spaghetti_args = list(colour = scales::alpha("#21908CFF", 0.2), linewidth=1)) + ylab("P(Stayed together)") + xlab("Relatedness")
plot2<-plot(conditional_effects(m2.females, "RankDiff", resp="Fissioned", spaghetti = T, ndraws=250), spaghetti_args = list(colour = scales::alpha("#21908CFF", 0.2), linewidth=1)) + xlab("Rank difference")
plot3<-plot(conditional_effects(m2.females, "groom.z", resp="Fissioned", spaghetti = T, ndraws=250), spaghetti_args = list(colour = scales::alpha("#21908CFF", 0.2), linewidth=1)) + xlab("Grooming relationship")

plot1<-conditional_effects(m2.females, "Pedigree.R", resp="Fissioned", prob=0.89)
plot1$Fissioned.Fissioned_Pedigree.R$estimate__
plot1$Fissioned.Fissioned_Pedigree.R$estimate__<-1-plot1$Fissioned.Fissioned_Pedigree.R$estimate__
plot1$Fissioned.Fissioned_Pedigree.R$lower__
plot1$Fissioned.Fissioned_Pedigree.R$lower__<-1-plot1$Fissioned.Fissioned_Pedigree.R$lower__
plot1$Fissioned.Fissioned_Pedigree.R$upper__
plot1$Fissioned.Fissioned_Pedigree.R$upper__<-1-plot1$Fissioned.Fissioned_Pedigree.R$upper__
viridis(3)
bison_fission_related<-plot(plot1, plot = FALSE,  line_args = list(colour="black", fill="#21908CFF"))[[1]] +
  ylim(0,1) + ylab("P(stayed together)") + xlab("Relatedness") 
bison_fission_related

plot2<-conditional_effects(m2.females, "RankDiff", resp="Fissioned", prob=0.89)
plot2$Fissioned.Fissioned_RankDiff$estimate__
plot2$Fissioned.Fissioned_RankDiff$estimate__<-1-plot2$Fissioned.Fissioned_RankDiff$estimate__
plot2$Fissioned.Fissioned_RankDiff$lower__
plot2$Fissioned.Fissioned_RankDiff$lower__<-1-plot2$Fissioned.Fissioned_RankDiff$lower__
plot2$Fissioned.Fissioned_RankDiff$upper__
plot2$Fissioned.Fissioned_RankDiff$upper__<-1-plot2$Fissioned.Fissioned_RankDiff$upper__
bison_fission_rank<-plot(plot2, plot = FALSE,  line_args = list(colour="black", fill="#21908CFF"))[[1]] +
  ylim(0,1) + ylab("") + xlab("Rank difference") 
bison_fission_rank

plot3<-conditional_effects(m2.females, "groom.z", resp="Fissioned", prob=0.89)
plot3$Fissioned.Fissioned_groom.z$estimate__
plot3$Fissioned.Fissioned_groom.z$estimate__<-1-plot3$Fissioned.Fissioned_groom.z$estimate__
plot3$Fissioned.Fissioned_groom.z$lower__
plot3$Fissioned.Fissioned_groom.z$lower__<-1-plot3$Fissioned.Fissioned_groom.z$lower__
plot3$Fissioned.Fissioned_groom.z$upper__
plot3$Fissioned.Fissioned_groom.z$upper__<-1-plot3$Fissioned.Fissioned_groom.z$upper__
bison_fission_groom<-plot(plot3, plot = FALSE,  line_args = list(colour="black", fill="#21908CFF"))[[1]] +
  ylim(0,1) + ylab("") + xlab("Grooming relationship\n(z-scored)")
bison_fission_groom

## FIGURE 2
bison_fission_related + bison_fission_rank + bison_fission_groom + plot_annotation(tag_levels = "a")
theme_set(theme_classic(base_size = 16))

## CALCULATE INDIRECT EFFECTS
## KINSHIP EFFECT
draws <- as_draws_df(m2.females)
draws <- draws %>%
  mutate(
    Indirect_Effect = b_groomz_scalePedigree.R * b_Fissioned_scalegroom.z, # Path a * Path b
    Direct_Effect = -b_Fissioned_scalePedigree.R, # Indirect + Path c'
    Total_Effect = Indirect_Effect + b_Fissioned_scalePedigree.R, # Indirect + Path c'
    Percent_Mediated = (Indirect_Effect / Total_Effect) * 100
  )

quantile(draws$Indirect_Effect, probs=c(0.025,0.975))
quantile(draws$Percent_Mediated, probs=c(0.025,0.975))
quantile(draws$Percent_Mediated, probs=c(0.055,0.5,0.945))

kinship.effects.plot<-ggplot(data=draws) + 
  geom_density(aes(x=-Indirect_Effect), fill="#21908CFF", alpha=0.5) +
  geom_density(aes(x=Direct_Effect), fill="#440154FF", alpha=0.5) +
  #geom_density(aes(x=-Total_Effect), fill="gray40", alpha=0.5) +
  xlim(-0.5,2) + xlab("Kinship effects") + ylim(0,10)

## RANK EFFECT
draws <- draws %>%
  mutate(
    Indirect_Effect = b_groomz_scaleRankDiff * b_Fissioned_scalegroom.z, # Path a * Path b
    Direct_Effect = b_Fissioned_scaleRankDiff, # Indirect + Path c'
    Total_Effect = Indirect_Effect + b_Fissioned_scaleRankDiff, # Indirect + Path c'
    Percent_Mediated = (Indirect_Effect / Total_Effect) * 100
  )

quantile(draws$Indirect_Effect, probs=c(0.025,0.975))
quantile(draws$Percent_Mediated, probs=c(0.025,0.975))
quantile(draws$Percent_Mediated, probs=c(0.055,0.5,0.945))

rank.effects.plot<-ggplot(data=draws) + 
  geom_density(aes(x=Indirect_Effect), fill="#21908CFF", alpha=0.5) +
  geom_density(aes(x=Direct_Effect), fill="#440154FF", alpha=0.5) +
  #geom_density(aes(x=Total_Effect), fill="gray40", alpha=0.5) +
  xlim(-0.5,2) + xlab("Rank effects") + ylim(0,10)
kinship.effects.plot+rank.effects.plot

## MEASURE FEMALE-MALE DYADS
## USING BISON MODELS
## RUN MALE MODELS BISON STYLE
fit_males_D <- bison_model(
  (Groom.Total|Focal.Total) ~ dyad(ID1, ID2), 
  data=fm.dyads[fm.dyads$Unit=="D",], 
  model_type="binary",
  priors=priors, directed = F)
plot_predictions(fit_males_D, num_draws=20, type="density")

fit_males_C <- bison_model(
  (Groom.Total|Focal.Total) ~ dyad(ID1, ID2), 
  data=fm.dyads[fm.dyads$Unit=="C",], 
  model_type="binary",
  priors=priors, directed = F)
plot_predictions(fit_males_C, num_draws=20, type="density")

fit_males_M <- bison_model(
  (Groom.Total|Focal.Total) ~ dyad(ID1, ID2), 
  data=fm.dyads[fm.dyads$Unit=="M",], 
  model_type="binary",
  priors=priors, directed = F)
plot_predictions(fit_males_M, num_draws=20, type="density")

fit_males_T <- bison_model(
  (Groom.Total|Focal.Total) ~ dyad(ID1, ID2), 
  data=fm.dyads[fm.dyads$Unit=="T",], 
  model_type="binary",
  priors=priors, directed = F)
plot_predictions(fit_males_T, num_draws=20, type="density")

fit_males_Z <- bison_model(
  (Groom.Total|Focal.Total) ~ dyad(ID1, ID2), 
  data=fm.dyads[fm.dyads$Unit=="Z",], 
  model_type="binary",
  priors=priors, directed = F)
plot_predictions(fit_males_Z, num_draws=20, type="density")

fit_males_B <- bison_model(
  (Groom.Total|Focal.Total) ~ dyad(ID1, ID2), 
  data=fm.dyads[fm.dyads$Unit=="B",], 
  model_type="binary",
  priors=priors, directed = F)
plot_predictions(fit_males_B, num_draws=20, type="density")

## GET ESTIMATES
bison.list.males<-list()
i<-1
for (i in 1:20) {
  edges<-c(extract_metric(fit_males_B,"edge_weight")[i,],
           extract_metric(fit_males_C,"edge_weight")[i,],
           extract_metric(fit_males_D,"edge_weight")[i,],
           extract_metric(fit_males_M,"edge_weight")[i,],
           extract_metric(fit_males_T,"edge_weight")[i,],
           extract_metric(fit_males_Z,"edge_weight")[i,])
  new.data<-fm.dyads
  for (j in 1:nrow(new.data)) {
    new.data$groom.value[j]<-edges[j]
  }
  new.data$groom.value<-group.center(new.data$groom.value, new.data$Unit)
  bison.list.males[[i]]<-new.data
}

## CONSTRUCT STRUCTURAL EQUATION MODEL
f1<-bf(groom.value ~ Infant + (1|mm(ID1,ID2)), family="student") 
f2<-bf(Fissioned~ Infant + scale(groom.value) + (1|mm(ID1,ID2)), family="bernoulli")

get_prior(f1+f2 + set_rescor(FALSE), data=bison.list.males[[1]])

bison.list.males[[1]]$Fissioned

m2.males <- brm_multiple(f1 + f2 + set_rescor(FALSE), 
                         prior = c(prior(normal(0, 2), class="b", resp="groomvalue"),
                                   prior(normal(0, 2), class="b", resp="Fissioned")), data = bison.list.males[1:5], cores = 4, 
                         control=list(adapt_delta=0.99))
summary(m2.males, prob=0.89)

## PART 2: AGGRESION DATA
agg.data$Period<-as.factor(agg.data$Period)
agg.data$Period<-factor(agg.data$Period, levels=c("Stable","Pre","Post"))

agg.model<-brm(data=agg.data,Total_Aggression~ Period + offset(log(N.Focal)) + (1||Subject), 
            family="poisson", 
            iter=4000, cores=2, prior = c(
              prior(normal(0, 2), "b"),
              prior(normal(0, 5), "Intercept"),
              prior(student_t(3, 0, 1), "sd")),
            control = list(adapt_delta=0.99))
summary(agg.model, prob=0.89)

## MAKE AGONISM PLOT
agg.plota<-conditional_effects(agg.model)

agg.plota$Period$estimate__<-agg.plota$Period$estimate__/agg.plota$Period$N.Focal
agg.plota$Period$lower__<-agg.plota$Period$lower__/agg.plota$Period$N.Focal
agg.plota$Period$upper__<-agg.plota$Period$upper__/agg.plota$Period$N.Focal

eff_data <- as.data.frame(conditional_effects(agg.model)$Period)

viridis(3, begin = 1, end = 0)

palette<-c("#FDE725FF","#21908CFF","#440154FF")
agg.plota<-ggplot(data=eff_data) +
  geom_errorbar(aes(x=Period, ymin=lower__/N.Focal, ymax=upper__/N.Focal, color=Period), width=0.3, cex=2) +
  geom_point(aes(x=Period, y=estimate__/N.Focal, color=Period), size=4) +
  scale_color_manual(values = palette) + ylab("Agonism (events/hour)") + xlab("Period") +
  scale_y_continuous(limits=c(0,0.5), breaks=seq(0,0.5,0.25)) + theme(legend.position = "none") +
  scale_x_discrete(labels=c("Stable", "Pre-fission", "Post-fission")) 
agg.plota

## PART 3: ACTIVITY BUDGETS
## VISUALLY SUMMARIZE ACTIVITY
state.data$Period<-as.factor(state.data$Period)
state.data$Period<-factor(state.data$Period, levels=c("Pre","Fission","Post"))

ID.states.total<-state.data %>%
  dplyr::group_by(Period) %>%
  dplyr::summarise(Feed.Prop=sum(Total.Feed)/sum(Total.Time),
                   Rest.Prop=sum(Total.Rest)/sum(Total.Time),
                   Move.Prop=sum(Total.Move)/sum(Total.Time),
                   Social.Prop=sum(Total.Social)/sum(Total.Time))

ID.states.total<-reshape2::melt(ID.states.total)

budget.plot<-ggplot(ID.states.total, aes(x = Period, y = value, fill = variable)) +
  geom_col(position = "fill", alpha=0.6, color="black") +
  scale_fill_viridis(discrete=T,labels=c("Feed","Rest","Move","Social"), begin = 0.7, end = 0) +
  scale_x_discrete(labels=c("Stable", "Pre-fission","Post-fission")) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "Period", y = "Proportion of activity budget", fill = "Category")

## FIGURE 3
agg.plota+budget.plot + plot_annotation(tag_levels = "a")

## CONSTRUCT MODELS FOR EACH BEHAVIORAL STATE
## MODEL INDIVIDUAL PROPORTIONS
activity.model.feed<-brm(data=state.data,
                         Total.Feed|trials(Total.Time) ~ Period + Season + (1|focalanimal_code),
                         chains=4, iter=4000, cores=2,
                         prior = c(
                           prior(normal(0, 2), "b"),
                           prior(normal(0, 5), "Intercept"),
                           prior(student_t(3, 0, 1), "sd")),
                         family="binomial")
summary(activity.model.feed, prob=0.89)

activity.model.move<-brm(data=state.data,
                         Total.Move|trials(Total.Time) ~ Period + Season + (1|focalanimal_code),
                         chains=4, iter=4000, cores=2,
                         prior = c(
                           prior(normal(0, 2), "b"),
                           prior(normal(0, 5), "Intercept"),
                           prior(student_t(3, 0, 1), "sd")),
                         family="binomial")
pp_check(activity.model.move)
summary(activity.model.move)

activity.model.rest<-brm(data=state.data,
                         Total.Rest|trials(Total.Time) ~ Period + Season  + (1|focalanimal_code),
                         chains=4, iter=8000, cores=2,
                         prior = c(
                           prior(normal(0, 2), "b"),
                           prior(normal(0, 5), "Intercept"),
                           prior(student_t(3, 0, 1), "sd")),
                         family="binomial")
pp_check(activity.model.rest)
summary(activity.model.rest)

activity.model.social<-brm(data=state.data,
                           Total.Social|trials(Total.Time) ~ Period + Season  + (1|focalanimal_code),
                           chains=4, iter=8000, cores=2,
                           prior = c(
                             prior(normal(0, 2), "b"),
                             prior(normal(0, 5), "Intercept"),
                             prior(student_t(3, 0, 1), "sd")),
                           family="binomial")
pp_check(activity.model.social)
summary(activity.model.social, prob=0.89)

## PART 4: INFANT SURVIVORSHIP
## GROUP SIZE COMPARISONS
## FOR VISUALIZATION
unit.sizes$N.Females[grepl("T", unit.sizes$Unit) & unit.sizes$N.Females<5]<-5
unit.sizes.rep<-data.frame(N.Females=rep(unit.sizes$N.Females, unit.sizes$N.Females), 
                           Fissioned=rep(unit.sizes$Category, unit.sizes$N.Females),
                           Unit=rep(unit.sizes$Unit, unit.sizes$N.Females),
                           Year=rep(unit.sizes$Year, unit.sizes$N.Females))

unit.sizes.year<-unit.sizes %>%
  group_by(Unit, Year, Category) %>%
  dplyr::summarise(N.Females=mean(N.Females))

ks.test(unit.sizes.year$N.Females[unit.sizes.year$Category=="Stable"], unit.sizes.year$N.Females[unit.sizes.year$Category=="Pre"])
ks.test(unit.sizes.year$N.Females[unit.sizes.year$Category=="Stable"], unit.sizes.year$N.Females[unit.sizes.year$Category=="Post"])

library(ggridges)
plot3a<-ggplot(unit.sizes.rep, aes(x=N.Females, y=Fissioned, fill=Fissioned)) + stat_density_ridges(bandwidth=1, scale=0.9, alpha=0.6) +
  scale_fill_viridis(discrete=T) + xlab("# females") + ylab("") +
  annotate("rect", xmin=4.5, xmax=7.5, ymin=0, ymax=Inf, alpha=0.4, fill="gray", color="black") + theme(legend.position = "none") +
  annotate("text", x=12.2, y=0.51, label="optimal range")  + geom_segment(aes(x = 8.8, y = 0.5, xend = 7.6, yend = 0.5), arrow = arrow(length = unit(0.3, "cm"))) +   coord_cartesian(clip = "off") + scale_x_continuous(breaks = seq(2,14,4))
plot3a

## SURVIVORSHIP MODELS
## ADD NOMINAL VALUE TO INFANTS WHO DIED AT BIRTH
inf.surv$End.Age[inf.surv$End.Age==0]<-0.0001

## MAKE SURVIVAL OBJECT
surv_object<-Surv(time=inf.surv$End.Age, event=inf.surv$Censored)
plot(surv_object)

survival.model <- rstanarm::stan_surv(surv_object ~  Fission + Sex + (1|Mom) + (1|Unit) + (1|Year), data = inf.surv, cores=2, basehaz = 'weibull')
summary(survival.model, digits=2, prob=c(0.055,0.945))

## VISUALIZE SURVIVAL CURVES
infant.fissions <- dplyr::select(inf.surv, Fission, Sex, Mom, Unit, Year) %>%
  distinct() %>%
  dplyr::mutate(id = 1:n(),
                TimeCat2 = Fission) %>%
  tidyr::nest_legacy(-c(TimeCat2)) %>%
  mutate(., preds = map(
    data,
    ~ posterior_survfit(
      survival.model, prob=0.89,
      newdata = .x,
      times = 0,
      standardise = TRUE,
      extrapolate = TRUE,
      dynamic = TRUE
    )
  ))

infant.fissions <- tidyr::unnest(infant.fissions, preds) %>%
  select(-data)

viridis(3)
library(plyr)
infant.fissions$TimeCat2<-revalue(infant.fissions$TimeCat2, c("Pre"="Pre-fission", "Post"="Post-fission", "Stable"="Stable"))
infant.fissions$TimeCat2<-factor(infant.fissions$TimeCat2, levels=c("Pre-fission","Post-fission", "Stable"))
#infant.fissions<-infant.fissions[!infant.fissions$TimeCat2=="Stable",]

plot3b<-ggplot(infant.fissions,
               aes(
                 time,
                 fill = TimeCat2,
                 group = TimeCat2
               )) + labs(fill="") + theme(legend.position = "none") +
  geom_ribbon(
    aes(time, ymin = ci_lb, ymax = ci_ub),
    col = NA,
    alpha = 0.6
  ) +   geom_line(aes(x=time, y = median)) +
  scale_fill_manual(values = c("#21908CFF","#440154FF", "#FDE725FF")) +
  ylim(c(0.5, 1)) + xlim(0,1.5) +
  xlab("infant age (years)") + ylab("proportion of infants surviving")
plot3b

## PART 5: INFANT MORTALITY CAUSES
model.infanticides<-brm(data=inf.mortalities, Infanticide ~ Fissioned + (1|Unit), iter=4000, family="bernoulli")
summary(model.infanticides, prob=0.89)
conditional_effects(model.infanticides)

sum(inf.mortalities$Infanticide[inf.mortalities$Fissioned=="Stable"])/nrow(inf.mortalities[inf.mortalities$Fissioned=="Stable",])
sum(inf.mortalities$Infanticide[inf.mortalities$Fissioned=="Pre-fission"])/nrow(inf.mortalities[inf.mortalities$Fission=="Pre-fission",])
sum(inf.mortalities$Infanticide[inf.mortalities$Fissioned=="Post-fission"])/nrow(inf.mortalities[inf.mortalities$Fission=="Post-fission",])

# Model 2: Maternal death
model.mom<-brm(data=inf.mortalities, MaternalDeath ~ Fissioned + (1|Unit), control=list(adapt_delta=0.99), iter=4000, family="bernoulli")
summary(model.mom)
conditional_effects(model.mom, prob=0.89)

sum(inf.mortalities$MaternalDeath[inf.mortalities$Fissioned=="Stable"])/nrow(inf.mortalities[inf.mortalities$Fissioned=="Stable",])
sum(inf.mortalities$MaternalDeath[inf.mortalities$Fissioned=="Pre-fission"])/nrow(inf.mortalities[inf.mortalities$Fission=="Pre-fission",])
sum(inf.mortalities$MaternalDeath[inf.mortalities$Fissioned=="Post-fission"])/nrow(inf.mortalities[inf.mortalities$Fission=="Post-fission",])

# Model 3: Unknown deaths
model.unk<-brm(data=inf.mortalities, Unknown ~ Fissioned + (1|Unit), iter=4000, family="bernoulli")
summary(model.unk, 0.89)
conditional_effects(model.unk, prob=0.89)

sum(inf.mortalities$Unknown[inf.mortalities$Fissioned=="Stable"])/nrow(inf.mortalities[inf.mortalities$Fissioned=="Stable",])
sum(inf.mortalities$Unknown[inf.mortalities$Fissioned=="Pre-fission"])/nrow(inf.mortalities[inf.mortalities$Fission=="Pre-fission",])
sum(inf.mortalities$Unknown[inf.mortalities$Fissioned=="Post-fission"])/nrow(inf.mortalities[inf.mortalities$Fission=="Post-fission",])

