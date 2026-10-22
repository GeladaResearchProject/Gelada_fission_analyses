# Data from: Social and reproductive benefits of group fission in wild geladas

This 15-year dataset was used to analyze the drivers of group choice, changes in aggression and activity budgets, and variation in infant survivorship across fissions in geladas. 

Behavioral data stem from individual focal observations that were divided into three periods: stable (4-2 years prior to fission), pre-fission (2-0 years prior to fission), and post-fission (0-2 years after fission) to track changes in behavior. We also used pre-fission grooming relationships to determine whether stronger grooming partners were more likely to join the same units post-fission.

Demographic data comprise infant survival records that were previously published in Tinsley Johnson et al. 2021.

Tinsley Johnson, E., Feder, J. A., Bergman, T. J., Lu, A., Snyder-Mackler, N., & Beehner, J. C. (2021). The Goldilocks effect: female geladas in mid-sized groups have higher fitness. Proc Biol Sci 288: 20210820

## Description of the data and file structure
“Fission_Dyadic_FF.csv” provides data on grooming interactions among female-female dyads in the two years prior to fission, noting also the dyad’s kinship, rank difference, age difference, and whether they joined the same group following fission. 
 
“Fission_Dyadic_FM.csv” provides data on grooming interactions among female-male dyads in the two years prior to fission, noting also whether the female had a dependent infant at the time of fission and whether the dyad joined the same group after fission.

“Fission_Agonism.csv” notes the number of aggressive interactions received by each focal female during the stable, pre-fission, and post-fission periods, as well as sampling efforts to account for variation in observation and calculate/model rates.

“Fission_States.csv” notes the proportion of time females spent in various forms of activity (feeding, moving, resting, socializing) in the stable, pre-fission, and post-fission periods.

“Fission_Units.csv” notes the number of females in each study unit in each month of the study. These values were used to visually compare unit sizes across stable, pre-fission, and post-fission periods.

“Fission_Infant_Survival.csv” notes the survivorship of infants conceived during stable, pre-fission, and post-fission periods.

“Fission_Infant_Survival.csv” notes the causes of infant mortality (infanticide, maternal death, unknown) of infants conceived during stable, pre-fission, and post-fission periods.

## Code/Software
All analyses were completed in R v4.5.2 using a Bayesian framework, 'rstan' compiling, and a suite of corresponding packages (‘bisonR’, ‘rstanarm’, ‘brms’).

“Gelada_Fissions_Code.R” provides the scripts to construct (i) Bayesian network models, (ii) structural equation models examining the drivers of group co-membership, (iii) mixed effects models focused on changes in aggression rates, (iv) mixed effects models examining changes in activity budgets, (v) Bayesian infant survival models, and (vi) binomial models examining variation in the causes of infant mortality. 
