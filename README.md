# Flexible biparental care in barn owls
Data and script to replicate analysis and figures of the preprint titled "Dynamic parental roles revealed by fine-scale hunting behaviour with concurrent pair tracking in the wild" in EcoEvoRxiv (submitted). 

See Method sections "Behavioural classification and variables" and "Statistical analysis" before running the scripts.

## Tables (.csv)
- Pair_params: It contains all variables relative to the male and female of each one of the 68 pairs. Movement/foraging variables were averaged or summed per individal
- nightly_params: Same as "Pair_params" but movement/foraging variables were averaged by night.
- pamm_table: Table arranged to run time-to-event model using PAMMs
- chicks_table: Table containing measuraments of chicks (several measuraments per chicks) annotated with their parents' biparentality and other movement/foraging parameters.

## Scripts (.R)
- 1.PredFlexBipCare: Script relative to the Result sub-sections: "Variation in biparental care", "Predictors of flexible biparental care" and "Relative and combined parental chick provisioning"
- 2.PAMM_chick-growth-survival: Script relative to the Result sub-sections: "Temporal dynamics of foraging probability", "Contributions of specific parental behaviours to nestling survival" and "Impact of flexible biparental care on chick growth"
