# Flexible biparental care in barn owls
Data and script to replicate analysis and figures of the unpubished study titled "Real-time coordination of parental provisioning revealed by high-resolution biologging in the wild" under review.

See Method sections before running the scripts.

Raw GPS data are on movebank.org under the project named “Barn owl (Tyto alba)”, ID 231741797.

Raw movement data annotated with behavioural classification from accelerometer data can be found at https://github.com/kimschalcher/data-availability-Schalcher-et-al-eLife (firstly used in Schalcher et al. 2024 eLife https://elifesciences.org/articles/87775).


## Tables (.csv)
- Pair_params: It contains all variables relative to the male and female of each one of the 68 pairs. Movement/foraging variables were averaged or summed per individal
- nightly_params: Same as "Pair_params" but movement/foraging variables were averaged by night.
- pamm_table: Table arranged to run time-to-event model using PAMMs
- chicks_table: Table containing measuraments of chicks (several measuraments per chicks) annotated with their parents' provisioning share and other movement/foraging parameters.

## Scripts (.R)
- 1.PredFlexBipCare: Script relative to the Result sub-sections: "Variation in biparental provisioning ", "Within-night adjustments of parental effort" and "Between-night carry-over effects"
- 2.PAMM_chick-growth-survival: Script relative to the results on piece‐wise exponential additive mixed model (PAMM), and Result sub-section "Contributions of specific parental behaviours to nestling survival".
