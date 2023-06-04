# Estimating the Spatial Dynamics of Plant Recruitment using Approximate Bayesian Computation
## Supporting Data for MSc Thesis
**Author:** Jennifer Freeman 
**Supervisor:** Ben Bolker
**Institution:** McMaster University, Hamilton, ON, Canada
**Date:** June 3, 2023

---

Directory structure and description of files:

`ABC_Slash_Pine`
 - `┣ data`
    - `┣ LingsUTM-3.txt` seedling data from original source
    - `┣ pine2.RData` seed and seedling data used in all analysis (transformed coordinates)
    - `┗ SeedsUTM-2.txt` seed data from original source
 - `┣ estimates`
    - `┣ abc.R` ABC estimates
    - `┗ mle.R` MLE seed estimates
 - `┣ general`
    - `┣ abc_local.R` local ABC function modified from `abc::abc()` to remove MAD scaling
    - `┣ functions.R` general functions for priors, data simulations and data prep
    - `┗ R_sessionInfo.txt` output from R `sessionInfo()` (all required packages with version numbers)
 - `┣ priors`
   - `┣ priors.html` output from `priors.Rmd`
   - `┗ priors.Rmd` prior visualization, PPD for mean seed count and mean establishment probability, exploring reparameterized Matérn scale
 - `┣ simulations`
    - `┣ farm` META-Farm shell scripts and case file for simulations, run in SIMPLE mode. No modifications were made to `single_case.sh`.
      - `┣ final.sh`
      - `┣ job_script.sh`
      - `┣ resubmit_script.sh`
      - `┗ table.dat`
    - `┣ output`
      - `┣ params_meta_splinecorr.RDS` 1e6 prior parameter draws from META-Farm run using files in `farm`, output from `simulate.R` and `combine.R`
      - `┗ summary_statistics_meta_splinecorr.RDS` 1e6 summary statistics from META-Farm run, output from `simulate.R` and `combine.R`. 
    - `┣ combine.R` post-processing script using META-Farm package to aggregate results
    - `┣ r_parallel.R` exploratory script to investigate R parallel packages to generate simulations
    - `┗ simulate.R` script to generate simulations using META-Farm package
 - `┣ summary_statistics`
    - `┣ output`
      - `┗ observed_stats.RDS` output from `observed_sumstats.R`
   - `┣ observed_sumstats.R` script to compute summary statistics of observed data
   - `┗ spline_correlogram.R` exploratory analysis of spline correlograms, setting degrees of freedom and testing distance breaks
 - `┣ thesis_plots`
    - `┣ output` TikZ thesis plot images, output from `thesis_plots.R`
      - `┣ abc_parameters.tex`
      - `┣ coverage.tex`
      - `┣ cv_quant.tex`
      - `┣ dataset.tex`
      - `┣ grf_shape.tex`
      - `┣ grf_shape_ras1.png`
      - `┣ grf_shape_ras2.png`
      - `┣ grf_shape_ras3.png`
      - `┣ grf_shape_ras4.png`
      - `┣ interpoint_distances.tex`
      - `┣ matern_reparam.tex`
      - `┣ profile_likelihood.tex`
      - `┣ profile_likelihood_nugget.tex`
      - `┣ rmsre.tex`
      - `┣ sbc.tex`
      - `┗ splinecorr.tex`
    - `┗ thesis_plots.R` - generates all plots in thesis
- `┣ validation`
  - `┣ farm` META-Farm shell scripts and case file for SBC validation, run in SIMPLE mode. No modifications were made to `single_case.sh`.
    - `┣ job_script.sh`
    - `┣ resubmit_script.sh`
    - `┗ table.dat`
  - `┣ output`
    - `┣ sbc` output from `sbc.R` and `sbc_plot.R` with MAD scaling
    - `┣ sbc_ac` output from `sbc.R` and `sbc_plot.R` with 0.5% acceptance rate
    - `┣ sbc_noscale` output from `sbc.R` and `sbc_plot.R` with no scaling
    - `┣ sbc_scale.default` output from `sbc.R` and `sbc_plot.R` with `scale.default()` (subtract mean, divide by standard deviation)
    - `┣ sbc_scaleonly` output from `sbc.R` and `sbc_plot.R` with scaling by standard deviation
    - `┣ cat_plot.RDS` output from `coverage.R`
    - `┣ check_coverage.RDS` output from `coverage.R`
    - `┣ knownparam_stats_splinecorr.RDS` output from `known_parameter.R`
    - `┗ rmsre_results.RDS` output from `rmsre.R`
  - `┣ coverage.R` percent coverage
  - `┣ known_parameter.R` ABC run using known parameter set
  - `┣ rmsre.R` validating acceptance rate using RMSRE
  - `┣ sbc.R` script to generate SBC results using META-Farm package
  - `┗ sbc_plot.R` plots of all SBC results
