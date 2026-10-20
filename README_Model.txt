WOLBACHIA ABM — CORE MODEL SCRIPTS
---------------------------------
This README describes the four Julia scripts implementing the agent-based mosquito population model. Together, they define the model parameters, individual adult mosquito and immature-cohort dynamics, reproductive interactions, and the main simulation workflow. The analysis scripts call these components to evaluate Wolbachia release strategies and establishment outcomes.

1. parameters.jl — MODEL PARAMETERS
----------------------------------
Defines the MosquitoParameters structure, including simulation settings, environmental input paths, adult survival and immature development parameters, mating and gonotrophic-cycle parameters, Wolbachia transmission, fecundity, hatching, and density-dependent recruitment. The carrying capacity (K) is supplied for each simulation or analysis rather than fixed by this file.

2. mosquito.jl — MOSQUITO STATES AND DEMOGRAPHY
----------------------------------------------
Defines female and male mosquito agents and immature cohorts. Implements population initialization, temperature-dependent immature development and adult mortality, emergence and sex allocation, and Wolbachia-infected mosquito releases. Released adults are assigned an age of one or two days, with physiological age initialized consistently using the previous day's temperature.

3. interaction.jl — MATING AND REPRODUCTION
------------------------------------------
Implements blood-meal acquisition, gonotrophic-cycle development, mating, and egg production. Offspring recruitment incorporates cross-specific fecundity and hatching, cytoplasmic incompatibility, maternal Wolbachia transmission, and density dependence. Reproductive states are updated after completed cycles.

4. main.jl — SIMULATION AND EVALUATION
-------------------------------------
Loads temperature inputs and coordinates the daily simulation sequence: mosquito releases, blood meals, gonotrophic-cycle progression, reproduction, adult mortality and aging, immature development, and emergence. Tracks mosquito populations and releases, and supports both new simulations and continuation from serialized burn-in states.

The analysis scripts (EffGain.jl, MosqPestTTE.jl and Pareto.jl) include main.jl and use load_burnin, posterior_parameters, and main_from_burnin to run the same core model for their independently constructed release schedules. They use the population counts returned by the simulation and shared rolling-mean and verification utilities to calculate establishment outcomes, while defining their own post-release verification windows, time-to-establishment calculations, replicate evaluations, and search procedures. The generic establishment and release-search helpers in main.jl are separate from these analysis-specific procedures.

SCRIPT ORGANIZATION
-------------------
The four scripts are intended to reside together in scripts/, with input datasets in data/. main.jl includes parameters.jl, mosquito.jl, and interaction.jl. The analysis scripts import this model functionality and manage their own release schedules, parameter configurations, replicates, and outputs. The companion README_analysis.txt documents the analysis scripts and their SLURM execution workflow.
