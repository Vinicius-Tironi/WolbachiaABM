WOLBACHIA ABM — ANALYSIS SCRIPTS
------------------
This README describes the four Julia analysis scripts: EffGain.jl, MosqPestTTE.jl, Pareto.jl, and SMD_Calibration.jl, together with their corresponding SLURM launchers. All analyses use the Zone 3 model, environmental inputs, and a calibrated burn-in state.. Each analysis uses 500 simulation replicates per fully evaluated candidate to assess establishment outcomes, release strategies, or model calibration.

COMMON DEFINITIONS
------------------
The analyses consider three mosquito-release strategies: daily, weekly, and fortnightly. The high-survival (HS) period corresponds to weeks 18–43 of the annual cycle, during which environmental conditions are expected to favor adult mosquito survival. The seasonal allocation parameter (alpha_S) determines the proportion of the total release effort allocated to this period.

Successful Wolbachia establishment is defined by the persistence of a 30-day moving average of adult Wolbachia prevalence at or above 0.60 throughout the verification period, which extends from 30 to 730 days after the final release. For replicates satisfying this criterion, the time to establishment (TTE) is measured from the first release to the earliest day after which the prevalence threshold is maintained continuously through the end of the verification period.

A release effort is considered to achieve reliable establishment when the estimated establishment probability (P_est) is at least 0.90. Establishment is further classified according to its speed: fast establishment occurs when the rounded mean TTE among successfully established replicates does not exceed 365 days for alpha_S < 2.0 or 182 days for alpha_S = 2.0, whereas slow establishment refers to cases satisfying the establishment probability criterion but exceeding the corresponding TTE threshold.

These definitions are applied across the analyses to distinguish establishment reliability from establishment speed. In particular, EffGain identifies the minimum release effort (M*) satisfying both criteria, corresponding to the first fast-establishment candidate identified by its search procedure, including downward verification. MosqPestTTE characterizes the transition between unsuccessful, slow, and fast establishment across the evaluated release-effort grid, while Pareto determines the minimum release effort based exclusively on the establishment probability criterion (P_est ≥ 0.90), after which the resulting establishment outcomes are classified as slow or fast according to the corresponding mean time to establishment.

1. EffGain.jl — MINIMUM RELEASE EFFORT
--------------------------------------

The EffGain analysis determines the minimum total release effort (M*) required to achieve fast Wolbachia establishment for each release strategy (daily, weekly, and fortnightly) and seasonal allocation parameter alpha_S ∈ {1.0, 1.25, 1.50, 1.75, 2.0}. The allocation of released mosquitoes to the high-survival period follows M_HS ≈ (alpha_S / 2) M, with integer adjustments to ensure the prescribed total release effort.

The search procedure first identifies the establishment threshold using coarse increments of 5,000 mosquitoes, followed by a refined search with increments of 1,000 mosquitoes. An additional downward verification examines every 1,000-mosquito increment within a 20,000-mosquito interval below the provisional M*, allowing the identification of qualifying release efforts that may have been missed during the initial search. Each candidate is evaluated using 500 independent simulation replicates, and the final M* must satisfy both the establishment probability and fast-establishment criteria.

The corresponding efficiency gain is calculated as E_g = 100 × (1 − M*/630,000), where 630,000 mosquitoes represents the fixed reference release effort used for comparison rather than an upper bound on the search. For the manuscript figures, the release effort obtained under uniform seasonal allocation (alpha_S = 1.0) is additionally used as a strategy-specific reference, M_ref = M*(alpha_S = 1.0), to express the relative reduction in mosquito releases achieved through seasonal allocation.

2. MosqPestTTE.jl — ESTABLISHMENT AND TIMING GRID
-------------------------------------------------

The MosqPestTTE analysis evaluates Wolbachia establishment probability and time to establishment (TTE) across a range of total release efforts for each release strategy (daily, weekly, and fortnightly) and seasonal allocation parameter alpha_S ∈ {1.0, 1.75, 2.0}. Starting from predefined release-effort values, the analysis systematically evaluates candidates at increments of 1,000 mosquitoes, with 500 independent simulation replicates per candidate.

For each release effort, establishment probability and TTE statistics are calculated to characterize the transition between unsuccessful, slow, and fast establishment. The grid advances until five consecutive candidates satisfy both the establishment probability and fast-establishment criteria, providing a consistent stopping rule across the investigated strategies and seasonal allocations.

3. Pareto.jl — RELEASE EVENTS AND EFFORT
----------------------------------------

The Pareto analysis investigates the relationship between the number of release events, total mosquito release effort, and time to Wolbachia establishment under high-survival (HS) allocation (alpha_S = 2.0). For each release strategy (daily, weekly, and fortnightly), the analysis determines the minimum total release effort M*(k) required to achieve an establishment probability of at least 0.90 for each number of release events (k). The procedure begins with the full HS release campaign and progressively evaluates shorter campaigns, using the minimum effort identified for the preceding configuration as the starting point for the next search.

The search procedure initially brackets the establishment threshold using coarse increments of 5,000 mosquitoes, followed by refinement at 1,000-mosquito increments. An additional downward verification evaluates every 1,000-mosquito increment within a 20,000-mosquito interval below the provisional M*(k), allowing the identification of qualifying release efforts within this range.

To improve computational efficiency, the analysis implements an exact early-classification procedure based on a fixed budget of 500 simulation replicates per candidate. A candidate is classified as successful upon reaching 450 established replicates, whereas 51 unsuccessful replicates are sufficient for rejection. These thresholds guarantee that the resulting classification is identical to that obtained by completing all 500 replicates, irrespective of the remaining outcomes. Every selected M*(k) is subsequently evaluated using the full 500 replicates to obtain establishment probability and TTE statistics. The resulting establishment outcomes are then classified as slow or fast according to the corresponding mean time to establishment.


4. SMD_Calibration.jl — CARRYING-CAPACITY CALIBRATION
-----------------------------------------------------

The SMD_Calibration analysis estimates the mosquito population carrying capacity (K) by comparing simulated monthly Wolbachia prevalence with field monitoring data from Zone 3 of Niterói, using the historical mosquito-release schedule. Model calibration is performed using a weighted squared-distance objective (Q) based on the Anscombe transformation of observed and simulated prevalence, accounting for variation in monthly sample sizes.

For each candidate carrying capacity, 500 independent simulation replicates are performed, and the calibration objective is calculated for each replicate. The selected carrying capacity corresponds to the value of K that minimizes the mean calibration score across replicates. The analysis exports replicate-level calibration results, monthly prevalence trajectories, and summary statistics, together with a burn-in snapshot corresponding to the selected carrying capacity.


CHECKPOINTS AND RESTARTS
------------------------
EffGain, MosqPestTTE, and Pareto store replicate and candidate records plus a checkpoint.jls in separate output directories for each release strategy. An interrupted job can reuse completed simulations after checking configuration compatibility; a first run starts with empty histories. Pareto also preserves early-classified candidates. SMD Calibration instead resumes by validating complete per-K replicate and trajectory CSV pairs; an interrupted K may need to be computed again.


SLURM EXECUTION
---------------
Each analysis has a matching .slurm launcher. EffGain, MosqPestTTE, and Pareto use a three-task SLURM array: task 0 runs daily, task 1 weekly, and task 2 fortnightly. Each task uses its own strategy output directory and a file lock to prevent two jobs from writing the same results simultaneously. SMD Calibration is a single job that processes its K grid, distributing replicate evaluations across Julia workers.

All four SLURM launchers request one exclusive node with 32 CPUs, use Julia 1.12.6 with 31 workers by default, and specify a seven-day wall-time limit. The number of Julia workers can be adjusted through the corresponding environment variables (EFFFREQ_WORKERS, EFFFREQ_B_WORKERS, and PARETO_WORKERS) for the release-strategy analyses.

The model and analysis Julia files are located in scripts/, with their corresponding SLURM launchers in slurm/. Environmental inputs are stored in data/, and the required burn-in snapshots must be available in data/BurnIn/zone3/. Simulation results are organized in outputs/ by analysis and release strategy.
