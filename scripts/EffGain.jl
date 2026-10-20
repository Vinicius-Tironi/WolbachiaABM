using Distributed
using CSV
using Serialization
using SHA
using Statistics
using Printf

const MODEL_ROOT = dirname(@__DIR__)
const PROJECT_ROOT = MODEL_ROOT
const MODEL_ENTRYPOINT = joinpath(@__DIR__, "main.jl")

@everywhere include($MODEL_ENTRYPOINT)

@everywhere const EFFFREQ_WORKER_CONTEXT = Ref{Any}(nothing)

@everywhere function initialize_efffreq_worker!(P, burnin)
    EFFFREQ_WORKER_CONTEXT[] = (P=P, burnin=burnin)
    return nothing
end

# EffGain: minimum release effort and efficiency gain experiment (fig1) 
# R=500 replicates per candidate, for each strategy and alpha.
# Seasonal allocation: alpha_S = 1.0, 1.25, 1.50, 1.75, 2.0
# Candidate release allocation: M_HS = alpha_S/2 * M; M_outside = M - M_HS
# Establishment: 30-day moving prevalence >= 0.60 from last release +30 through +730 days.
# Search uses 5,000-mosquito coarse and 1,000-mosquito fine steps, followed by a 20,000-mosquito downward verification; M_ref = 630,000 is a reference value, not a search ceiling
# For each strategy and alpha_S, search total release effort M, M*(alpha_S) = minimum M such that both:
#(i) P_est(M) >= 0.90; (ii) mean TTE among persistently established replicates <= T_max(alpha_S), T_max = {365 days, alpha<2; 182 days, alpha=2}
# Afterwards, efficiency_gain = 100.0 * (1 - M_star / M_ref)
# Finally, the selected M* must satisfy both the establishment (P_est >= 0.90) and mean-TTE (rounded mean TTE <= T_max) criteria. A downward verification
# evaluates all candidates at 1,000-mosquito increments within 20,000 mosquitoes below the provisional M*, retaining the smallest qualifying effort. The
# immediate predecessor (M* - 1,000) is also explicitly verified to be nonqualifying.


# Release design

const VALID_STRATEGIES = [:daily, :weekly, :fortnightly]
const ALPHAS = Float64[1.0, 1.25, 1.50, 1.75, 2.0]

const REFERENCE_TOTAL_M = 630_000
const MIN_TOTAL_M = 1_000

const RELEASE_FEMALE_PROBABILITY = 0.5

const FAVORABLE_FIRST_WEEK = 18
const FAVORABLE_LAST_WEEK = 43
const FAVORABLE_FIRST_DAY = 1 + 7 * (FAVORABLE_FIRST_WEEK - 1)
const FAVORABLE_LAST_DAY = 7 * FAVORABLE_LAST_WEEK

const STRATEGY_SPECS = Dict(
    :daily => (label="Daily", interval=1, events=365, favorable_events=182),
    :weekly => (label="Weekly", interval=7, events=52, favorable_events=26),
    :fortnightly => (label="Fortnightly", interval=14, events=26, favorable_events=13),
)

length(ARGS) == 1 || error("Usage: julia -p <workers> EffGain.jl <daily|weekly|fortnightly>")

const STRATEGY = Symbol(lowercase(strip(ARGS[1])))
STRATEGY in VALID_STRATEGIES || error("Invalid strategy '$STRATEGY'. Choose daily, weekly, or fortnightly.")

# Replicate / search settings

const NUMBER_SIMULATIONS = 500
const FIRST_SIMULATIONNUMBER = 1

const REQUIRED_SUCCESS_PROBABILITY = 0.90
const REQUIRED_SUCCESSES = ceil(Int64, REQUIRED_SUCCESS_PROBABILITY * NUMBER_SIMULATIONS - 1e-12)

const COARSE_STEP = 5_000
const FINE_STEP = 1_000
const DOWNWARD_VERIFY_WINDOW = 20_000

const PARALLEL_BATCH_SIZE = max(1, min(NUMBER_SIMULATIONS, max(nworkers(), 1)))

# Establishment settings

const ESTABLISHMENT_THRESHOLD = 0.60
const ROLLING_WINDOW = 30
const WASHOUT_DAYS = 30
const POST_RELEASE_FOLLOWUP_DAYS = 730
const SIMULATION_END_DAY = 1095

tte_deadline_days(alpha::Real) = isapprox(Float64(alpha), 2.0; atol=1e-12, rtol=0.0) ? Int64(182) : Int64(365)

# Model settings
const ZONE = 3
const K_VALUE = 25000.0

const BURNIN_DIR = joinpath(MODEL_ROOT, "data", "BurnIn")
const ENVIRONMENTAL_DATA = joinpath(MODEL_ROOT, "data", "ambiental.csv")

# Output / checkpoint
const RUN_ID = "EffGain_R500"

const OUTPUT_ROOT = get(ENV, "EFFGAIN_OUTPUT_ROOT", joinpath(MODEL_ROOT, "outputs", "EffGain"))

const OUTPUT_DIR = joinpath(OUTPUT_ROOT, string(STRATEGY))
const REPLICATE_CSV = joinpath(OUTPUT_DIR, "replicates.csv")
const CANDIDATE_CSV = joinpath(OUTPUT_DIR, "candidates.csv")
const REQUIRED_EFFORT_CSV = joinpath(OUTPUT_DIR, "Mstar.csv")
const CONFIGURATION_CSV = joinpath(OUTPUT_DIR, "config.csv")
const CHECKPOINT_JLS = joinpath(OUTPUT_DIR, "checkpoint.jls")

mkpath(OUTPUT_DIR)

# Basic helpers
percentage_string(x::Real) = isnan(Float64(x)) ? "NA" : "$(round(100.0 * Float64(x), digits=2))%"

csv_value(x) = ismissing(x) || isnothing(x) ? "NA" : string(x)

alpha_tag(alpha::Real) = replace(@sprintf("%.2f", Float64(alpha)), "." => "p")
alpha_code(alpha::Real) = round(Int64, 100 * Float64(alpha))

function mean_or_missing(values)
    valid = Float64[Float64(x) for x in values if !ismissing(x) && !isnothing(x) && !isnan(Float64(x))]
    return isempty(valid) ? missing : mean(valid)
end

function median_or_missing(values)
    valid = Float64[Float64(x) for x in values if !ismissing(x) && !isnothing(x) && !isnan(Float64(x))]
    return isempty(valid) ? missing : median(valid)
end

function successful_tte_stats(records)
    values = Float64[Float64(r.time_to_establishment) for r in records if r.established && !ismissing(r.time_to_establishment)]

    isempty(values) && return (
        n=0,
        mean=missing,
        sd=missing,
        se=missing,
        lower95=missing,
        upper95=missing,
        q025=missing,
        median=missing,
        q975=missing,
   )

    n = length(values)
    mu = mean(values)
    sd_value = n > 1 ? std(values) : 0.0
    se_value = sd_value / sqrt(n)
    margin = 1.96 * se_value

    return (
        n=n,
        mean=mu,
        sd=sd_value,
        se=se_value,
        lower95=mu - margin,
        upper95=mu + margin,
        q025=quantile(values, 0.025),
        median=median(values),
        q975=quantile(values, 0.975),
   )
end

function split_integer_total(total::Int, n::Int)
    total >= 0 || error("total must be nonnegative")
    n >= 1 || error("n must be at least 1")

    base, remainder = divrem(total, n)
    values = fill(Int64(base), n)

    for i in 1:remainder
        values[i] += 1
    end

    return values
end

round_up_to_step(value::Int64, step::Int64) = value <= 0 ? step : Int64(cld(value, step) * step)

# Release schedule
function release_days(strategy::Symbol)
    spec = STRATEGY_SPECS[strategy]
    return Int64[1 + spec.interval * i for i in 0:(spec.events - 1)]
end

is_hs_day(day::Int) = FAVORABLE_FIRST_DAY <= day <= FAVORABLE_LAST_DAY

# Allocate the total release effort between HS and non-HS periods.
function alpha_release_schedule(strategy::Symbol, alpha::Float64, total_M::Int64)
    0.0 <= alpha <= 2.0 || error("alpha must be between 0 and 2")
    total_M >= MIN_TOTAL_M || error("total_M=$total_M is below MIN_TOTAL_M=$MIN_TOTAL_M")

    days = release_days(strategy)
    hs = findall(day -> is_hs_day(day), days)
    outside = findall(day -> !is_hs_day(day), days)

    length(hs) == STRATEGY_SPECS[strategy].favorable_events || error("Unexpected number of HS events for $strategy")

    hs_total = isapprox(alpha, 2.0; atol=1e-12, rtol=0.0) ? total_M : round(Int64, total_M * alpha / 2.0)
    outside_total = total_M - hs_total

    sizes = zeros(Int64, length(days))
    sizes[hs] = split_integer_total(hs_total, length(hs))
    sizes[outside] = split_integer_total(outside_total, length(outside))

    schedule = Dict(days[i] => sizes[i] for i in eachindex(days))

    return schedule
end

function release_design_metrics(strategy::Symbol, alpha::Float64, total_M::Int64)
    schedule = alpha_release_schedule(strategy, alpha, total_M)
    days = release_days(strategy)
    positive_days = sort(Int64[day for day in days if schedule[day] > 0])

    hs_total = sum(schedule[day] for day in days if is_hs_day(day))
    outside_total = total_M - hs_total

    first_release_day = first(positive_days)
    last_release_day = last(positive_days)
    positive_release_events = length(positive_days)

    return (
        release_schedule=schedule,
        total_M=total_M,
        hs_total=hs_total,
        outside_total=outside_total,
        hs_fraction=hs_total / total_M,
        outside_fraction=outside_total / total_M,
        scheduled_release_events=length(days),
        positive_release_events=positive_release_events,
        release_interval_days=STRATEGY_SPECS[strategy].interval,
        first_release_day=first_release_day,
        last_release_day=last_release_day,
        campaign_duration_days=last_release_day - first_release_day,
        mosquitoes_per_positive_release=total_M / positive_release_events,
   )
end

# Load model / burn-in
burnin, burnin_path = load_burnin(BURNIN_DIR, Int64(ZONE), Float64(K_VALUE))

P = posterior_parameters(burnin, Int64(SIMULATION_END_DAY), false, ENVIRONMENTAL_DATA)

# Worker-side metrics
@everywhere function efffreq_safe_prevalence(infected, total)
    return Float64[total[i] > 0 ? infected[i] / total[i] : 0.0 for i in eachindex(total)]
end

@everywhere function efffreq_persistent_day(days, rolling, threshold::Float64, first_day::Int64, last_day::Int64)
    first_index = findfirst(day -> day >= first_day, days)
    last_index = findlast(day -> day <= last_day, days)

    isnothing(first_index) && return missing
    isnothing(last_index) && return missing

    for i in first_index:last_index
        if !isnan(rolling[i]) &&
           rolling[i] >= threshold &&
           all(!isnan(rolling[j]) && rolling[j] >= threshold for j in i:last_index)
            return Int64(days[i])
        end
    end

    return missing
end

# For established replicates, record the first day from which the prevalence threshold is maintained through the verification endpoint.
@everywhere function execute_efffreq_replicate(
    simulationnumber::Int64,
    release_schedule,
    threshold::Float64,
    rolling_window::Int64,
    washout_days::Int64,
    post_release_followup_days::Int64,
    simulation_end_day::Int64,
)
    context = EFFFREQ_WORKER_CONTEXT[]
    isnothing(context) && error("worker context is not initialized, check slurm")

    result = main_from_burnin(simulationnumber, context.P, context.burnin; release_schedule=release_schedule)

    positive_release_days = sort(Int64[day for (day, size) in release_schedule if size > 0])

    first_release_day = first(positive_release_days)
    last_release_day = last(positive_release_days)

    state_indices = findall(day -> 1 <= day <= simulation_end_day, result.state_calendar_days)

    days = Int64.(result.state_calendar_days[state_indices])

    adult_free = Int64.(result.free_ctr[state_indices])
    adult_wolb = Int64.(result.wolb_ctr[state_indices])
    adult_total = adult_free .+ adult_wolb

    adult_prevalence = efffreq_safe_prevalence(adult_wolb, adult_total)
    rolling_prevalence = Float64.(rolling_mean(adult_prevalence, rolling_window))

    verification_start_day = last_release_day + washout_days
    verification_end_day = last_release_day + post_release_followup_days

    verification_end_day <= simulation_end_day || error("Verification endpoint $verification_end_day exceeds simulation endpoint " * "$simulation_end_day for seed=$simulationnumber")
    # TTE is measured from the first release, rather than from the end of the release campaign.
    established = period_stable(days, rolling_prevalence, verification_start_day, verification_end_day, threshold)

    verification_minimum = period_minimum(days, rolling_prevalence, verification_start_day, verification_end_day)

    establishment_day = established ? efffreq_persistent_day(days, rolling_prevalence, threshold, first_release_day, verification_end_day) : missing

    time_to_establishment = ismissing(establishment_day) ? missing : establishment_day - first_release_day

    establishment_lag_after_campaign = ismissing(establishment_day) ? missing : establishment_day - last_release_day

    realized_released_females = sum(result.released_female_ctr)
    realized_released_males = sum(result.released_male_ctr)
    realized_released_total = realized_released_females + realized_released_males

    realized_released_female_fraction = realized_released_total > 0 ? realized_released_females / realized_released_total : NaN

    return (
        established=established,
        verification_start_day=verification_start_day,
        verification_end_day=verification_end_day,
        verification_minimum=verification_minimum,
        establishment_day=establishment_day,
        time_to_establishment=time_to_establishment,
        establishment_lag_after_campaign=establishment_lag_after_campaign,
        realized_released_females=realized_released_females,
        realized_released_males=realized_released_males,
        realized_released_total=realized_released_total,
        realized_released_female_fraction=realized_released_female_fraction,
   )
end

# Input fingerprints / checkpoint signature
function file_sha256(path::String)
    return open(path, "r") do io
        bytes2hex(sha256(io))
    end
end

optional_file_sha256(path::String) = isfile(path) ? file_sha256(path) : "not_present"

const SCRIPT_PATH = abspath(@__FILE__)

const INPUT_FINGERPRINTS = (
    script=file_sha256(SCRIPT_PATH),
    main=file_sha256(joinpath(@__DIR__, "main.jl")),
    parameters=file_sha256(joinpath(@__DIR__, "parameters.jl")),
    mosquito=file_sha256(joinpath(@__DIR__, "mosquito.jl")),
    interaction=file_sha256(joinpath(@__DIR__, "interaction.jl")),
    project=optional_file_sha256(joinpath(dirname(PROJECT_ROOT), "Project.toml")),
    manifest=optional_file_sha256(joinpath(dirname(PROJECT_ROOT), "Manifest.toml")),
    environmental_data=file_sha256(ENVIRONMENTAL_DATA),
    burnin=file_sha256(burnin_path),
)

const CHECKPOINT_SIGNATURE = (
    run_id=RUN_ID,
    julia_version=string(VERSION),
    strategy=STRATEGY,
    alphas=Tuple(ALPHAS),
    reference_total_M=REFERENCE_TOTAL_M,
    min_total_M=MIN_TOTAL_M,
    release_female_probability=RELEASE_FEMALE_PROBABILITY,
    number_simulations=NUMBER_SIMULATIONS,
    first_simulationnumber=FIRST_SIMULATIONNUMBER,
    required_success_probability=REQUIRED_SUCCESS_PROBABILITY,
    required_successes=REQUIRED_SUCCESSES,
    qualification_metric=:persistent_probability_and_rounded_mean_tte,
    tte_rounding=:julia_round_nearest_integer,
    tte_deadline=:nominal_intervention_days,
    coarse_step=COARSE_STEP,
    fine_step=FINE_STEP,
    downward_verify_window=DOWNWARD_VERIFY_WINDOW,
    establishment_threshold=ESTABLISHMENT_THRESHOLD,
    rolling_window=ROLLING_WINDOW,
    washout_days=WASHOUT_DAYS,
    post_release_followup_days=POST_RELEASE_FOLLOWUP_DAYS,
    simulation_end_day=SIMULATION_END_DAY,
    favorable_first_week=FAVORABLE_FIRST_WEEK,
    favorable_last_week=FAVORABLE_LAST_WEEK,
    zone=ZONE,
    K=K_VALUE,
    input_fingerprints=INPUT_FINGERPRINTS,
)

# Changes to this analysis script are permitted when resuming a run.
# The model, parameters, environmental data, and burn-in fingerprints must still match the saved checkpoint.
function checkpoint_input_fingerprints_compatible(saved, current)
    for field in propertynames(current)
        # The search script fingerprint is intentionally allowed to differ so that we can resume an existing checkpoint
        # All model, parameter, data, burn-in, Project, and Manifest fingerprints remain when resumed.
        # Once the checkpoint is saved again, it receives the current script fingerprint.
        field == :script && continue

        hasproperty(saved, field) || return false
        getproperty(saved, field) == getproperty(current, field) || return false
    end

    return true
end

# Before resuming a previous run, check that the saved settings and model inputs are compatible with the current analysis.
# This prevents simulation results from different configurations from being combined.
function checkpoint_signature_compatible(saved, current)
    for field in propertynames(current)
        hasproperty(saved, field) || return false

        if field == :input_fingerprints
            checkpoint_input_fingerprints_compatible(getproperty(saved, field), getproperty(current, field)) || return false
        elseif field == :alphas
            saved_alphas = Set(Float64.(getproperty(saved, field)))
            current_alphas = Set(Float64.(getproperty(current, field)))
            issubset(saved_alphas, current_alphas) || return false
        else
            getproperty(saved, field) == getproperty(current, field) || return false
        end
    end

    return true
end

# CSV fields
# Note: The current outputs contain an exhaustive list of metrics and quantities used in the analysis.
# The final version will be filtered to retain only the relevant data.
const REPLICATE_FIELDS = Symbol[
    :run_id,
    :strategy,
    :alpha,
    :total_M,
    :reference_M,
    :hs_total,
    :outside_total,
    :hs_fraction,
    :outside_fraction,
    :scheduled_release_events,
    :positive_release_events,
    :release_interval_days,
    :first_release_day,
    :last_release_day,
    :campaign_duration_days,
    :mosquitoes_per_positive_release,
    :simulationnumber,
    :release_female_probability,
    :realized_released_females,
    :realized_released_males,
    :realized_released_total,
    :realized_released_female_fraction,
    :established,
    :establishment_threshold,
    :rolling_window_days,
    :washout_days,
    :post_release_followup_days,
    :verification_start_day,
    :verification_end_day,
    :verification_minimum,
    :establishment_day,
    :time_to_establishment,
    :establishment_lag_after_campaign,
]

const CANDIDATE_FIELDS = Symbol[
    :run_id,
    :evaluation_order,
    :strategy,
    :alpha,
    :total_M,
    :reference_M,
    :hs_total,
    :outside_total,
    :hs_fraction,
    :outside_fraction,
    :scheduled_release_events,
    :positive_release_events,
    :release_interval_days,
    :first_release_day,
    :last_release_day,
    :campaign_duration_days,
    :mosquitoes_per_positive_release,
    :replicates_completed,
    :persistent_successes,
    :persistent_failures,
    :persistent_establishment_probability,
    :persistent_criterion_met,
    :mean_tte_criterion_met,
    :qualifies,
    :required_success_probability,
    :required_successes,
    :qualification_metric,
    :tte_deadline_days,
    :mean_verification_minimum,
    :successful_time_replicates,
    :mean_time_to_establishment,
    :rounded_mean_time_to_establishment,
    :std_time_to_establishment,
    :se_time_to_establishment,
    :lower_95ci_time_to_establishment,
    :upper_95ci_time_to_establishment,
    :q025_time_to_establishment,
    :median_time_to_establishment,
    :q975_time_to_establishment,
    :mean_establishment_day,
    :median_establishment_day,
    :verification_start_day,
    :verification_end_day,
]

const REQUIRED_EFFORT_FIELDS = Symbol[
    :run_id,
    :strategy,
    :alpha,
    :status,
    :reference_M,
    :starting_anchor_M,
    :starting_anchor_source,
    :M_star,
    :efficiency_gain_percent,
    :hs_total_at_M_star,
    :outside_total_at_M_star,
    :hs_fraction_at_M_star,
    :outside_fraction_at_M_star,
    :persistent_successes_at_M_star,
    :persistent_failures_at_M_star,
    :persistent_establishment_probability_at_M_star,
    :persistent_criterion_met_at_M_star,
    :mean_tte_criterion_met_at_M_star,
    :number_simulations,
    :tte_deadline_days_at_M_star,
    :mean_verification_minimum_at_M_star,
    :successful_time_replicates_at_M_star,
    :mean_time_to_establishment_at_M_star,
    :rounded_mean_time_to_establishment_at_M_star,
    :std_time_to_establishment_at_M_star,
    :se_time_to_establishment_at_M_star,
    :lower_95ci_time_to_establishment_at_M_star,
    :upper_95ci_time_to_establishment_at_M_star,
    :q025_time_to_establishment_at_M_star,
    :median_time_to_establishment_at_M_star,
    :q975_time_to_establishment_at_M_star,
    :mean_establishment_day_at_M_star,
    :median_establishment_day_at_M_star,
    :verification_start_day_at_M_star,
    :verification_end_day_at_M_star,
    :lower_failed_M,
    :lower_failed_persistent_successes,
    :lower_failed_persistent_establishment_probability,
    :lower_failed_persistent_criterion_met,
    :lower_failed_mean_time_to_establishment,
    :lower_failed_rounded_mean_time_to_establishment,
    :lower_failed_mean_tte_criterion_met,
]

# Checkpoint state
# Allows to resume only if the settings are compatible
if isfile(CHECKPOINT_JLS)
    checkpoint = deserialize(CHECKPOINT_JLS)

    checkpoint_signature_compatible(checkpoint.signature, CHECKPOINT_SIGNATURE) || error("Existing checkpoint uses different fixed settings or input files. " * "Use a new output directory or remove the old checkpoint.")

    replicate_history = checkpoint.replicate_history
    candidate_history = checkpoint.candidate_history
    required_effort_results = checkpoint.required_effort_results
else
    replicate_history = NamedTuple[]
    candidate_history = NamedTuple[]
    required_effort_results = NamedTuple[]
end

const replicate_cache = Dict{Tuple{Symbol,Int64,Int64,Int64},NamedTuple}()

for r in replicate_history
    replicate_cache[(r.strategy, alpha_code(r.alpha), r.total_M, r.simulationnumber)] = r
end

# CSV / checkpoint writing

function write_namedtuple_csv(path::String, fields::Vector{Symbol}, records)
    mkpath(dirname(path))
    temporary_path = path * ".tmp"

    open(temporary_path, "w") do io
        println(io, join(string.(fields), ","))

        for r in records
            println(io, join([csv_value(getproperty(r, field)) for field in fields], ","))
        end
    end

    mv(temporary_path, path; force=true)
    return nothing
end

function write_csvs()
    write_namedtuple_csv(REPLICATE_CSV, REPLICATE_FIELDS, replicate_history)
    write_namedtuple_csv(CANDIDATE_CSV, CANDIDATE_FIELDS, candidate_history)
    write_namedtuple_csv(REQUIRED_EFFORT_CSV, REQUIRED_EFFORT_FIELDS, required_effort_results)
    return nothing
end

function save_checkpoint(; write_output_csvs::Bool=false)
    temporary_path = CHECKPOINT_JLS * ".tmp"

    serialize(temporary_path, (signature=CHECKPOINT_SIGNATURE, replicate_history=replicate_history, candidate_history=candidate_history, required_effort_results=required_effort_results))

    mv(temporary_path, CHECKPOINT_JLS; force=true)
    write_output_csvs && write_csvs()
    return nothing
end

function write_configuration_csv()
    row = (
        run_id=RUN_ID,
        strategy=string(STRATEGY),
        release_female_probability=RELEASE_FEMALE_PROBABILITY,
        alpha_first=first(ALPHAS),
        alpha_last=last(ALPHAS),
        alpha_values=join(ALPHAS, ";"),
        number_alpha_values=length(ALPHAS),
        reference_M=REFERENCE_TOTAL_M,
        minimum_M=MIN_TOTAL_M,
        coarse_step=COARSE_STEP,
        fine_step=FINE_STEP,
        downward_verify_window=DOWNWARD_VERIFY_WINDOW,
        number_simulations=NUMBER_SIMULATIONS,
        first_simulationnumber=FIRST_SIMULATIONNUMBER,
        required_success_probability=REQUIRED_SUCCESS_PROBABILITY,
        required_successes=REQUIRED_SUCCESSES,
        qualification_metric="persistent_probability_and_rounded_mean_tte",
        tte_rounding="julia_round_nearest_integer",
        tte_deadline="365_if_alpha_lt_2__182_if_alpha_eq_2",
        establishment_threshold=ESTABLISHMENT_THRESHOLD,
        rolling_window_days=ROLLING_WINDOW,
        washout_days=WASHOUT_DAYS,
        post_release_followup_days=POST_RELEASE_FOLLOWUP_DAYS,
        simulation_end_day=SIMULATION_END_DAY,
        favorable_first_week=FAVORABLE_FIRST_WEEK,
        favorable_last_week=FAVORABLE_LAST_WEEK,
        favorable_first_day=FAVORABLE_FIRST_DAY,
        favorable_last_day=FAVORABLE_LAST_DAY,
        zone=ZONE,
        K=K_VALUE,
        burnin_path=burnin_path,
        environmental_data=ENVIRONMENTAL_DATA,
        julia_version=string(VERSION),
        parallel_batch_size=PARALLEL_BATCH_SIZE,
   )

    temporary_path = CONFIGURATION_CSV * ".tmp"
    CSV.write(temporary_path, [row])
    mv(temporary_path, CONFIGURATION_CSV; force=true)
    return nothing
end

# Resume lookups
function existing_candidate_record(alpha::Float64, total_M::Int64)
    matches = [r for r in candidate_history if r.strategy == STRATEGY && alpha_code(r.alpha) == alpha_code(alpha) && r.total_M == total_M && r.replicates_completed == NUMBER_SIMULATIONS]

    return isempty(matches) ? nothing : last(matches)
end

function remove_existing_candidate_record!(alpha::Float64, total_M::Int64)
    filter!(r -> !(r.strategy == STRATEGY && alpha_code(r.alpha) == alpha_code(alpha) && r.total_M == total_M), candidate_history)
    return nothing
end

function completed_result(alpha::Float64)
    matches = [r for r in required_effort_results if r.strategy == STRATEGY && alpha_code(r.alpha) == alpha_code(alpha)]

    return isempty(matches) ? nothing : last(matches)
end

# Replicate recording
function build_replicate_record(alpha::Float64, total_M::Int64, simulationnumber::Int64, design, metrics)
    return (
        run_id=RUN_ID,
        strategy=STRATEGY,
        alpha=alpha,
        total_M=total_M,
        reference_M=REFERENCE_TOTAL_M,
        hs_total=design.hs_total,
        outside_total=design.outside_total,
        hs_fraction=design.hs_fraction,
        outside_fraction=design.outside_fraction,
        scheduled_release_events=design.scheduled_release_events,
        positive_release_events=design.positive_release_events,
        release_interval_days=design.release_interval_days,
        first_release_day=design.first_release_day,
        last_release_day=design.last_release_day,
        campaign_duration_days=design.campaign_duration_days,
        mosquitoes_per_positive_release=design.mosquitoes_per_positive_release,
        simulationnumber=simulationnumber,
        release_female_probability=RELEASE_FEMALE_PROBABILITY,
        realized_released_females=metrics.realized_released_females,
        realized_released_males=metrics.realized_released_males,
        realized_released_total=metrics.realized_released_total,
        realized_released_female_fraction=metrics.realized_released_female_fraction,
        established=metrics.established,
        establishment_threshold=ESTABLISHMENT_THRESHOLD,
        rolling_window_days=ROLLING_WINDOW,
        washout_days=WASHOUT_DAYS,
        post_release_followup_days=POST_RELEASE_FOLLOWUP_DAYS,
        verification_start_day=metrics.verification_start_day,
        verification_end_day=metrics.verification_end_day,
        verification_minimum=metrics.verification_minimum,
        establishment_day=metrics.establishment_day,
        time_to_establishment=metrics.time_to_establishment,
        establishment_lag_after_campaign=metrics.establishment_lag_after_campaign,
   )
end

function store_replicate_record!(record)
    key = (
        record.strategy,
        alpha_code(record.alpha),
        record.total_M,
        record.simulationnumber,
   )

    if !haskey(replicate_cache, key)
        replicate_cache[key] = record
        push!(replicate_history, record)
    end

    return record
end

# Full R=500 candidate evaluation...
# Evaluate a release-effort candidate using 500 simulation replicates.
# Previously completed replicates are reused when available.
function evaluate_candidate(alpha::Float64, total_M::Int64)
    existing = existing_candidate_record(alpha, total_M)
    !isnothing(existing) && return existing

    design = release_design_metrics(STRATEGY, alpha, total_M)
    seeds = FIRST_SIMULATIONNUMBER:(FIRST_SIMULATIONNUMBER + NUMBER_SIMULATIONS - 1)

    records = NamedTuple[replicate_cache[(STRATEGY, alpha_code(alpha), total_M, Int64(seed))] for seed in seeds if haskey(replicate_cache, (STRATEGY, alpha_code(alpha), total_M, Int64(seed)))]
    sort!(records, by=r -> r.simulationnumber)

    !isempty(records) && println("Resuming alpha=$alpha, M=$total_M: $(length(records))/$NUMBER_SIMULATIONS replicates")

    pending_seeds = Int64[Int64(seed) for seed in seeds if !haskey(replicate_cache, (STRATEGY, alpha_code(alpha), total_M, Int64(seed)))]

    while !isempty(pending_seeds)
        batch_length = min(PARALLEL_BATCH_SIZE, length(pending_seeds))
        batch_seeds = pending_seeds[1:batch_length]
        deleteat!(pending_seeds, 1:batch_length)

        release_schedule = design.release_schedule

        run_seed = simulationnumber -> execute_efffreq_replicate(simulationnumber, release_schedule, ESTABLISHMENT_THRESHOLD, ROLLING_WINDOW, WASHOUT_DAYS, POST_RELEASE_FOLLOWUP_DAYS, SIMULATION_END_DAY)

        batch_metrics = isempty(workers()) ?
            map(run_seed, batch_seeds) :
            pmap(run_seed, WorkerPool(workers()), batch_seeds; batch_size=1)

        for (simulationnumber, metrics) in zip(batch_seeds, batch_metrics)
            record = build_replicate_record(alpha, total_M, simulationnumber, design, metrics)
            store_replicate_record!(record)
            push!(records, record)
        end

        sort!(records, by=r -> r.simulationnumber)
        save_checkpoint(write_output_csvs=false)

    end

    length(records) == NUMBER_SIMULATIONS || error("Candidate alpha=$alpha M=$total_M did not complete all $NUMBER_SIMULATIONS replicates")
    # Estimate establishment probability using all 500 replicates.
    # Candidates with P_est < 0.90 are classified as unsuccessful.
    # Among those meeting P_est >= 0.90, the rounded mean TTE determines whether establishment is slow or fast.
    persistent_successes = count(r -> r.established, records)
    persistent_failures = NUMBER_SIMULATIONS - persistent_successes
    persistent_establishment_probability = persistent_successes / NUMBER_SIMULATIONS
    persistent_criterion_met = persistent_successes >= REQUIRED_SUCCESSES

    tte = successful_tte_stats(records)
    tte_deadline = tte_deadline_days(alpha)
    rounded_mean_time_to_establishment = ismissing(tte.mean) ? missing : round(Int64, tte.mean)
    mean_tte_criterion_met = !ismissing(rounded_mean_time_to_establishment) &&
        rounded_mean_time_to_establishment <= tte_deadline

    qualifies = persistent_criterion_met && mean_tte_criterion_met
    establishment_days = Float64[Float64(r.establishment_day) for r in records if !ismissing(r.establishment_day)]

    remove_existing_candidate_record!(alpha, total_M)

    candidate_record = (
        run_id=RUN_ID,
        evaluation_order=length(candidate_history) + 1,
        strategy=STRATEGY,
        alpha=alpha,
        total_M=total_M,
        reference_M=REFERENCE_TOTAL_M,
        hs_total=design.hs_total,
        outside_total=design.outside_total,
        hs_fraction=design.hs_fraction,
        outside_fraction=design.outside_fraction,
        scheduled_release_events=design.scheduled_release_events,
        positive_release_events=design.positive_release_events,
        release_interval_days=design.release_interval_days,
        first_release_day=design.first_release_day,
        last_release_day=design.last_release_day,
        campaign_duration_days=design.campaign_duration_days,
        mosquitoes_per_positive_release=design.mosquitoes_per_positive_release,
        replicates_completed=NUMBER_SIMULATIONS,
        persistent_successes=persistent_successes,
        persistent_failures=persistent_failures,
        persistent_establishment_probability=persistent_establishment_probability,
        persistent_criterion_met=persistent_criterion_met,
        mean_tte_criterion_met=mean_tte_criterion_met,
        qualifies=qualifies,
        required_success_probability=REQUIRED_SUCCESS_PROBABILITY,
        required_successes=REQUIRED_SUCCESSES,
        qualification_metric="persistent_probability_and_rounded_mean_tte",
        tte_deadline_days=tte_deadline,
        mean_verification_minimum=mean(r.verification_minimum for r in records),
        successful_time_replicates=tte.n,
        mean_time_to_establishment=tte.mean,
        rounded_mean_time_to_establishment=rounded_mean_time_to_establishment,
        std_time_to_establishment=tte.sd,
        se_time_to_establishment=tte.se,
        lower_95ci_time_to_establishment=tte.lower95,
        upper_95ci_time_to_establishment=tte.upper95,
        q025_time_to_establishment=tte.q025,
        median_time_to_establishment=tte.median,
        q975_time_to_establishment=tte.q975,
        mean_establishment_day=isempty(establishment_days) ? missing : mean(establishment_days),
        median_establishment_day=isempty(establishment_days) ? missing : median(establishment_days),
        verification_start_day=design.last_release_day + WASHOUT_DAYS,
        verification_end_day=design.last_release_day + POST_RELEASE_FOLLOWUP_DAYS,
   )

    push!(candidate_history, candidate_record)
    save_checkpoint(write_output_csvs=true)

    println("Candidate alpha=$alpha | M=$total_M | P_est=$(percentage_string(persistent_establishment_probability)) | rounded TTE=$rounded_mean_time_to_establishment | qualifies=$qualifies")

    return candidate_record
end

# Required-effort results
function qualified_required_effort_result(alpha::Float64, starting_anchor_M::Int64, starting_anchor_source::Symbol, minimum_record, lower_failed_record)
    M_star = minimum_record.total_M
    efficiency_gain = M_star <= REFERENCE_TOTAL_M ? 100.0 * (1.0 - M_star / REFERENCE_TOTAL_M) : missing

    return (
        run_id=RUN_ID,
        strategy=STRATEGY,
        alpha=alpha,
        status=:qualified,
        reference_M=REFERENCE_TOTAL_M,
        starting_anchor_M=starting_anchor_M,
        starting_anchor_source=starting_anchor_source,
        M_star=M_star,
        efficiency_gain_percent=efficiency_gain,
        hs_total_at_M_star=minimum_record.hs_total,
        outside_total_at_M_star=minimum_record.outside_total,
        hs_fraction_at_M_star=minimum_record.hs_fraction,
        outside_fraction_at_M_star=minimum_record.outside_fraction,
        persistent_successes_at_M_star=minimum_record.persistent_successes,
        persistent_failures_at_M_star=minimum_record.persistent_failures,
        persistent_establishment_probability_at_M_star=minimum_record.persistent_establishment_probability,
        persistent_criterion_met_at_M_star=minimum_record.persistent_criterion_met,
        mean_tte_criterion_met_at_M_star=minimum_record.mean_tte_criterion_met,
        number_simulations=NUMBER_SIMULATIONS,
        tte_deadline_days_at_M_star=minimum_record.tte_deadline_days,
        mean_verification_minimum_at_M_star=minimum_record.mean_verification_minimum,
        successful_time_replicates_at_M_star=minimum_record.successful_time_replicates,
        mean_time_to_establishment_at_M_star=minimum_record.mean_time_to_establishment,
        rounded_mean_time_to_establishment_at_M_star=minimum_record.rounded_mean_time_to_establishment,
        std_time_to_establishment_at_M_star=minimum_record.std_time_to_establishment,
        se_time_to_establishment_at_M_star=minimum_record.se_time_to_establishment,
        lower_95ci_time_to_establishment_at_M_star=minimum_record.lower_95ci_time_to_establishment,
        upper_95ci_time_to_establishment_at_M_star=minimum_record.upper_95ci_time_to_establishment,
        q025_time_to_establishment_at_M_star=minimum_record.q025_time_to_establishment,
        median_time_to_establishment_at_M_star=minimum_record.median_time_to_establishment,
        q975_time_to_establishment_at_M_star=minimum_record.q975_time_to_establishment,
        mean_establishment_day_at_M_star=minimum_record.mean_establishment_day,
        median_establishment_day_at_M_star=minimum_record.median_establishment_day,
        verification_start_day_at_M_star=minimum_record.verification_start_day,
        verification_end_day_at_M_star=minimum_record.verification_end_day,
        lower_failed_M=isnothing(lower_failed_record) ? missing : lower_failed_record.total_M,
        lower_failed_persistent_successes=isnothing(lower_failed_record) ? missing : lower_failed_record.persistent_successes,
        lower_failed_persistent_establishment_probability=isnothing(lower_failed_record) ? missing : lower_failed_record.persistent_establishment_probability,
        lower_failed_persistent_criterion_met=isnothing(lower_failed_record) ? missing : lower_failed_record.persistent_criterion_met,
        lower_failed_mean_time_to_establishment=isnothing(lower_failed_record) ? missing : lower_failed_record.mean_time_to_establishment,
        lower_failed_rounded_mean_time_to_establishment=isnothing(lower_failed_record) ? missing : lower_failed_record.rounded_mean_time_to_establishment,
        lower_failed_mean_tte_criterion_met=isnothing(lower_failed_record) ? missing : lower_failed_record.mean_tte_criterion_met,
   )
end

# Verify all candidates in the 20,000-mosquito interval below the provisional minimum.
# This allows isolated qualifying candidates even when the qualification curve is nonmonotonic.
function verify_downward_minimum(alpha::Float64, provisional_record)
    provisional_record.qualifies || error("Provisional M* does not qualify at alpha=$alpha")
    selected = provisional_record
    lower_bound = max(MIN_TOTAL_M, Int64(provisional_record.total_M) - DOWNWARD_VERIFY_WINDOW)

    for candidate_M in lower_bound:FINE_STEP:(Int64(provisional_record.total_M) - FINE_STEP)
        record = evaluate_candidate(alpha, Int64(candidate_M))
        if record.qualifies && record.total_M < selected.total_M
            selected = record
        end
    end

    lower_failed_record = nothing
    if selected.total_M > MIN_TOTAL_M
        lower_failed_record = evaluate_candidate(alpha, Int64(selected.total_M) - FINE_STEP)
        lower_failed_record.qualifies && error("Immediate predecessor of verified M* qualifies: alpha=$alpha M*=$(selected.total_M)")
    end
    println("Downward verification alpha=$alpha | provisional=$(provisional_record.total_M) | selected=$(selected.total_M) | window=$DOWNWARD_VERIFY_WINDOW")
    return selected, lower_failed_record
end

# Search for the minimum qualifying release effort.
# Starting from an initial anchor, move downward if it qualifies or upward otherwise, using 5,000-mosquito steps to bracket the threshold.
# Refine the bracket at 1,000-mosquito intervals, then perform the 20,000-mosquito downward verification.
function find_minimum_M(alpha::Float64, starting_anchor_M::Int64, starting_anchor_source::Symbol)
    starting_anchor_M >= MIN_TOTAL_M || error("Starting anchor M=$starting_anchor_M is below MIN_TOTAL_M=$MIN_TOTAL_M")

    anchor = max(round_up_to_step(starting_anchor_M, FINE_STEP), MIN_TOTAL_M)

    println("Searching alpha=$alpha | anchor M=$anchor ($(starting_anchor_source))")

    evaluate_M(M::Int64) = evaluate_candidate(alpha, M)

    anchor_record = evaluate_M(anchor)

    successful_M = nothing
    failed_M = nothing
    minimum_record = nothing
    lower_failed_record = nothing

    # Anchor passes -> search downward

    if anchor_record.qualifies
        successful_M = anchor

        if anchor == MIN_TOTAL_M
            minimum_record = anchor_record
        else
            candidate_M = max(anchor - COARSE_STEP, MIN_TOTAL_M)

            while true
                record = evaluate_M(Int64(candidate_M))

                if record.qualifies
                    successful_M = Int64(candidate_M)

                    if successful_M == MIN_TOTAL_M
                        minimum_record = record
                        break
                    end

                    candidate_M = max(Int64(successful_M) - COARSE_STEP, MIN_TOTAL_M)
                else
                    failed_M = Int64(candidate_M)
                    break
                end
            end
        end

    # Anchor fails -> search upward, without using M_ref as a ceiling

    else
        failed_M = anchor

        while true
            candidate_M = Int64(failed_M) + COARSE_STEP
            record = evaluate_M(candidate_M)

            if record.qualifies
                successful_M = candidate_M
                break
            end

            failed_M = candidate_M
        end
    end

    # Fine search within failed/successful bracket

    if isnothing(minimum_record)
        isnothing(failed_M) && error("Search did not produce a lower failed bound for alpha=$alpha")
        isnothing(successful_M) && error("Search did not produce an upper successful bound for alpha=$alpha")

        previous_record = evaluate_M(Int64(failed_M))
        previous_record.qualifies && error("Expected lower bracket candidate M=$(failed_M) to fail for alpha=$alpha")

        for total_M in (Int64(failed_M) + FINE_STEP):FINE_STEP:Int64(successful_M)
            record = evaluate_M(Int64(total_M))

            if record.qualifies
                minimum_record = record
                lower_failed_record = previous_record
                break
            end

            previous_record = record
        end
    end

    isnothing(minimum_record) && error("Fine search failed despite an existing successful candidate for alpha=$alpha")

    minimum_record.replicates_completed == NUMBER_SIMULATIONS || error("Accepted M* was not completed to all $NUMBER_SIMULATIONS replicates")
    minimum_record.qualifies || error("Accepted M* does not satisfy the persistent-establishment and rounded-mean-TTE criteria")

    # Scan the entire local verification window, not just the immediate predecessor.
    minimum_record, lower_failed_record = verify_downward_minimum(alpha, minimum_record)

    return qualified_required_effort_result(alpha, anchor, starting_anchor_source, minimum_record, lower_failed_record)
end

# Initialize workers after checkpoint validation

initialize_efffreq_worker!(P, burnin)
for pid in workers()
    remotecall_fetch(initialize_efffreq_worker!, pid, P, burnin)
end

write_configuration_csv()
save_checkpoint(write_output_csvs=true)

# Run alpha search (each subsequent alpha starts from the previous accepted M*)
println("EffFreq: $(STRATEGY_SPECS[STRATEGY].label) | R=$NUMBER_SIMULATIONS | output=$OUTPUT_DIR")

# Process the seasonal-allocation scenarios sequentially.
# The minimum effort obtained for one alpha_S serves as the starting anchor for the next, reducing unnecessary candidate evaluations.
# Completed scenarios are recovered from the checkpoint and their downward verification is repeated using cached results where possible.
function run_alpha_search!()
    previous_M_star = nothing

    for (alpha_index, alpha_value) in enumerate(ALPHAS)
        alpha = Float64(alpha_value)

        saved = completed_result(alpha)
        if !isnothing(saved)
            # Checkpoint-compatible re-verification: reuse all previously evaluated
            # candidate replicates, simulate only missing points in the downward window.
            provisional_record = evaluate_candidate(alpha, Int64(saved.M_star))
            minimum_record, lower_failed_record = verify_downward_minimum(alpha, provisional_record)
            result = qualified_required_effort_result(alpha, Int64(saved.starting_anchor_M), Symbol(saved.starting_anchor_source), minimum_record, lower_failed_record)
            filter!(r -> !(r.strategy == STRATEGY && alpha_code(r.alpha) == alpha_code(alpha)), required_effort_results)
            push!(required_effort_results, result)
            save_checkpoint(write_output_csvs=true)
        else
            if alpha_index == 1 || isnothing(previous_M_star)
                starting_anchor_M = REFERENCE_TOTAL_M
                starting_anchor_source = :reference_M
            else
                starting_anchor_M = Int64(previous_M_star)
                starting_anchor_source = :previous_alpha_M_star
            end

            result = find_minimum_M(alpha, Int64(starting_anchor_M), starting_anchor_source)
            push!(required_effort_results, result)
            save_checkpoint(write_output_csvs=true)
        end

        previous_M_star = Int64(result.M_star)

        println("Accepted alpha=$alpha | M*=$(result.M_star) | P_est=$(result.persistent_establishment_probability_at_M_star) | rounded TTE=$(result.rounded_mean_time_to_establishment_at_M_star)")
    end

    return nothing
end

run_alpha_search!()

write_csvs()
save_checkpoint(write_output_csvs=false)

println("EffFreq complete: $(STRATEGY_SPECS[STRATEGY].label) | Mstar=$REQUIRED_EFFORT_CSV")
