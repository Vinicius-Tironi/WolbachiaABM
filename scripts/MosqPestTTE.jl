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

@everywhere const EFFFREQ_B_WORKER_CONTEXT = Ref{Any}(nothing)

@everywhere function initialize_efffreq_b_worker!(P, burnin)
    EFFFREQ_B_WORKER_CONTEXT[] = (P=P, burnin=burnin)
    return nothing
end

# Establishment-yield grid (fig2)
# R=500 replicates per candidate; alpha_S = 1.0, 1.75, 2.0; daily, weekly, fortnightly.
# Allocate M_HS = alpha_S/2 * M, with integer release sizes; evaluate every 1,000 mosquitoes.
# Establishment requires 30-day moving prevalence >= 0.60 from last release +30 through +730 days.
# Establishment is classified as unsuccessful when P_est < 0.90.
# For P_est >= 0.90, establishment is classified as fast when the rounded
# mean TTE <= 365 days (alpha_S < 2.0) or <= 182 days (alpha_S = 2.0), and as slow when the corresponding TTE threshold is exceeded.
# The grid starts at GRID_START_M and stops after five consecutive qualifying candidates.
# M_ref = 630,000 is a reference, not a ceiling.

# Release design

const VALID_STRATEGIES = [:daily, :weekly, :fortnightly]
const ALPHAS = [1.0, 1.75, 2.0]

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

length(ARGS) == 1 || error("Usage: julia -p <workers> MosqPestTTE.jl <daily|weekly|fortnightly>")

const STRATEGY = Symbol(lowercase(strip(ARGS[1])))
STRATEGY in VALID_STRATEGIES || error("Invalid strategy '$STRATEGY'. Choose daily, weekly, or fortnightly.")

# Replicate / search settings

const NUMBER_SIMULATIONS = 500
const FIRST_SIMULATIONNUMBER = 1

const REQUIRED_SUCCESS_PROBABILITY = 0.90
const REQUIRED_SUCCESSES = ceil(Int64, REQUIRED_SUCCESS_PROBABILITY * NUMBER_SIMULATIONS - 1e-12)

const GRID_STEP = 1_000

const PARALLEL_BATCH_SIZE = max(1, min(NUMBER_SIMULATIONS, max(nworkers(), 1)))

# Establishment settings

const ESTABLISHMENT_THRESHOLD = 0.60
const ROLLING_WINDOW = 30
const WASHOUT_DAYS = 30
const POST_RELEASE_FOLLOWUP_DAYS = 730
const SIMULATION_END_DAY = 1095

# Model settings

const ZONE = 3
const K_VALUE = 25000.0
const BURNIN_DIR = joinpath(MODEL_ROOT, "data", "BurnIn")
const ENVIRONMENTAL_DATA = joinpath(MODEL_ROOT, "data", "ambiental.csv")

# Output / checkpoint

const RUN_ID = "MosqPestTTE_R500"

const OUTPUT_ROOT = get(ENV, "MOSQPESTTTE_OUTPUT_ROOT", joinpath(MODEL_ROOT, "outputs", "MosqPestTTE"))

const OUTPUT_DIR = joinpath(OUTPUT_ROOT, string(STRATEGY))
const REPLICATE_CSV = joinpath(OUTPUT_DIR, "replicates.csv")
const CANDIDATE_CSV = joinpath(OUTPUT_DIR, "candidates.csv")
const CONFIGURATION_CSV = joinpath(OUTPUT_DIR, "config.csv")
const CHECKPOINT_JLS = joinpath(OUTPUT_DIR, "checkpoint.jls")

mkpath(OUTPUT_DIR)

# Basic helpers

percentage_string(x::Real) = isnan(Float64(x)) ? "NA" : "$(round(100.0 * Float64(x), digits=2))%"

csv_value(x) = ismissing(x) || isnothing(x) ? "NA" : string(x)

alpha_tag(alpha::Real) = replace(@sprintf("%.2f", Float64(alpha)), "." => "p")
alpha_code(alpha::Real) = round(Int64, 100 * Float64(alpha))

tte_deadline_days(alpha::Real) = isapprox(Float64(alpha), 2.0; atol=1e-12, rtol=0.0) ? Int64(182) : Int64(365)

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
    1.0 <= alpha <= 2.0 || error("alpha must be between 1 and 2")
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

# For established replicates, determine the first day from which the prevalence threshold remains satisfied until verification ends.
@everywhere function execute_efffreq_b_replicate(simulationnumber::Int64, release_schedule, threshold::Float64, rolling_window::Int64, washout_days::Int64, post_release_followup_days::Int64, simulation_end_day::Int64)
    context = EFFFREQ_B_WORKER_CONTEXT[]
    isnothing(context) && error("Efficiency-frequency worker context is not initialized")

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
    grid_probability_metric=:persistent_establishment,
    qualification_metric=:persistent_probability_and_rounded_mean_tte,
    tte_deadline=:nominal_intervention_duration_days,
    grid_step=GRID_STEP,
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

function checkpoint_input_fingerprints_compatible(saved, current)
    for field in propertynames(current)
        # The search script fingerprint is intentionally allowed to differ so a
        # pure orchestration/runtime fix can resume an existing checkpoint. All
        # model, parameter, data, burn-in, Project, and Manifest fingerprints
        # remain protected. Once the checkpoint is saved again, it receives the
        # current script fingerprint.
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

# Csv fields

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
    :persistent_establishment_yield_per_mosquito,
    :persistent_criterion_met,
    :mean_tte_criterion_met,
    :qualifies,
    :required_success_probability,
    :required_successes,
    :grid_probability_metric,
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

# Checkpoint state

if isfile(CHECKPOINT_JLS)
    checkpoint = deserialize(CHECKPOINT_JLS)

    checkpoint_signature_compatible(
        checkpoint.signature,
        CHECKPOINT_SIGNATURE,
    ) || error(
        "Existing checkpoint uses different fixed settings or input files. " *
        "Use a new output directory or remove the old checkpoint."
    )

    replicate_history = checkpoint.replicate_history
    candidate_history = checkpoint.candidate_history
else
    replicate_history = NamedTuple[]
    candidate_history = NamedTuple[]
end

const replicate_cache = Dict{Tuple{Symbol,Int64,Int64,Int64},NamedTuple}()

for r in replicate_history
    replicate_cache[(r.strategy, alpha_code(r.alpha), r.total_M, r.simulationnumber)] = r
end

# Csv / checkpoint writing

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
    return nothing
end

function save_checkpoint(; write_output_csvs::Bool=false)
    temporary_path = CHECKPOINT_JLS * ".tmp"

    serialize(
        temporary_path,
        (
            signature=CHECKPOINT_SIGNATURE,
            replicate_history=replicate_history,
            candidate_history=candidate_history,
        ),
    )

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
        alpha_step="selected scenarios: 1.0, 1.75, 2.0",
        number_alpha_values=length(ALPHAS),
        reference_M=REFERENCE_TOTAL_M,
        minimum_M=MIN_TOTAL_M,
        yield_grid_step=GRID_STEP,
        number_simulations=NUMBER_SIMULATIONS,
        first_simulationnumber=FIRST_SIMULATIONNUMBER,
        required_success_probability=REQUIRED_SUCCESS_PROBABILITY,
        required_successes=REQUIRED_SUCCESSES,
        grid_probability_metric="persistent_establishment",
        qualification_metric="persistent_probability_and_rounded_mean_tte",
        tte_deadline="365_if_alpha_lt_2__182_if_alpha_eq_2",
        grid_has_no_upper_M_ceiling=true,
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
    matches = [
        r
        for r in candidate_history
        if r.strategy == STRATEGY &&
           alpha_code(r.alpha) == alpha_code(alpha) &&
           r.total_M == total_M &&
           r.replicates_completed == NUMBER_SIMULATIONS
    ]

    return isempty(matches) ? nothing : last(matches)
end

function remove_existing_candidate_record!(alpha::Float64, total_M::Int64)
    filter!(
        r -> !(
            r.strategy == STRATEGY &&
            alpha_code(r.alpha) == alpha_code(alpha) &&
            r.total_M == total_M
        ),
        candidate_history,
    )
    return nothing
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

    records = NamedTuple[
        replicate_cache[(STRATEGY, alpha_code(alpha), total_M, Int64(seed))]
        for seed in seeds
        if haskey(
            replicate_cache,
            (STRATEGY, alpha_code(alpha), total_M, Int64(seed)),
        )
    ]
    sort!(records, by=r -> r.simulationnumber)

    !isempty(records) && println("Resuming alpha=$alpha, M=$total_M: $(length(records))/$NUMBER_SIMULATIONS replicates")

    pending_seeds = Int64[
        Int64(seed)
        for seed in seeds
        if !haskey(
            replicate_cache,
            (STRATEGY, alpha_code(alpha), total_M, Int64(seed)),
        )
    ]

    while !isempty(pending_seeds)
        batch_length = min(PARALLEL_BATCH_SIZE, length(pending_seeds))
        batch_seeds = pending_seeds[1:batch_length]
        deleteat!(pending_seeds, 1:batch_length)

        release_schedule = design.release_schedule

        run_seed = simulationnumber -> execute_efffreq_b_replicate(simulationnumber, release_schedule, ESTABLISHMENT_THRESHOLD, ROLLING_WINDOW, WASHOUT_DAYS, POST_RELEASE_FOLLOWUP_DAYS, SIMULATION_END_DAY)

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

    length(records) == NUMBER_SIMULATIONS || error(
        "Candidate alpha=$alpha M=$total_M did not complete all $NUMBER_SIMULATIONS replicates"
    )

    persistent_successes = count(r -> r.established, records)
    persistent_failures = NUMBER_SIMULATIONS - persistent_successes
    persistent_establishment_probability = persistent_successes / NUMBER_SIMULATIONS

    tte = successful_tte_stats(records)
    tte_deadline = tte_deadline_days(alpha)

    persistent_criterion_met = persistent_successes >= REQUIRED_SUCCESSES
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
        persistent_establishment_yield_per_mosquito=persistent_establishment_probability / total_M,
        persistent_criterion_met=persistent_criterion_met,
        mean_tte_criterion_met=mean_tte_criterion_met,
        qualifies=qualifies,
        required_success_probability=REQUIRED_SUCCESS_PROBABILITY,
        required_successes=REQUIRED_SUCCESSES,
        grid_probability_metric="persistent_establishment",
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

# Evaluate establishment probability and timing across a regular release-effort grid, using a predefined starting value for each strategy and seasonal allocation.
# The grid advances in 1,000-mosquito increments until five consecutive candidates achieve fast establishment.
const GRID_START_M = Dict(
    (:daily,       100) => 540_000,
    (:weekly,      100) => 550_000,
    (:fortnightly, 100) => 590_000,

    (:daily,       175) => 400_000,
    (:weekly,      175) => 410_000,
    (:fortnightly, 175) => 450_000,

    (:daily,       200) => 400_000,
    (:weekly,      200) => 420_000,
    (:fortnightly, 200) => 470_000,
)

function run_yield_grid!(alpha::Float64)
    key = (STRATEGY, alpha_code(alpha))
    haskey(GRID_START_M, key) || error("No grid start configured for $key")

    total_M = Int64(GRID_START_M[key])
    qualifying_streak = 0

    println("Grid alpha=$alpha | starting M=$total_M")

    while true
        record = evaluate_candidate(alpha, total_M)

        qualifying_streak = record.qualifies ? qualifying_streak + 1 : 0
        # Break when five consecutive candidates achieve fast establishment.
        if qualifying_streak >= 5
            println("Grid complete: alpha=$alpha | last M=$total_M")
            return nothing
        end

        total_M += GRID_STEP
    end
end

# Initialize workers after checkpoint validation

initialize_efffreq_b_worker!(P, burnin)
for pid in workers()
    remotecall_fetch(initialize_efffreq_b_worker!, pid, P, burnin)
end

write_configuration_csv()
save_checkpoint(write_output_csvs=true)

println("EffFreq_B: $(STRATEGY_SPECS[STRATEGY].label) | R=$NUMBER_SIMULATIONS | output=$OUTPUT_DIR")

for alpha in ALPHAS
    run_yield_grid!(Float64(alpha))
end

write_csvs()
save_checkpoint(write_output_csvs=false)

println("EffFreq_B complete: $(STRATEGY_SPECS[STRATEGY].label) | candidates=$CANDIDATE_CSV")
