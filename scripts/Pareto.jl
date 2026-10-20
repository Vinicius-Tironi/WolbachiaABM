using Distributed
using Serialization
using SHA
using Statistics

const MODEL_ROOT = dirname(@__DIR__)
const PROJECT_ROOT = MODEL_ROOT
const MODEL_ENTRYPOINT = joinpath(@__DIR__, "main.jl")

@everywhere include($MODEL_ENTRYPOINT)

@everywhere const PARETO_WORKER_CONTEXT = Ref{Any}(nothing)

@everywhere function initialize_pareto_worker!(P, burnin, use_burnin::Bool)
    PARETO_WORKER_CONTEXT[] = (P = P, burnin = burnin, use_burnin = use_burnin)
    return nothing
end

# Pareto analysis (fig3): release frequency, required effort, and establishment time.
# R=500 replicates; daily, weekly, and fortnightly strategies; HS-only releases (alpha=2.0).
# For each number of release events k, search M*(k) such that P_est >= 0.90.
# Candidate classification uses establishment probability only, independently of TTE.
# Early stopping: 450 successes for acceptance or 51 failures for rejection.
# Accepted M*(k) is evaluated with all 500 replicates.
# Two-sided search: 5,000-mosquito coarse and 1,000-mosquito fine steps, followed by 20,000-mosquito downward verification.
# Each k starts from the previously accepted M*(k+1), except k_max, which uses a fixed anchor.
# Establishment: moving prevalence >= 0.60, verified from last release + washout through +730 days.
# Fast and slow establishment are distinguished by mean TTE among established replicates.
# Fast establishment: rounded mean TTE <= 182 days; slow establishment: > 182 days.


# Release design

const RELEASE_YEAR_DAYS = 365

const FAVORABLE_FIRST_WEEK = 18
const FAVORABLE_LAST_WEEK = 43
const FAVORABLE_FIRST_DAY = 1 + 7 * (FAVORABLE_FIRST_WEEK - 1)
const FAVORABLE_LAST_DAY = 7 * FAVORABLE_LAST_WEEK

const STRATEGY_SPECS = Dict(
    :daily => (label = "Daily", interval_days = 1, events_per_year = 365),
    :weekly => (label = "Weekly", interval_days = 7, events_per_year = 52),
    :fortnightly => (label = "Fortnightly", interval_days = 14, events_per_year = 26)
)

length(ARGS) == 1 || error("Usage: julia -p <workers> Pareto.jl <daily|weekly|fortnightly>")

const SELECTED_STRATEGY = Symbol(lowercase(ARGS[1]))
haskey(STRATEGY_SPECS, SELECTED_STRATEGY) || error("Unknown strategy '$(ARGS[1])'. Choose daily, weekly, or fortnightly.")
const STRATEGIES = Symbol[SELECTED_STRATEGY]

const ALLOCATION = :hs
const ALPHA = 2.0

const INITIAL_FULL_CAMPAIGN_ANCHORS = Dict(
    :daily => 396_000,
    :weekly => 404_000,
    :fortnightly => 451_000
)

# Replicate and search settings

const NUMBER_SIMULATIONS = 500
const FIRST_SIMULATIONNUMBER = 1

const REQUIRED_SUCCESS_PROBABILITY = 0.90
const REQUIRED_SUCCESSES = ceil(Int64, REQUIRED_SUCCESS_PROBABILITY * NUMBER_SIMULATIONS - 1e-12)
const REJECTION_FAILURES = NUMBER_SIMULATIONS - REQUIRED_SUCCESSES + 1

const COARSE_STEP = 5_000
const FINE_STEP = 1_000
const DOWNWARD_VERIFY_WINDOW = 20_000

const MIN_TOTAL_N = 1_000

const ROLLING_WINDOW = parse(Int, get(ENV, "PARETO_DMW", "30"))
const WASHOUT_DAYS = parse(Int, get(ENV, "PARETO_WASHOUT", "30"))
const MAX_TOTAL_N = parse(Int, get(ENV, "PARETO_MAX_N", "700000"))
const RUN_SERIES = get(ENV, "PARETO_RUN_SERIES", "washoutcmp_v1")

ROLLING_WINDOW >= 1 || error("PARETO_DMW must be at least 1")
WASHOUT_DAYS >= 0 || error("PARETO_WASHOUT must be non-negative")
MAX_TOTAL_N >= MIN_TOTAL_N || error("PARETO_MAX_N must be at least $MIN_TOTAL_N")

const RUN_ID = "Pareto_R500"

const PARALLEL_BATCH_SIZE = max(1, min(NUMBER_SIMULATIONS, max(nworkers(), 1), parse(Int, get(ENV, "PARETO_BATCH_SIZE", string(max(nworkers(), 1))))))

# Establishment settings

const POST_RELEASE_FOLLOWUP_DAYS = 730
const ESTABLISHMENT_THRESHOLD = 0.60

# Simulation settings

use_burnin = true
zone = 3
K = 25000.0
burn_in_days = 364
calendar_days = RELEASE_YEAR_DAYS + POST_RELEASE_FOLLOWUP_DAYS
initial_mosq_free = 10_000

burnin_dir = joinpath(MODEL_ROOT, "data", "BurnIn")
burnin_environmental_data_path = joinpath(MODEL_ROOT, "data", "ambiental16.csv")
environmental_data_path = joinpath(MODEL_ROOT, "data", "ambiental.csv")

# Output and checkpoint

output_dir = joinpath(MODEL_ROOT, "outputs", "Pareto", string(SELECTED_STRATEGY))
mkpath(output_dir)

checkpoint_path = joinpath(output_dir, "checkpoint.jls")
replicate_csv_path = joinpath(output_dir, "replicates.csv")
candidate_csv_path = joinpath(output_dir, "candidates.csv")
required_effort_csv_path = joinpath(output_dir, "Mstar.csv")
configuration_csv_path = joinpath(output_dir, "config.csv")

# Basic helpers

percentage_string(x::Float64) = isnan(x) ? "NA" : "$(round(100.0 * x, digits = 2))%"
csv_value(x) = ismissing(x) || isnothing(x) ? "NA" : string(x)

function mean_or_nan(values)
    valid = Float64[Float64(x) for x in values if !ismissing(x) && !isnothing(x) && !isnan(Float64(x))]
    return isempty(valid) ? NaN : mean(valid)
end

function successful_tte_stats(records)
    values = Float64[Float64(r.time_to_establishment) for r in records if r.established && !ismissing(r.time_to_establishment) && !isnan(Float64(r.time_to_establishment))]

    isempty(values) && return (
        n = 0,
        mean = NaN,
        rounded_mean = missing,
        sd = NaN,
        se = NaN,
        lower95 = NaN,
        upper95 = NaN,
        q025 = NaN,
        median = NaN,
        q975 = NaN,
    )

    n = length(values)
    mu = mean(values)
    sd_value = n > 1 ? std(values) : 0.0
    se_value = sd_value / sqrt(n)
    margin = 1.96 * se_value

    return (
        n = n,
        mean = mu,
        rounded_mean = round(Int64, mu),
        sd = sd_value,
        se = se_value,
        lower95 = mu - margin,
        upper95 = mu + margin,
        q025 = quantile(values, 0.025),
        median = median(values),
        q975 = quantile(values, 0.975),
    )
end

function split_integer_total(total::Int64, number_parts::Int64)
    total >= 0 || error("total must be non-negative")
    number_parts >= 1 || error("number_parts must be at least 1")

    base_size, remainder = divrem(total, number_parts)
    sizes = fill(Int64(base_size), number_parts)

    for i in 1:remainder
        sizes[i] += 1
    end

    return sizes
end

round_up_to_step(value::Int64, step::Int64) = value <= 0 ? step : Int64(cld(value, step) * step)

function fmt_int(x)
    s = string(Int(round(x)))
    parts = String[]

    while lastindex(s) > 3
        pushfirst!(parts, s[end-2:end])
        s = s[1:end-3]
    end

    pushfirst!(parts, s)
    return join(parts, ",")
end

# Release schedule

function annual_release_days(strategy::Symbol)
    spec = STRATEGY_SPECS[strategy]
    return Int64[1 + spec.interval_days * i for i in 0:(spec.events_per_year - 1)]
end

function full_campaign_release_days(strategy::Symbol)
    return Int64[day for day in annual_release_days(strategy) if FAVORABLE_FIRST_DAY <= day <= FAVORABLE_LAST_DAY]
end

# Construct a release schedule using the first k events of the full HS campaign.
function campaign_release_schedule(strategy::Symbol, number_events::Int64, total_N::Int64)
    total_N > 0 || error("total_N must be greater than zero")

    all_release_days = full_campaign_release_days(strategy)
    1 <= number_events <= length(all_release_days) || error("number_events must be between 1 and $(length(all_release_days))")

    release_days = all_release_days[1:number_events]
    release_sizes = split_integer_total(total_N, number_events)
    release_schedule = Dict(release_days[i] => release_sizes[i] for i in eachindex(release_days))

    return release_schedule
end

function release_design_metrics(strategy::Symbol, k::Int64, total_N::Int64)
    release_schedule = campaign_release_schedule(strategy, k, total_N)
    release_days = sort(collect(keys(release_schedule)))

    first_release_day = first(release_days)
    last_release_day = last(release_days)
    campaign_duration_days = last_release_day - first_release_day
    number_release_events = k
    release_interval_days = STRATEGY_SPECS[strategy].interval_days
    mosquitoes_per_release = total_N / number_release_events
    average_release_intensity = total_N / (campaign_duration_days + 1)

    return (
        release_schedule = release_schedule,
        first_release_day = first_release_day,
        last_release_day = last_release_day,
        campaign_duration_days = campaign_duration_days,
        number_release_events = number_release_events,
        release_interval_days = release_interval_days,
        mosquitoes_per_release = mosquitoes_per_release,
        average_release_intensity = average_release_intensity
    )
end

# Model parameters and burn-in

if use_burnin
    burnin, burnin_path = load_burnin(burnin_dir, zone, K)
    P = posterior_parameters(burnin, calendar_days, false, environmental_data_path)
else
    burnin = nothing

    P = MosquitoParameters(
        initial_mosq_free = initial_mosq_free,
        zone = zone,
        burn_in_days = burn_in_days,
        calendar_days = calendar_days,
        weekly_release = false,
        K = K,
        burnin_environmental_data_path = burnin_environmental_data_path,
        environmental_data_path = environmental_data_path
    )
end

# Establishment and timing metrics

@everywhere function persistent_establishment_day_until(calendar_days, rolling_prevalence, threshold::Float64, first_release_day::Int64, end_day::Int64)
    start_index = findfirst(day -> day >= first_release_day, calendar_days)
    end_index = findlast(day -> day <= end_day, calendar_days)

    isnothing(start_index) && return nothing
    isnothing(end_index) && return nothing
    start_index > end_index && return nothing

    for i in start_index:end_index
        if !isnan(rolling_prevalence[i]) &&
           rolling_prevalence[i] >= threshold &&
           all(!isnan(rolling_prevalence[j]) && rolling_prevalence[j] >= threshold for j in i:end_index)
            return Int64(calendar_days[i])
        end
    end

    return nothing
end

# Evaluate persistent Wolbachia establishment after the release campaign.
# Establishment requires the 30-day moving prevalence to remain above the threshold throughout the post-release verification period.
# For established replicates, TTE is measured from the first release.
@everywhere function evaluate_post_release_establishment(result, release_schedule; threshold::Float64 = ESTABLISHMENT_THRESHOLD, rolling_window::Int64 = ROLLING_WINDOW, washout_days::Int64 = WASHOUT_DAYS, post_release_followup_days::Int64 = POST_RELEASE_FOLLOWUP_DAYS)
    positive_release_days = sort(Int64[day for (day, size) in release_schedule if size > 0])

    first_release_day = first(positive_release_days)
    last_release_day = last(positive_release_days)

    calendar_indices = findall(day -> day >= 1, result.state_calendar_days)
    calendar_days_local = Int64.(result.state_calendar_days[calendar_indices])
    adult_free = result.free_ctr[calendar_indices]
    adult_wolb = result.wolb_ctr[calendar_indices]
    adult_total = adult_free .+ adult_wolb

    adult_wolb_prevalence = [adult_total[i] > 0 ? adult_wolb[i] / adult_total[i] : 0.0 for i in eachindex(adult_total)]
    rolling_prevalence = rolling_mean(adult_wolb_prevalence, rolling_window)

    full_post_release_start_day = last_release_day + 1
    full_post_release_end_day = last_release_day + post_release_followup_days
    verification_start_day = last_release_day + washout_days
    verification_end_day = full_post_release_end_day

    last(calendar_days_local) >= verification_end_day || error(
        "Simulation ends on day $(last(calendar_days_local)), but verification requires day $verification_end_day"
    )

    verification_minimum = period_minimum(calendar_days_local, rolling_prevalence, verification_start_day, verification_end_day)

    full_post_release_minimum = period_minimum(calendar_days_local, rolling_prevalence, full_post_release_start_day, full_post_release_end_day)

    established = period_stable(calendar_days_local, rolling_prevalence, verification_start_day, verification_end_day, threshold)
    # Determine the earliest day from which the prevalence threshold remains satisfied through the verification endpoint.
    #This is calculated only for replicates that meet the persistent-establishment criterion.
    establishment_day = if established
        persistent_establishment_day_until(calendar_days_local, rolling_prevalence, threshold, first_release_day, verification_end_day)
    else
        nothing
    end

    time_to_establishment = isnothing(establishment_day) ? missing : establishment_day - first_release_day
    establishment_lag_after_campaign = isnothing(establishment_day) ? missing : establishment_day - last_release_day

    return (
        first_release_day = first_release_day,
        last_release_day = last_release_day,
        campaign_duration_days = last_release_day - first_release_day,
        verification_start_day = verification_start_day,
        verification_end_day = verification_end_day,
        verification_minimum = verification_minimum,
        full_post_release_minimum = full_post_release_minimum,
        established = established,
        establishment_day = establishment_day,
        time_to_establishment = time_to_establishment,
        establishment_lag_after_campaign = establishment_lag_after_campaign
    )
end

@everywhere function execute_pareto_replicate(simulationnumber::Int64, release_schedule, threshold::Float64, rolling_window::Int64, washout_days::Int64, post_release_followup_days::Int64)
    context = PARETO_WORKER_CONTEXT[]
    isnothing(context) && error("Pareto worker context was not initialized")

    result = context.use_burnin ?
        main_from_burnin(
            simulationnumber,
            context.P,
            context.burnin;
            release_schedule = release_schedule
        ) :
        main(
            simulationnumber,
            context.P;
            release_schedule = release_schedule
        )

    return evaluate_post_release_establishment(result, release_schedule; threshold = threshold, rolling_window = rolling_window, washout_days = washout_days, post_release_followup_days = post_release_followup_days)
end

# Checkpoint signature

function file_sha256(path::String)
    return open(path, "r") do io
        bytes2hex(sha256(io))
    end
end

optional_file_sha256(path::String) = isfile(path) ? file_sha256(path) : "not_present"

input_fingerprints = (
    main = file_sha256(joinpath(@__DIR__, "main.jl")),
    parameters = file_sha256(joinpath(@__DIR__, "parameters.jl")),
    mosquito = file_sha256(joinpath(@__DIR__, "mosquito.jl")),
    interaction = file_sha256(joinpath(@__DIR__, "interaction.jl")),
    project = optional_file_sha256(joinpath(dirname(PROJECT_ROOT), "Project.toml")),
    manifest = optional_file_sha256(joinpath(dirname(PROJECT_ROOT), "Manifest.toml")),
    environmental_data = file_sha256(environmental_data_path),
    burnin_environmental_data = file_sha256(burnin_environmental_data_path),
    burnin = use_burnin ? file_sha256(burnin_path) : "not_used"
)

checkpoint_signature = (
    run_id = RUN_ID,
    julia_version = string(VERSION),
    zone = zone,
    K = K,
    calendar_days = calendar_days,
    strategies = Tuple(STRATEGIES),
    allocation = ALLOCATION,
    alpha = ALPHA,
    number_simulations = NUMBER_SIMULATIONS,
    first_simulationnumber = FIRST_SIMULATIONNUMBER,
    required_success_probability = REQUIRED_SUCCESS_PROBABILITY,
    establishment_threshold = ESTABLISHMENT_THRESHOLD,
    rolling_window = ROLLING_WINDOW,
    washout_days = WASHOUT_DAYS,
    post_release_followup_days = POST_RELEASE_FOLLOWUP_DAYS,
    favorable_first_week = FAVORABLE_FIRST_WEEK,
    favorable_last_week = FAVORABLE_LAST_WEEK,
    coarse_step = COARSE_STEP,
    fine_step = FINE_STEP,
    downward_verify_window = DOWNWARD_VERIFY_WINDOW,
    min_total_N = MIN_TOTAL_N,
    input_fingerprints = input_fingerprints
)

# Before resuming a previous run, check that the saved settings and model inputs are compatible with the current analysis.
# This prevents simulation results from different configurations from being combined.
function checkpoint_signature_compatible(saved_signature, current_signature)
    return all(hasproperty(saved_signature, field) && getproperty(saved_signature, field) == getproperty(current_signature, field) for field in propertynames(current_signature))
end

if isfile(checkpoint_path)
    checkpoint = deserialize(checkpoint_path)

    checkpoint_signature_compatible(checkpoint.signature, checkpoint_signature) || error(
        "Existing checkpoint uses different fixed settings. Change RUN_ID/output directory or remove the old checkpoint."
    )

    saved_max_total_N = Int64(checkpoint.max_total_N)
    MAX_TOTAL_N >= saved_max_total_N || error(
        "MAX_TOTAL_N cannot be decreased within a run " *
        "(checkpoint = $saved_max_total_N, current = $MAX_TOTAL_N)."
    )

    replicate_history = checkpoint.replicate_history
    candidate_history = checkpoint.candidate_history
    duration_release_results = checkpoint.duration_release_results
else
    replicate_history = NamedTuple[]
    candidate_history = NamedTuple[]
    duration_release_results = NamedTuple[]
end

replicate_cache = Dict{Tuple{Symbol,Int64,Int64,Int64},NamedTuple}()

for r in replicate_history
    replicate_cache[(r.strategy, r.k, r.total_N, r.simulationnumber)] = r
end

initialize_pareto_worker!(P, burnin, use_burnin)
for worker_id in workers()
    remotecall_fetch(initialize_pareto_worker!, worker_id, P, burnin, use_burnin)
end

# CSV fields

const REPLICATE_FIELDS = Symbol[
    :strategy,
    :allocation,
    :alpha,
    :k,
    :max_events,
    :total_N,
    :number_release_events,
    :release_interval_days,
    :mosquitoes_per_release,
    :first_release_day,
    :last_release_day,
    :campaign_duration_days,
    :average_release_intensity,
    :simulationnumber,
    :established,
    :verification_start_day,
    :verification_end_day,
    :verification_minimum,
    :full_post_release_minimum,
    :establishment_day,
    :time_to_establishment,
    :establishment_lag_after_campaign
]

const CANDIDATE_FIELDS = Symbol[
    :strategy,
    :allocation,
    :alpha,
    :k,
    :max_events,
    :total_N,
    :number_release_events,
    :release_interval_days,
    :mosquitoes_per_release,
    :first_release_day,
    :last_release_day,
    :campaign_duration_days,
    :average_release_intensity,
    :replicate_budget,
    :classification,
    :completed_replicates,
    :successes,
    :failures,
    :required_successes,
    :rejection_failures,
    :observed_success_fraction,
    :establishment_probability,
    :verification_start_day,
    :verification_end_day,
    :mean_verification_minimum,
    :mean_full_post_release_minimum,
    :mean_establishment_day,
    :mean_time_to_establishment,
    :rounded_mean_time_to_establishment,
    :successful_time_replicates,
    :std_time_to_establishment,
    :se_time_to_establishment,
    :lower_95ci_time_to_establishment,
    :upper_95ci_time_to_establishment,
    :q025_time_to_establishment,
    :median_time_to_establishment,
    :q975_time_to_establishment,
    :mean_establishment_lag_after_campaign
]

const REQUIRED_EFFORT_FIELDS = Symbol[
    :strategy,
    :allocation,
    :alpha,
    :k,
    :max_events,
    :starting_anchor_N,
    :starting_anchor_source,
    :total_N,
    :number_release_events,
    :release_interval_days,
    :mosquitoes_per_release,
    :first_release_day,
    :last_release_day,
    :campaign_duration_days,
    :average_release_intensity,
    :successes,
    :failures,
    :number_simulations,
    :establishment_probability,
    :verification_start_day,
    :verification_end_day,
    :mean_verification_minimum,
    :mean_full_post_release_minimum,
    :mean_establishment_day,
    :mean_time_to_establishment,
    :rounded_mean_time_to_establishment,
    :successful_time_replicates,
    :std_time_to_establishment,
    :se_time_to_establishment,
    :lower_95ci_time_to_establishment,
    :upper_95ci_time_to_establishment,
    :q025_time_to_establishment,
    :median_time_to_establishment,
    :q975_time_to_establishment,
    :mean_establishment_lag_after_campaign,
    :lower_failed_N,
    :lower_failed_completed_replicates,
    :lower_failed_successes,
    :lower_failed_failures,
    :lower_failed_observed_success_fraction
]

# CSV and checkpoint writing

function write_namedtuple_csv(path::String, fields::Vector{Symbol}, records)
    temporary_path = path * ".tmp"

    open(temporary_path, "w") do io
        println(io, join(string.(fields), ","))

        for r in records
            println(io, join([csv_value(getproperty(r, field)) for field in fields], ","))
        end
    end

    mv(temporary_path, path; force = true)
end

function write_csvs()
    write_namedtuple_csv(replicate_csv_path, REPLICATE_FIELDS, replicate_history)
    write_namedtuple_csv(candidate_csv_path, CANDIDATE_FIELDS, candidate_history)
    write_namedtuple_csv(required_effort_csv_path, REQUIRED_EFFORT_FIELDS, duration_release_results)
end

function write_configuration_csv()
    row = (
        run_id = RUN_ID,
        strategy = string(SELECTED_STRATEGY),
        allocation = string(ALLOCATION),
        alpha = ALPHA,
        number_simulations = NUMBER_SIMULATIONS,
        first_simulationnumber = FIRST_SIMULATIONNUMBER,
        required_success_probability = REQUIRED_SUCCESS_PROBABILITY,
        required_successes = REQUIRED_SUCCESSES,
        rejection_failures = REJECTION_FAILURES,
        search_criterion = "persistent_establishment_probability_only",
        coarse_step = COARSE_STEP,
        fine_step = FINE_STEP,
        downward_verify_window = DOWNWARD_VERIFY_WINDOW,
        minimum_N = MIN_TOTAL_N,
        maximum_N = MAX_TOTAL_N,
        establishment_threshold = ESTABLISHMENT_THRESHOLD,
        rolling_window_days = ROLLING_WINDOW,
        washout_days = WASHOUT_DAYS,
        post_release_followup_days = POST_RELEASE_FOLLOWUP_DAYS,
        simulation_end_day = calendar_days,
        favorable_first_week = FAVORABLE_FIRST_WEEK,
        favorable_last_week = FAVORABLE_LAST_WEEK,
        favorable_first_day = FAVORABLE_FIRST_DAY,
        favorable_last_day = FAVORABLE_LAST_DAY,
        zone = zone,
        K = K,
        burnin_path = use_burnin ? burnin_path : "not_used",
        environmental_data = environmental_data_path,
        julia_version = string(VERSION),
        parallel_batch_size = PARALLEL_BATCH_SIZE,
    )

    write_namedtuple_csv(configuration_csv_path, Symbol[propertynames(row)...], [row])
    return nothing
end

function save_checkpoint(; write_output_csvs::Bool = false)
    temporary_path = checkpoint_path * ".tmp"

    serialize(
        temporary_path,
        (
            signature = checkpoint_signature,
            max_total_N = MAX_TOTAL_N,
            replicate_history = replicate_history,
            candidate_history = candidate_history,
            duration_release_results = duration_release_results
        )
    )

    mv(temporary_path, checkpoint_path; force = true)
    write_output_csvs && write_csvs()
end

# Resume lookups

function existing_candidate_record(strategy::Symbol, k::Int64, total_N::Int64, replicate_budget::Int64 = NUMBER_SIMULATIONS)
    matches = [
        r
        for r in candidate_history
        if r.strategy == strategy &&
           r.k == k &&
           r.total_N == total_N &&
           r.replicate_budget == replicate_budget
    ]

    return isempty(matches) ? nothing : last(matches)
end

function remove_existing_candidate_record!(strategy::Symbol, k::Int64, total_N::Int64, replicate_budget::Int64 = NUMBER_SIMULATIONS)
    filter!(
        r -> !(
            r.strategy == strategy &&
            r.k == k &&
            r.total_N == total_N &&
            r.replicate_budget == replicate_budget
        ),
        candidate_history
    )
end

function already_completed(strategy::Symbol, k::Int64, replicate_budget::Int64 = NUMBER_SIMULATIONS)
    return any(
        r.strategy == strategy &&
        r.k == k &&
        r.number_simulations == replicate_budget
        for r in duration_release_results
    )
end

function completed_result(strategy::Symbol, k::Int64, replicate_budget::Int64 = NUMBER_SIMULATIONS)
    matches = [
        r
        for r in duration_release_results
        if r.strategy == strategy &&
           r.k == k &&
           r.number_simulations == replicate_budget
    ]

    return isempty(matches) ? nothing : last(matches)
end

# Replicate records

function build_replicate_record(strategy::Symbol, k::Int64, total_N::Int64, simulationnumber::Int64, design, metrics)
    return (
        strategy = strategy,
        allocation = ALLOCATION,
        alpha = ALPHA,
        k = k,
        max_events = length(full_campaign_release_days(strategy)),
        total_N = total_N,
        number_release_events = design.number_release_events,
        release_interval_days = design.release_interval_days,
        mosquitoes_per_release = design.mosquitoes_per_release,
        first_release_day = design.first_release_day,
        last_release_day = design.last_release_day,
        campaign_duration_days = design.campaign_duration_days,
        average_release_intensity = design.average_release_intensity,
        simulationnumber = simulationnumber,
        established = metrics.established,
        verification_start_day = metrics.verification_start_day,
        verification_end_day = metrics.verification_end_day,
        verification_minimum = metrics.verification_minimum,
        full_post_release_minimum = metrics.full_post_release_minimum,
        establishment_day = isnothing(metrics.establishment_day) ? missing : metrics.establishment_day,
        time_to_establishment = metrics.time_to_establishment,
        establishment_lag_after_campaign = metrics.establishment_lag_after_campaign
    )
end

function store_replicate_record!(record)
    key = (record.strategy, record.k, record.total_N, record.simulationnumber)

    if !haskey(replicate_cache, key)
        replicate_cache[key] = record
        push!(replicate_history, record)
    end

    return record
end

# Evaluate candidates using a maximum of 500 simulation replicates.
# Early-classification thresholds are calculated automatically from the replicate budget and required establishment probability (P_est >= 0.90).
# With R=500, 450 successes guarantee acceptance and 51 failures guarantee rejection, regardless of the remaining replicate outcomes.
# Each selected M*(k) is subsequently completed to all 500 replicates to obtain full establishment and TTE statistics.
function evaluate_candidate(strategy::Symbol, k::Int64, total_N::Int64; force_complete::Bool = false)
    existing = existing_candidate_record(strategy, k, total_N)

    if !isnothing(existing) && (!force_complete || existing.completed_replicates == NUMBER_SIMULATIONS)
        return existing
    end

    design = release_design_metrics(strategy, k, total_N)
    seeds = FIRST_SIMULATIONNUMBER:(FIRST_SIMULATIONNUMBER + NUMBER_SIMULATIONS - 1)

    records = NamedTuple[replicate_cache[(strategy, k, total_N, Int64(seed))] for seed in seeds if haskey(replicate_cache, (strategy, k, total_N, Int64(seed)))]
    sort!(records, by = r -> r.simulationnumber)

    successes = count(r -> r.established, records)
    failures = length(records) - successes
    cached = length(records)

    cached > 0 && println("Resuming k=$k, N=$total_N: $cached/$NUMBER_SIMULATIONS replicates")

    pending_seeds = Int64[Int64(seed) for seed in seeds if !haskey(replicate_cache, (strategy, k, total_N, Int64(seed)))]

    while !isempty(pending_seeds)
        if !force_complete &&
           (successes >= REQUIRED_SUCCESSES || failures >= REJECTION_FAILURES)
            break
        end

        batch_length = min(PARALLEL_BATCH_SIZE, length(pending_seeds))
        batch_seeds = pending_seeds[1:batch_length]
        deleteat!(pending_seeds, 1:batch_length)

        release_schedule = design.release_schedule
        threshold = ESTABLISHMENT_THRESHOLD
        rolling_window = ROLLING_WINDOW
        washout_days = WASHOUT_DAYS
        post_release_followup_days = POST_RELEASE_FOLLOWUP_DAYS

        run_seed = simulationnumber -> execute_pareto_replicate(simulationnumber, release_schedule, threshold, rolling_window, washout_days, post_release_followup_days)

        batch_metrics = isempty(workers()) ?
            map(run_seed, batch_seeds) :
            pmap(run_seed, WorkerPool(workers()), batch_seeds; batch_size = 1)

        for (simulationnumber, metrics) in zip(batch_seeds, batch_metrics)
            record = build_replicate_record(strategy, k, total_N, simulationnumber, design, metrics)

            store_replicate_record!(record)
            push!(records, record)
            successes += record.established ? 1 : 0
            failures += record.established ? 0 : 1

        end

        sort!(records, by = r -> r.simulationnumber)
        save_checkpoint(write_output_csvs = true)
    end

    completed_replicates = length(records)
    classification = successes >= REQUIRED_SUCCESSES ? :pass : :fail
    observed_success_fraction = successes / completed_replicates
    establishment_probability = completed_replicates == NUMBER_SIMULATIONS ? successes / NUMBER_SIMULATIONS : NaN
    tte = successful_tte_stats(records)

    remove_existing_candidate_record!(strategy, k, total_N)

    candidate_record = (
        strategy = strategy,
        allocation = ALLOCATION,
        alpha = ALPHA,
        k = k,
        max_events = length(full_campaign_release_days(strategy)),
        total_N = total_N,
        number_release_events = design.number_release_events,
        release_interval_days = design.release_interval_days,
        mosquitoes_per_release = design.mosquitoes_per_release,
        first_release_day = design.first_release_day,
        last_release_day = design.last_release_day,
        campaign_duration_days = design.campaign_duration_days,
        average_release_intensity = design.average_release_intensity,
        replicate_budget = NUMBER_SIMULATIONS,
        classification = classification,
        completed_replicates = completed_replicates,
        successes = successes,
        failures = failures,
        required_successes = REQUIRED_SUCCESSES,
        rejection_failures = REJECTION_FAILURES,
        observed_success_fraction = observed_success_fraction,
        establishment_probability = establishment_probability,
        verification_start_day = design.last_release_day + WASHOUT_DAYS,
        verification_end_day = design.last_release_day + POST_RELEASE_FOLLOWUP_DAYS,
        mean_verification_minimum = mean_or_nan([r.verification_minimum for r in records]),
        mean_full_post_release_minimum = mean_or_nan([r.full_post_release_minimum for r in records]),
        mean_establishment_day = mean_or_nan([r.establishment_day for r in records]),
        mean_time_to_establishment = tte.mean,
        rounded_mean_time_to_establishment = tte.rounded_mean,
        successful_time_replicates = tte.n,
        std_time_to_establishment = tte.sd,
        se_time_to_establishment = tte.se,
        lower_95ci_time_to_establishment = tte.lower95,
        upper_95ci_time_to_establishment = tte.upper95,
        q025_time_to_establishment = tte.q025,
        median_time_to_establishment = tte.median,
        q975_time_to_establishment = tte.q975,
        mean_establishment_lag_after_campaign = mean_or_nan([r.establishment_lag_after_campaign for r in records])
    )

    push!(candidate_history, candidate_record)
    save_checkpoint(write_output_csvs = true)

    println("Candidate k=$k | N=$total_N | classification=$classification | successes=$successes/$completed_replicates")

    return candidate_record
end

# Verify every 1,000-mosquito candidate within 20,000 below the provisional M*, retaining the smallest qualifying effort.
# Previously evaluated candidates are recovered from the checkpoint.
# The immediate predecessor of the selected M* must also fail.
function verify_downward_minimum_N(strategy::Symbol, k::Int64, provisional_record)
    provisional_record.classification == :pass || error("Provisional M* does not pass for k=$k")
    selected = provisional_record
    lower_bound = max(MIN_TOTAL_N, Int64(provisional_record.total_N) - DOWNWARD_VERIFY_WINDOW)

    for candidate_N in lower_bound:FINE_STEP:(Int64(provisional_record.total_N) - FINE_STEP)
        record = evaluate_candidate(strategy, k, Int64(candidate_N))
        if record.classification == :pass && record.total_N < selected.total_N
            selected = record
        end
    end

    lower_failed_record = nothing
    if selected.total_N > MIN_TOTAL_N
        lower_failed_record = evaluate_candidate(strategy, k, Int64(selected.total_N) - FINE_STEP)
        lower_failed_record.classification == :fail || error(
            "Immediate predecessor of verified M* passes for k=$k M*=$(selected.total_N)"
        )
    end

    println("Downward verification k=$k | provisional=$(provisional_record.total_N) | selected=$(selected.total_N) | window=$DOWNWARD_VERIFY_WINDOW")
    return selected, lower_failed_record
end

# Two-sided M* search: force complete accepted M*

function find_replicated_minimum_N(strategy::Symbol, k::Int64, starting_anchor_N::Int64, starting_anchor_source::Symbol)
    starting_anchor_N > 0 || error("starting_anchor_N must be greater than zero")
    starting_anchor_N <= MAX_TOTAL_N || error(
        "Starting anchor N = $starting_anchor_N exceeds MAX_TOTAL_N = $MAX_TOTAL_N"
    )

    anchor = max(round_up_to_step(starting_anchor_N, FINE_STEP), MIN_TOTAL_N)
    anchor <= MAX_TOTAL_N || error("Rounded starting anchor exceeds MAX_TOTAL_N")

    evaluate_N(N::Int64) = evaluate_candidate(strategy, k, N)

    println("Searching k=$k | anchor N=$anchor ($(starting_anchor_source))")

    anchor_record = evaluate_N(anchor)

    successful_N = nothing
    failed_N = nothing
    minimum_record = nothing
    lower_failed_record = nothing

    if anchor_record.classification == :pass
        successful_N = anchor

        if anchor == MIN_TOTAL_N
            minimum_record = anchor_record
        else
            candidate_N = max(anchor - COARSE_STEP, MIN_TOTAL_N)

            while true
                record = evaluate_N(Int64(candidate_N))

                if record.classification == :pass
                    successful_N = Int64(candidate_N)

                    if successful_N == MIN_TOTAL_N
                        minimum_record = record
                        break
                    end

                    candidate_N = max(Int64(successful_N) - COARSE_STEP, MIN_TOTAL_N)
                else
                    failed_N = Int64(candidate_N)
                    break
                end
            end
        end

    else
        failed_N = anchor

        while Int64(failed_N) < MAX_TOTAL_N
            candidate_N = min(Int64(failed_N) + COARSE_STEP, MAX_TOTAL_N)
            record = evaluate_N(Int64(candidate_N))

            if record.classification == :pass
                successful_N = Int64(candidate_N)
                break
            end

            failed_N = Int64(candidate_N)
            candidate_N == MAX_TOTAL_N && break
        end

        isnothing(successful_N) && return nothing
    end

    if isnothing(minimum_record)
        isnothing(failed_N) && error("Two-sided search did not produce a lower failed bound")
        isnothing(successful_N) && error("Two-sided search did not produce an upper successful bound")

        previous_record = evaluate_N(Int64(failed_N))
        previous_record.classification == :fail || error("Expected lower bracket candidate to be failed")

        for total_N in (Int64(failed_N) + FINE_STEP):FINE_STEP:Int64(successful_N)
            record = evaluate_N(Int64(total_N))

            if record.classification == :pass
                minimum_record = record
                lower_failed_record = previous_record
                break
            end

            previous_record = record
        end
    end

    isnothing(minimum_record) && error("Fine search failed despite an existing successful candidate")

    # Verify below provisional M* before completing the accepted candidate to R=500.
    minimum_record, lower_failed_record = verify_downward_minimum_N(strategy, k, minimum_record)
    
    # Complete the selected minimum to all 500 replicates, ensuring that its reported establishment probability and TTE are based on the full simulation sample.
    minimum_record = evaluate_candidate(strategy, k, minimum_record.total_N; force_complete = true)

    minimum_record.classification == :pass || error("Accepted candidate did not remain successful after completion")
    minimum_record.completed_replicates == NUMBER_SIMULATIONS || error("Accepted candidate was not completed to all replicates")

    design = release_design_metrics(strategy, k, minimum_record.total_N)

    return (
        strategy = strategy,
        allocation = ALLOCATION,
        alpha = ALPHA,
        k = k,
        max_events = length(full_campaign_release_days(strategy)),
        starting_anchor_N = anchor,
        starting_anchor_source = starting_anchor_source,
        total_N = minimum_record.total_N,
        number_release_events = design.number_release_events,
        release_interval_days = design.release_interval_days,
        mosquitoes_per_release = design.mosquitoes_per_release,
        first_release_day = design.first_release_day,
        last_release_day = design.last_release_day,
        campaign_duration_days = design.campaign_duration_days,
        average_release_intensity = design.average_release_intensity,
        successes = minimum_record.successes,
        failures = minimum_record.failures,
        number_simulations = NUMBER_SIMULATIONS,
        establishment_probability = minimum_record.establishment_probability,
        verification_start_day = minimum_record.verification_start_day,
        verification_end_day = minimum_record.verification_end_day,
        mean_verification_minimum = minimum_record.mean_verification_minimum,
        mean_full_post_release_minimum = minimum_record.mean_full_post_release_minimum,
        mean_establishment_day = minimum_record.mean_establishment_day,
        mean_time_to_establishment = minimum_record.mean_time_to_establishment,
        rounded_mean_time_to_establishment = minimum_record.rounded_mean_time_to_establishment,
        successful_time_replicates = minimum_record.successful_time_replicates,
        std_time_to_establishment = minimum_record.std_time_to_establishment,
        se_time_to_establishment = minimum_record.se_time_to_establishment,
        lower_95ci_time_to_establishment = minimum_record.lower_95ci_time_to_establishment,
        upper_95ci_time_to_establishment = minimum_record.upper_95ci_time_to_establishment,
        q025_time_to_establishment = minimum_record.q025_time_to_establishment,
        median_time_to_establishment = minimum_record.median_time_to_establishment,
        q975_time_to_establishment = minimum_record.q975_time_to_establishment,
        mean_establishment_lag_after_campaign = minimum_record.mean_establishment_lag_after_campaign,
        lower_failed_N = isnothing(lower_failed_record) ? missing : lower_failed_record.total_N,
        lower_failed_completed_replicates = isnothing(lower_failed_record) ? missing : lower_failed_record.completed_replicates,
        lower_failed_successes = isnothing(lower_failed_record) ? missing : lower_failed_record.successes,
        lower_failed_failures = isnothing(lower_failed_record) ? missing : lower_failed_record.failures,
        lower_failed_observed_success_fraction = isnothing(lower_failed_record) ? missing : lower_failed_record.observed_success_fraction,
        downward_verified_window = DOWNWARD_VERIFY_WINDOW
    )
end

write_configuration_csv()
save_checkpoint(write_output_csvs = true)

println("Pareto: $(STRATEGY_SPECS[SELECTED_STRATEGY].label) | R=$NUMBER_SIMULATIONS | output=$output_dir")

# Begin with the full HS campaign and progressively reduce the number of release events. The accepted minimum effort for each campaign
# provides the starting anchor for the next, shorter campaign.
# Completed configurations are recovered directly from the checkpoint.
for strategy in STRATEGIES
    spec = STRATEGY_SPECS[strategy]
    release_days = full_campaign_release_days(strategy)
    maximum_k = length(release_days)
    previous_replicated_N = nothing

    for k in maximum_k:-1:1

        if already_completed(strategy, Int64(k))
            saved = completed_result(strategy, Int64(k))
            previous_replicated_N = saved.total_N

            println("Already complete: k=$k | M*=$(saved.total_N)")

            continue
        end

        if k == maximum_k
            starting_anchor_N = INITIAL_FULL_CAMPAIGN_ANCHORS[strategy]
            starting_anchor_source = :initial_full_campaign
        elseif !isnothing(previous_replicated_N)
            starting_anchor_N = Int64(previous_replicated_N)
            starting_anchor_source = :previous_replicated_k
        else
            error(
                "No previous replicated result available for $(spec.label), k = $k"
            )
        end

        result = find_replicated_minimum_N(strategy, Int64(k), Int64(starting_anchor_N), starting_anchor_source)

        if isnothing(result)
            println(
                "No successful N found up to MAX_TOTAL_N = $MAX_TOTAL_N for $(spec.label) | HS | k = $k"
            )

            save_checkpoint(write_output_csvs = true)
            break
        end

        push!(duration_release_results, result)
        previous_replicated_N = result.total_N
        save_checkpoint(write_output_csvs = true)

        println("Accepted k=$k | M*=$(result.total_N) | P_est=$(result.establishment_probability) | mean TTE=$(result.mean_time_to_establishment)")
    end
end

save_checkpoint(write_output_csvs = true)

println("Pareto complete: $(STRATEGY_SPECS[SELECTED_STRATEGY].label) | Mstar=$required_effort_csv_path")
println("Checkpoint = $checkpoint_path")
