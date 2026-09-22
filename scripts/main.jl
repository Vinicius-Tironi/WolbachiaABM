using CSV
using Parameters
using Distributions
using Random
using Serialization
using Statistics
include("parameters.jl")
include("mosquito.jl")
include("interaction.jl")
const FIRST_VERIFICATION_DAY = 366
const LAST_VERIFICATION_DAY = 1095

function load_tmean_period(path::String, zone::Int64, number_days::Int64; from_end::Bool = false)
    number_days == 0 && return Float64[]
    zone_name = "Zone $zone"
    rows = [(date = string(row.date), Tmean = Float64(row.Tmean)) for row in CSV.File(path) if String(row.zone) == zone_name]
    sort!(rows, by = row -> row.date)
    selected_rows = from_end ? rows[(end - number_days + 1):end] : rows[1:number_days]
    return [row.Tmean for row in selected_rows]
end

function load_tmean(P::MosquitoParameters)
    burnin_tmean = load_tmean_period(P.burnin_environmental_data_path, P.zone, P.burn_in_days; from_end = true)
    calendar_tmean = load_tmean_period(P.environmental_data_path, P.zone, P.calendar_days)
    return vcat(burnin_tmean, calendar_tmean)
end

function population_counts(mosq_fem::Vector{MosqFem}, mosq_male::Vector{MosqMale})
    female_wolb = 0
    not_yet_eligible = 0
    eligible_unmated = 0
    mated = 0
    for female in mosq_fem
        female.wolb && (female_wolb += 1)
        if female.mated
            mated += 1
        elseif female.age < female.mating_age
            not_yet_eligible += 1
        else
            eligible_unmated += 1
        end
    end
    male_wolb = count(male -> male.wolb, mosq_male)
    return (
        length(mosq_fem) - female_wolb,
        female_wolb,
        length(mosq_male) - male_wolb,
        male_wolb,
        not_yet_eligible,
        eligible_unmated,
        mated
    )
end

function simulate_mosquito_population(
    P::MosquitoParameters,
    Tmean::Vector{Float64},
    mosq_fem::Vector{MosqFem},
    mosq_male::Vector{MosqMale},
    immature_free::Vector{ImmatureCohort},
    immature_wolb::Vector{ImmatureCohort};
    release_schedule = Dict{Int64,Int64}(),
    first_simulation_day::Int64 = 1,
    number_days::Int64 = length(Tmean),
    calendar_origin_day::Int64 = P.burn_in_days,
    burnin_snapshot_day::Int64 = P.burn_in_days
)
    last_simulation_day = first_simulation_day + number_days - 1
    length(Tmean) >= last_simulation_day || error(
        "Temperature vector does not cover the requested " *
        "simulation period"
    )
    cumulative_Tmean = cumsum(Tmean)
    female_free_ctr = zeros(Int64, number_days + 1)
    female_wolb_ctr = zeros(Int64, number_days + 1)
    male_free_ctr = zeros(Int64, number_days + 1)
    male_wolb_ctr = zeros(Int64, number_days + 1)
    not_yet_eligible_female_ctr = zeros(Int64, number_days + 1)
    eligible_unmated_female_ctr = zeros(Int64, number_days + 1)
    mated_female_ctr = zeros(Int64, number_days + 1)
    bloodmeal_ctr = zeros(Int64, number_days)
    released_female_ctr = zeros(Int64, number_days)
    released_male_ctr = zeros(Int64, number_days)
    #ci_dead_offspring_ctr = zeros(Int64, number_days)
    unhatched_eggs_ctr = zeros(Int64, number_days)
    density_dead_offspring_ctr = zeros(Int64, number_days)
    female_free_ctr[1],
    female_wolb_ctr[1],
    male_free_ctr[1],
    male_wolb_ctr[1],
    not_yet_eligible_female_ctr[1],
    eligible_unmated_female_ctr[1],
    mated_female_ctr[1] = population_counts(mosq_fem, mosq_male)
    burnin_snapshot = nothing
    for step in 1:number_days
        simulation_day = first_simulation_day + step - 1
        calendar_day = simulation_day - calendar_origin_day
        temperature_today = Tmean[simulation_day]
        release_size = calendar_day > 0 ? get(release_schedule, calendar_day, 0) : 0
        if release_size > 0
            released_female_ctr[step],
            released_male_ctr[step] = release_wolbachia_mosquitoes!(
                mosq_fem,
                mosq_male,
                release_size,
                P
            )
        end
        bloodmeal_ctr[step] = bloodmeal_interaction!(
            mosq_fem,
            simulation_day,
            P
        )
        gonotrophic_cycle_interaction!(mosq_fem, temperature_today, P)
        total_immatures = total_immature_population(immature_free) + total_immature_population(immature_wolb)
        new_free,
        new_wolb,
        unhatched_eggs_ctr[step],
        density_dead_offspring_ctr[step] = reproductive_interaction!(
            mosq_fem,
            mosq_male,
            total_immatures,
            simulation_day,
            cumulative_Tmean,
            P
        )
        increase_mosquito_age!(mosq_fem, temperature_today, P)
        increase_mosquito_age!(mosq_male, temperature_today, P)
        emerging_free = advance_immatures!(immature_free, temperature_today, P)
        emerging_wolb = advance_immatures!(immature_wolb, temperature_today, P)
        projected_population =
            length(mosq_fem) +
            length(mosq_male) +
            emerging_free +
            emerging_wolb +
            total_immature_population(immature_free) +
            total_immature_population(immature_wolb) +
            new_free +
            new_wolb
        projected_population <= 5_000_000 || error(
            "Mosquito population exceeds 5,000,000 " *
            "on simulation day $simulation_day"
        )
        emerging_free_females,
        emerging_free_males = split_emerging_mosquitoes(
            emerging_free,
            P.natural_female_probability
        )
        emerging_wolb_females,
        emerging_wolb_males = split_emerging_mosquitoes(
            emerging_wolb,
            P.natural_female_probability
        )
        emerging_free_females > 0 && append!(mosq_fem, setup_mosq_fem(emerging_free_females, false, P))
        emerging_free_males > 0 && append!(mosq_male, setup_mosq_male(emerging_free_males, false))
        emerging_wolb_females > 0 && append!(mosq_fem, setup_mosq_fem(emerging_wolb_females, true, P))
        emerging_wolb_males > 0 && append!(mosq_male, setup_mosq_male(emerging_wolb_males, true))
        add_immature_cohort!(immature_free, new_free)
        add_immature_cohort!(immature_wolb, new_wolb)
        if burnin_snapshot_day > 0 && simulation_day == burnin_snapshot_day
            burnin_snapshot = (
                mosq_fem = deepcopy(mosq_fem),
                mosq_male = deepcopy(mosq_male),
                immature_free = deepcopy(immature_free),
                immature_wolb = deepcopy(immature_wolb)
            )
        end
        female_free_ctr[step + 1],
        female_wolb_ctr[step + 1],
        male_free_ctr[step + 1],
        male_wolb_ctr[step + 1],
        not_yet_eligible_female_ctr[step + 1],
        eligible_unmated_female_ctr[step + 1],
        mated_female_ctr[step + 1] = population_counts(mosq_fem, mosq_male)
    end
    free_ctr = female_free_ctr .+ male_free_ctr
    wolb_ctr = female_wolb_ctr .+ male_wolb_ctr
    female_ctr = female_free_ctr .+ female_wolb_ctr
    male_ctr = male_free_ctr .+ male_wolb_ctr
    state_calendar_days = collect((first_simulation_day - 1 - calendar_origin_day):(last_simulation_day - calendar_origin_day))
    return (
        state_calendar_days = state_calendar_days,
        free_ctr = free_ctr,
        wolb_ctr = wolb_ctr,
        female_ctr = female_ctr,
        male_ctr = male_ctr,
        female_free_ctr = female_free_ctr,
        female_wolb_ctr = female_wolb_ctr,
        male_free_ctr = male_free_ctr,
        male_wolb_ctr = male_wolb_ctr,
        not_yet_eligible_female_ctr = not_yet_eligible_female_ctr,
        eligible_unmated_female_ctr = eligible_unmated_female_ctr,
        mated_female_ctr = mated_female_ctr,
        bloodmeal_ctr = bloodmeal_ctr,
        released_female_ctr = released_female_ctr,
        released_male_ctr = released_male_ctr,
    #    ci_dead_offspring_ctr = ci_dead_offspring_ctr,
        unhatched_eggs_ctr = unhatched_eggs_ctr,
        density_dead_offspring_ctr = density_dead_offspring_ctr,
        burnin_state = burnin_snapshot,
        mosq_fem = mosq_fem,
        mosq_male = mosq_male,
        immature_free = immature_free,
        immature_wolb = immature_wolb
    )
end

function main(simulationnumber::Int64, P::MosquitoParameters; release_schedule = Dict{Int64,Int64}())
    Random.seed!(simulationnumber)
    P.K > 0.0 || error("K must be greater than zero")
    P.initial_mosq_free >= 0 || error("initial_mosq_free must be non-negative")
    Tmean = load_tmean(P)
    mosq_fem, mosq_male = setup_initial_free_mosquitoes(P.initial_mosq_free, P)
    return simulate_mosquito_population(
        P,
        Tmean,
        mosq_fem,
        mosq_male,
        ImmatureCohort[],
        ImmatureCohort[];
        release_schedule = release_schedule,
        first_simulation_day = 1,
        number_days = P.burn_in_days + P.calendar_days,
        calendar_origin_day = P.burn_in_days,
        burnin_snapshot_day = P.burn_in_days
    )
end

function load_burnin(burnin_dir::String, zone::Int64, K::Float64)
    path = joinpath(burnin_dir, "zone$(zone)", "K$(Int(round(K))).jls")
    burnin = deserialize(path)
    burnin.zone == zone || error("Burn-in belongs to Zone $(burnin.zone), not Zone $zone")
    isapprox(Float64(burnin.K), K; atol = 0.0, rtol = 1e-12) || error("Burn-in K = $(burnin.K), requested K = $K")
    burnin.absolute_day == burnin.burn_in_days || error("Burn-in absolute clock is inconsistent")
    length(burnin.burnin_tmean) == burnin.burn_in_days || error("Burn-in temperature history is inconsistent")
    return burnin, path
end

function posterior_parameters(burnin, calendar_days::Int64, weekly_release::Bool, environmental_data_path::String)
    kwargs = Dict{Symbol,Any}(name => getfield(burnin.parameters, name) for name in fieldnames(MosquitoParameters))
    current_parameters = MosquitoParameters(initial_mosq_free = 0, zone = burnin.zone)
    kwargs[:initial_mosq_free] = 0
    kwargs[:burn_in_days] = 0
    kwargs[:calendar_days] = calendar_days
    kwargs[:weekly_release] = weekly_release
    kwargs[:delta] = current_parameters.delta
 #   kwargs[:eta] = current_parameters.eta
    kwargs[:wolb_lifespan_multiplier] = current_parameters.wolb_lifespan_multiplier
    kwargs[:release_female_probability] = current_parameters.release_female_probability
    kwargs[:K] = Float64(burnin.K)
    kwargs[:environmental_data_path] = environmental_data_path
    return MosquitoParameters(; kwargs...)
end

function main_from_burnin(
    simulationnumber::Int64,
    P::MosquitoParameters,
    burnin;
    release_schedule = Dict{Int64,Int64}()
)
    Random.seed!(simulationnumber)
    P.K > 0.0 || error("K must be greater than zero")
    P.zone == burnin.zone || error("Simulation zone does not match burn-in zone")
    isapprox(P.K, Float64(burnin.K); atol = 0.0, rtol = 1e-12) || error("Simulation K does not match burn-in K")
    calendar_tmean = load_tmean_period(P.environmental_data_path, P.zone, P.calendar_days)
    Tmean = vcat(burnin.burnin_tmean, calendar_tmean)
    return simulate_mosquito_population(
        P,
        Tmean,
        deepcopy(burnin.mosq_fem),
        deepcopy(burnin.mosq_male),
        deepcopy(burnin.immature_free),
        deepcopy(burnin.immature_wolb);
        release_schedule = release_schedule,
        first_simulation_day = burnin.absolute_day + 1,
        number_days = P.calendar_days,
        calendar_origin_day = burnin.absolute_day,
        burnin_snapshot_day = 0
    )
end

function rolling_mean(values::Vector{Float64}, window::Int64)
    output = fill(NaN, length(values))
    for i in window:length(values)
        output[i] = mean(@view values[(i - window + 1):i])
    end
    return output
end

function period_indices(days::Vector{Int64}, first_day::Int64, last_day::Int64)
    return findall(i -> first_day <= days[i] <= last_day, eachindex(days))
end

function period_minimum(days::Vector{Int64}, values::Vector{Float64}, first_day::Int64, last_day::Int64)
    indices = [i for i in period_indices(days, first_day, last_day) if !isnan(values[i])]
    return isempty(indices) ? NaN : Base.minimum(values[indices])
end

function period_stable(days::Vector{Int64}, values::Vector{Float64}, first_day::Int64, last_day::Int64, threshold::Float64)
    indices = [i for i in period_indices(days, first_day, last_day) if !isnan(values[i])]
    return !isempty(indices) && all(values[i] >= threshold for i in indices)
end

function evaluate_2019_establishment(
    result,
    release_schedule;
    threshold::Float64 = 0.60,
    rolling_window::Int64 = 30,
    washout_days::Int64 = 30
)
    positive_release_days = sort(Int64[day for (day, size) in release_schedule if size > 0])
    isempty(positive_release_days) && error("Release schedule contains no positive releases")
    first_release_day = first(positive_release_days)
    last_release_day = last(positive_release_days)
    calendar_indices = findall(day -> day >= 1, result.state_calendar_days)
    calendar_days = Int64.(result.state_calendar_days[calendar_indices])
    adult_free = result.free_ctr[calendar_indices]
    adult_wolb = result.wolb_ctr[calendar_indices]
    adult_total = adult_free .+ adult_wolb
    adult_wolb_prevalence = [adult_total[i] > 0 ? adult_wolb[i] / adult_total[i] : 0.0 for i in eachindex(adult_total)]
    rolling_prevalence = rolling_mean(adult_wolb_prevalence, rolling_window)
    washout_end_day = last_release_day + washout_days
    verification_start_day = max(FIRST_VERIFICATION_DAY, washout_end_day)
    verification_end_day = min(LAST_VERIFICATION_DAY, last(calendar_days))
    verification_minimum = period_minimum(calendar_days, rolling_prevalence, verification_start_day, verification_end_day)
    full_2019_minimum = period_minimum(calendar_days, rolling_prevalence, FIRST_VERIFICATION_DAY, verification_end_day)
    established_2019 = period_stable(
        calendar_days,
        rolling_prevalence,
        verification_start_day,
        verification_end_day,
        threshold
    )
    return (
        calendar_days = calendar_days,
        adult_wolb_prevalence = adult_wolb_prevalence,
        rolling_prevalence = rolling_prevalence,
        first_release_day = first_release_day,
        last_release_day = last_release_day,
        washout_end_day = washout_end_day,
        verification_start_day = verification_start_day,
        verification_end_day = verification_end_day,
        verification_minimum = verification_minimum,
        full_2019_minimum = full_2019_minimum,
        established_2019 = established_2019
    )
end

function mean_without_nan(values::Vector{Float64})
    valid = [value for value in values if !isnan(value)]
    return isempty(valid) ? NaN : mean(valid)
end

function evaluate_release_probability(
    P::MosquitoParameters,
    release_schedule;
    use_burnin::Bool,
    burnin = nothing,
    number_simulations::Int64 = 1,
    first_simulationnumber::Int64 = 1,
    establishment_threshold::Float64 = 0.60,
    rolling_window::Int64 = 30,
    washout_days::Int64 = 30
)
    establishment_successes = 0
    verification_minima = Float64[]
    full_2019_minima = Float64[]
    representative_result = nothing
    representative_metrics = nothing
    for replicate in 1:number_simulations
        simulationnumber = first_simulationnumber + replicate - 1
        result = use_burnin ?
            main_from_burnin(simulationnumber, P, burnin; release_schedule = release_schedule) :
            main(simulationnumber, P; release_schedule = release_schedule)
        metrics = evaluate_2019_establishment(
            result,
            release_schedule;
            threshold = establishment_threshold,
            rolling_window = rolling_window,
            washout_days = washout_days
        )
        establishment_successes += metrics.established_2019 ? 1 : 0
        push!(verification_minima, metrics.verification_minimum)
        push!(full_2019_minima, metrics.full_2019_minimum)
        if replicate == 1
            representative_result = result
            representative_metrics = metrics
        end
    end
    return (
        establishment_probability = establishment_successes / number_simulations,
        mean_verification_minimum = mean_without_nan(verification_minima),
        mean_full_2019_minimum = mean_without_nan(full_2019_minima),
        representative_result = representative_result,
        representative_metrics = representative_metrics
    )
end

function search_minimum_annual_release(
    P::MosquitoParameters,
    candidate_annual_totals,
    schedule_builder::Function;
    use_burnin::Bool,
    burnin = nothing,
    number_simulations::Int64 = 1,
    first_simulationnumber::Int64 = 1,
    required_success_probability::Float64 = 0.90,
    establishment_threshold::Float64 = 0.60,
    rolling_window::Int64 = 30,
    washout_days::Int64 = 30
)
    records = NamedTuple[]
    minimum = nothing
    for annual_total in candidate_annual_totals
        release_schedule = schedule_builder(Int64(annual_total))
        positive_release_days = sort(Int64[day for (day, size) in release_schedule if size > 0])
        evaluation = evaluate_release_probability(
            P,
            release_schedule;
            use_burnin = use_burnin,
            burnin = burnin,
            number_simulations = number_simulations,
            first_simulationnumber = first_simulationnumber,
            establishment_threshold = establishment_threshold,
            rolling_window = rolling_window,
            washout_days = washout_days
        )
        record = (
            annual_total = Int64(annual_total),
            total_released = sum(values(release_schedule)),
            number_release_events = length(positive_release_days),
            first_release_day = first(positive_release_days),
            last_release_day = last(positive_release_days),
            establishment_probability = evaluation.establishment_probability,
            mean_verification_minimum = evaluation.mean_verification_minimum,
            mean_full_2019_minimum = evaluation.mean_full_2019_minimum
        )
        push!(records, record)
        if evaluation.establishment_probability >= required_success_probability
            minimum = (
                annual_total = Int64(annual_total),
                release_schedule = release_schedule,
                record = record,
                representative_result = evaluation.representative_result,
                representative_metrics = evaluation.representative_metrics
            )
            break
        end
    end
    return (records = records, minimum = minimum)
end
