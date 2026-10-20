# Individual adult states. Chronological age is counted in days, while thermal age accumulates temperature-dependent increments used in survival calculations.
mutable struct MosqFem
    wolb::Bool
    age::Int64
    thermal_age::Float64
    mating_age::Int64
    mated::Bool
    mate_wolb::Bool
    gonotrophic_cycle::Int64
    gc_progress::Float64
    gc_threshold::Float64
    gc_complete::Bool
    gc_start_day::Int64
end

mutable struct MosqMale
    wolb::Bool
    age::Int64
    thermal_age::Float64
end

# Immature mosquitoes are grouped into cohorts with a shared developmental age.
mutable struct ImmatureCohort
    number::Int64
    developmental_age::Float64
end

# Calculate the daily physiological-age increment from temperature and sex-specific survival parameters, rescaled from laboratory to field conditions.
function mosquito_weibull_beta(Tmean::Float64, mosqs, P::MosquitoParameters)
    female = eltype(mosqs) === MosqFem

    a = female ? P.weibull_a_fem : P.weibull_a_male
    b = female ? P.weibull_b_fem : P.weibull_b_male
    c = female ? P.weibull_c_fem : P.weibull_c_male

    beta_lab = exp(a + b * Tmean + c * Tmean^2)
    return beta_lab / P.field_lifespan_multiplier
end

# Conditional probability of death over the next daily Weibull-age increment.
function mosquito_mortality_probability(thermal_age::Float64, beta_today::Float64, k::Float64)
    next_thermal_age = thermal_age + beta_today
    return clamp(1.0 - exp(-(next_thermal_age^k - thermal_age^k)), 0.0, 1.0)
end

# Newly emerged adults start at chronological age 1 and thermal age 0.
# Females also receive an individual mating-eligibility age.
function setup_mosq_fem(number::Int64, wolb::Bool, P::MosquitoParameters)
    mosqs = Vector{MosqFem}(undef, number)
    for i in eachindex(mosqs)
        mosqs[i] = MosqFem(wolb, 1, 0.0, rand(P.mating_age_min:P.mating_age_max), false, false, 0, 0.0, 0.0, false, 0)
    end
    return mosqs
end

function setup_mosq_male(number::Int64, wolb::Bool)
    mosqs = Vector{MosqMale}(undef, number)
    for i in eachindex(mosqs)
        mosqs[i] = MosqMale(wolb, 1, 0.0)
    end
    return mosqs
end

# Split the initial Wolbachia-free adult population into females and males.
function setup_initial_free_mosquitoes(number::Int64, P::MosquitoParameters)
    number_females = rand(Binomial(number, P.natural_female_probability))
    return setup_mosq_fem(number_females, false, P), setup_mosq_male(number - number_females, false)
end

# Temperature-dependent developmental increment; development stops outside the configured lower and upper temperature limits.
function immature_development_beta(Tmean::Float64, P::MosquitoParameters)
    (Tmean <= P.immature_Tmin || Tmean >= P.immature_Tmax) && return 0.0
    median_rate = P.immature_logan_a * exp(P.immature_logan_b * Tmean) * (P.immature_Tmax - Tmean)
    return log(2.0)^(1.0 / P.immature_weibull_k) * median_rate
end

# Probability that a cohort member completes development during the daily increment.
function immature_emergence_probability(developmental_age::Float64, beta_today::Float64, k::Float64)
    next_developmental_age = developmental_age + beta_today
    return clamp(1.0 - exp(-(next_developmental_age^k - developmental_age^k)), 0.0, 1.0)
end

function total_immature_population(immatures::Vector{ImmatureCohort})
    total = 0
    for cohort in immatures
        total += cohort.number
    end
    return total
end

# New recruits enter a fresh cohort and begin developing on subsequent days.
function add_immature_cohort!(immatures::Vector{ImmatureCohort}, number::Int64)
    number == 0 && return nothing
    push!(immatures, ImmatureCohort(number, 0.0))
    return nothing
end

# Advance existing cohorts, sample adult emergence, and compact the cohort vector in place by discarding cohorts with no individuals remaining.
function advance_immatures!(immatures::Vector{ImmatureCohort}, Tmean::Float64, P::MosquitoParameters)
    beta_today = immature_development_beta(Tmean, P)
    k = P.immature_weibull_k
    total_emerging = 0
    write_index = 1

    for read_index in eachindex(immatures)
        @inbounds cohort = immatures[read_index]
        emerging = rand(Binomial(cohort.number, immature_emergence_probability(cohort.developmental_age, beta_today, k)))
        cohort.number -= emerging
        cohort.developmental_age += beta_today
        total_emerging += emerging

        if cohort.number > 0
            @inbounds immatures[write_index] = cohort
            write_index += 1
        end
    end

    resize!(immatures, write_index - 1)
    return total_emerging
end

function split_emerging_mosquitoes(number::Int64, female_probability::Float64)
    number_females = rand(Binomial(number, female_probability))
    return number_females, number - number_females
end

# Apply daily adult mortality before aging the survivors.
# Infected mosquitoes use the additional Wolbachia lifespan multiplier for their thermal-age increment.
# Here, wolb_lifespan_multiplier = 1.0 implies no lifespan fitness is applied.
# We consider fitness costs in fecundity and fertility for the wMel strain -- for other strains, implementing lifespan fitness cost is straightforward. 
function increase_mosquito_age!(mosqs, Tmean::Float64, P::MosquitoParameters)
    beta_today = mosquito_weibull_beta(Tmean, mosqs, P)
    beta_wolb_today = beta_today / P.wolb_lifespan_multiplier
    k = eltype(mosqs) === MosqFem ? P.weibull_k_fem : P.weibull_k_male
    write_index = 1

    for read_index in eachindex(mosqs)
        @inbounds mosquito = mosqs[read_index]
        beta = mosquito.wolb ? beta_wolb_today : beta_today
        death_probability = mosquito_mortality_probability(mosquito.thermal_age, beta, k)

        if rand() >= death_probability
            mosquito.thermal_age += beta
            mosquito.age += 1
            @inbounds mosqs[write_index] = mosquito
            write_index += 1
        end
    end

    resize!(mosqs, write_index - 1)
    return nothing
end

# Released adults are assigned age 1 or 2. Age-2 individuals receive one physiological-age increment based on the previous day's temperature.
# Released females start unmated and without an active gonotrophic cycle.
function setup_released_wolb_random_age!(mosqs::Vector{MosqFem}, P::MosquitoParameters, previous_Tmean::Float64)
    previous_beta = mosquito_weibull_beta(previous_Tmean, mosqs, P) / P.wolb_lifespan_multiplier
    for mosquito in mosqs
        mosquito.age = rand(1:2)
        mosquito.thermal_age = mosquito.age == 1 ? 0.0 : previous_beta
        mosquito.mating_age = rand(P.mating_age_min:P.mating_age_max)
        mosquito.mated = false
        mosquito.mate_wolb = false
        mosquito.gonotrophic_cycle = 0
        mosquito.gc_progress = 0.0
        mosquito.gc_threshold = 0.0
        mosquito.gc_complete = false
        mosquito.gc_start_day = 0
    end
    return nothing
end

function setup_released_wolb_random_age!(mosqs::Vector{MosqMale}, P::MosquitoParameters, previous_Tmean::Float64)
    previous_beta = mosquito_weibull_beta(previous_Tmean, mosqs, P) / P.wolb_lifespan_multiplier
    for mosquito in mosqs
        mosquito.age = rand(1:2)
        mosquito.thermal_age = mosquito.age == 1 ? 0.0 : previous_beta
    end
    return nothing
end

# Sample the sex composition of each release, initialize infected adults, and append them to the existing population without changing the total release size.
function release_wolbachia_mosquitoes!(mosq_fem::Vector{MosqFem}, mosq_male::Vector{MosqMale}, release_size::Int64, P::MosquitoParameters, previous_Tmean::Float64)
    released_females = rand(Binomial(release_size, P.release_female_probability))
    released_males = release_size - released_females
    new_females = setup_mosq_fem(released_females, true, P)
    new_males = setup_mosq_male(released_males, true)
    setup_released_wolb_random_age!(new_females, P, previous_Tmean)
    setup_released_wolb_random_age!(new_males, P, previous_Tmean)
    append!(mosq_fem, new_females)
    append!(mosq_male, new_males)
    return released_females, released_males
end
