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

mutable struct ImmatureCohort
    number::Int64
    developmental_age::Float64
end


function mosquito_weibull_beta(Tmean::Float64, mosqs, P::MosquitoParameters)
    female = eltype(mosqs) === MosqFem

    a = female ? P.weibull_a_fem : P.weibull_a_male
    b = female ? P.weibull_b_fem : P.weibull_b_male
    c = female ? P.weibull_c_fem : P.weibull_c_male

    beta_lab = exp(a + b * Tmean + c * Tmean^2)
    return beta_lab / P.field_lifespan_multiplier
end

function mosquito_mortality_probability(thermal_age::Float64, beta_today::Float64, k::Float64)
    next_thermal_age = thermal_age + beta_today
    return clamp(1.0 - exp(-(next_thermal_age^k - thermal_age^k)), 0.0, 1.0)
end

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

function setup_initial_free_mosquitoes(number::Int64, P::MosquitoParameters)
    number_females = rand(Binomial(number, P.natural_female_probability))
    return setup_mosq_fem(number_females, false, P), setup_mosq_male(number - number_females, false)
end


function immature_development_beta(Tmean::Float64, P::MosquitoParameters)
    Tmean >= P.immature_Tmax && return 0.0
    median_rate = P.immature_logan_a * exp(P.immature_logan_b * Tmean) * (P.immature_Tmax - Tmean)
    return log(2.0)^(1.0 / P.immature_weibull_k) * median_rate
end

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

function add_immature_cohort!(immatures::Vector{ImmatureCohort}, number::Int64)
    number == 0 && return nothing
    push!(immatures, ImmatureCohort(number, 0.0))
    return nothing
end

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


function setup_released_wolb_random_age!(mosqs::Vector{MosqFem}, P::MosquitoParameters)
    for mosquito in mosqs
        mosquito.age = rand(1:2)
        mosquito.thermal_age = 0.0
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

function setup_released_wolb_random_age!(mosqs::Vector{MosqMale})
    for mosquito in mosqs
        mosquito.age = rand(1:2)
        mosquito.thermal_age = 0.0
    end
    return nothing
end

function release_wolbachia_mosquitoes!(mosq_fem::Vector{MosqFem}, mosq_male::Vector{MosqMale}, release_size::Int64, P::MosquitoParameters)
    released_females = rand(Binomial(release_size, P.release_female_probability))
    released_males = release_size - released_females
    new_females = setup_mosq_fem(released_females, true, P)
    new_males = setup_mosq_male(released_males, true)
    setup_released_wolb_random_age!(new_females, P)
    setup_released_wolb_random_age!(new_males)
    append!(mosq_fem, new_females)
    append!(mosq_male, new_males)
    return released_females, released_males
end