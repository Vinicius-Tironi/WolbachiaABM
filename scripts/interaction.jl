# Temperature-dependent daily progression of an active gonotrophic cycle.
function gonotrophic_cycle_rate(Tmean::Float64, P::MosquitoParameters)
    if P.gc_briere_T0 < Tmean < P.gc_briere_Tm
        return P.gc_briere_c * Tmean * (Tmean - P.gc_briere_T0) * sqrt(P.gc_briere_Tm - Tmean)
    end
    return 0.0
end

# Scale egg production according to temperature, relative to the reference temperature.
function fecundity_temperature_multiplier(Tmean::Float64, P::MosquitoParameters)
    if P.fecundity_briere_T0 < Tmean < P.fecundity_briere_Tm
        numerator = Tmean * (Tmean - P.fecundity_briere_T0) * (P.fecundity_briere_Tm - Tmean)^P.fecundity_briere_b
        Tref = P.fecundity_reference_temperature
        denominator = Tref * (Tref - P.fecundity_briere_T0) * (P.fecundity_briere_Tm - Tref)^P.fecundity_briere_b
        return numerator / denominator
    end
    return 0.0
end

# Average temperature from the initiating blood meal through the current day, using cumulative sums to avoid traversing the temperature history.
function mean_gc_temperature(female::MosqFem, simulation_day::Int64, cumulative_Tmean::Vector{Float64})
    start_day = female.gc_start_day
    temperature_sum = cumulative_Tmean[simulation_day] - (start_day > 1 ? cumulative_Tmean[start_day - 1] : 0.0)
    return temperature_sum / (simulation_day - start_day + 1)
end


# A mated female can take a blood meal when no cycle is active or awaiting egg laying.
# The meal begins a new cycle if the reproductive-cycle limit has not been reached; subsequent meals are still counted but do not start further cycles.
function bloodmeal_interaction!(mosq_fem::Vector{MosqFem}, simulation_day::Int64, P::MosquitoParameters)
    total_bloodmeals_today = 0
    bloodmeal_distribution = Bernoulli(P.bite_probability)

    for female in mosq_fem
        !female.mated && continue
        (female.gc_threshold > 0.0 || female.gc_complete) && continue
        rand(bloodmeal_distribution) == 0 && continue

        total_bloodmeals_today += 1
        if female.gonotrophic_cycle < P.max_gonotrophic_cycles
            female.gonotrophic_cycle += 1
            female.gc_progress = 0.0
            female.gc_threshold = rand(LogNormal(-0.5 * P.gc_sigma^2, P.gc_sigma))
            female.gc_complete = false
            female.gc_start_day = simulation_day
        end
    end

    return total_bloodmeals_today
end

# Advance each active cycle with today's temperature. Completed cycles become eligible for egg laying in the reproductive interaction later that day.
function gonotrophic_cycle_interaction!(mosq_fem::Vector{MosqFem}, Tmean::Float64, P::MosquitoParameters)
    development_today = gonotrophic_cycle_rate(Tmean, P)

    for female in mosq_fem
        female.gc_threshold <= 0.0 && continue
        female.gc_complete && continue
        female.gc_progress += development_today
        female.gc_progress >= female.gc_threshold && (female.gc_complete = true)
    end

    return nothing
end


# Process mating and egg laying for the current adult population, then return newly recruited offspring for addition to the immature cohorts.
function reproductive_interaction!(
    mosq_fem::Vector{MosqFem},
    mosq_male::Vector{MosqMale},
    total_immatures::Int64,
    simulation_day::Int64,
    cumulative_Tmean::Vector{Float64},
    P::MosquitoParameters
)
    new_free = 0
    new_wolb = 0
    unhatched_eggs = 0
    density_dead_offspring = 0

    mating_distribution = Bernoulli(P.daily_mating_probability)
    # Existing immature abundance determines recruitment for all eggs laid today.
    # Maynard Smith–Slatkin/Shepherd functional form
    recruitment_probability = P.s0 / (1.0 + (total_immatures / P.K)^P.q)

    for female in mosq_fem
        # Unmated females can mate once they reach mating age. The selected mate determines their reproductive cross.
        if !female.mated &&
           female.mating_age > 0 &&
           female.age >= female.mating_age &&
           !isempty(mosq_male) &&
           rand(mating_distribution) == 1

            female.mated = true
            female.mate_wolb = rand(mosq_male).wolb
        end

        # Egg laying requires a completed cycle. Later cycles have their own fecundity and fertility multipliers.
        (!female.mated ||
         female.gonotrophic_cycle <= 0 ||
         female.gonotrophic_cycle > P.max_gonotrophic_cycles ||
         !female.gc_complete) && continue

        k = female.gonotrophic_cycle
        alpha_k = P.fecundity_scaling[k]
        beta_k = P.fertility_scaling[k]

        fecundity_temperature = fecundity_temperature_multiplier(mean_gc_temperature(female, simulation_day, cumulative_Tmean), P)

        # Egg production and hatching depend on both parents' Wolbachia status.
        # Cross notation is female × male (w: infected; f: free).
        # Cytoplasmic incompatibility in f × w is captured by hatch_fw.
        if female.wolb
            if female.mate_wolb
                # w × w
                eggs_r = P.eggs_r_ww
                eggs_mu = P.eggs_mu_ww
                hatch_probability = beta_k * P.hatch_ww
            else
                # w × f
                eggs_r = P.eggs_r_wf
                eggs_mu = P.eggs_mu_wf
                hatch_probability = beta_k * P.hatch_wf
            end
        else
            if female.mate_wolb
                # f × w
                eggs_r = P.eggs_r_fw
                eggs_mu = P.eggs_mu_fw
                hatch_probability = beta_k * P.hatch_fw
            else
                # f × f
                eggs_r = P.eggs_r_ff
                eggs_mu = P.eggs_mu_ff
                hatch_probability = beta_k * P.hatch_ff
            end
        end

        # Sample egg production from the cross-specific Negative Binomial distribution, adjusted for cycle number and the mean temperature of the cycle.
        temperature_adjusted_mean = alpha_k * fecundity_temperature * eggs_mu
        n_eggs = temperature_adjusted_mean > 0.0 ?
            rand(NegativeBinomial(eggs_r, eggs_r / (eggs_r + temperature_adjusted_mean))) : 0


        # Apply maternal transmission and cross-specific hatching before recruitment.
        # Count unhatched eggs separately from hatchlings lost to density dependence.
        if female.wolb
            # Infected females can produce both infected and uninfected offspring.
            viable_wolb = rand(Binomial(n_eggs, P.delta))
            viable_free = n_eggs - viable_wolb
            hatched_wolb = rand(Binomial(viable_wolb, hatch_probability))
            hatched_free = rand(Binomial(viable_free, hatch_probability))
            unhatched_eggs += n_eggs - hatched_wolb - hatched_free

            recruited_wolb = rand(Binomial(hatched_wolb, recruitment_probability))
            recruited_free = rand(Binomial(hatched_free, recruitment_probability))
            density_dead_offspring += hatched_wolb + hatched_free - recruited_wolb - recruited_free

            new_wolb += recruited_wolb
            new_free += recruited_free

        else
            # Uninfected females produce only Wolbachia-free offspring;
            # the mating cross still determines their hatching probability.
            hatched_free = rand(Binomial(n_eggs, hatch_probability))
            unhatched_eggs += n_eggs - hatched_free

            recruited_free = rand(Binomial(hatched_free, recruitment_probability))
            density_dead_offspring += hatched_free - recruited_free
            new_free += recruited_free
        end

        # Clear the completed cycle so another blood meal can initiate the next one.
        female.gc_progress = 0.0
        female.gc_threshold = 0.0
        female.gc_complete = false
        female.gc_start_day = 0
    end

    return (new_free, new_wolb, unhatched_eggs, density_dead_offspring)
end
