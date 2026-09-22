@with_kw struct MosquitoParameters
    initial_mosq_free::Int64

    zone::Int64
    burn_in_days::Int64 = 364
    calendar_days::Int64 = 1186
    weekly_release_interval::Int64 = 7
    weekly_release::Bool = true

    natural_female_probability::Float64 = 0.5
    release_female_probability::Float64 = 0.5   #3.0 / 7.0

    burnin_environmental_data_path::String = joinpath(dirname(@__DIR__), "data", "ambiental16.csv")
    environmental_data_path::String = joinpath(dirname(@__DIR__), "data", "ambiental.csv")

    # adult mosquito survival: Doeurk et al. (2025)
    # beta_s(T) = exp(a_s + b_s*Tmean + c_s*Tmean^2)

    # females
    weibull_a_fem::Float64 = 6.78985436
    weibull_b_fem::Float64 = -0.92326243
    weibull_c_fem::Float64 = 0.01953523
    weibull_k_fem::Float64 = 1.0478

    # males
    weibull_a_male::Float64 = 3.13665639
    weibull_b_male::Float64 = -0.53812540
    weibull_c_male::Float64 = 0.01141217
    weibull_k_male::Float64 = 1.2055

    # laboratory-to-field lifespan rescaling
    # m_adj(T) = field_lifespan_multiplier * m_lab(T)
    field_lifespan_multiplier::Float64 = 0.14

    wolb_lifespan_multiplier::Float64 = 1.0

    delta::Float64 = 0.95
#    eta::Float64 = 0.95

    K::Float64 = 0.0
    s0::Float64 = 0.87
    q::Float64 = 1.0 # q=1:contest, q>1:scramble

    # temperature-dependent immature development: Doeurk et al. (2025)
    # Weibull developmental requirement with Logan-10 thermal accumulation
    immature_weibull_k::Float64 = 6.21738199
    immature_logan_a::Float64 = 9.07832e-5
    immature_logan_b::Float64 = 0.15839717
    immature_Tmax::Float64 = 40.0

    # females become available for mating after emergence
    mating_age_min::Int64 = 2
    mating_age_max::Int64 = 2

    # daily mating probability for an eligible unmated female
    daily_mating_probability::Float64 = 0.5

    # probability that an eligible mated female obtains a blood meal on a given day
    bite_probability::Float64 = 0.5

    # Carrington
    gc_sigma::Float64 = 0.321
    gc_briere_c::Float64 = 4.0127e-4
    gc_briere_T0::Float64 = 15.8402
    gc_briere_Tm::Float64 = 35.7721
    max_gonotrophic_cycles::Int64 = 3

    # temperature-dependent fecundity: generalized Briere response
    # normalized so that phi_E(Tref) = 1 at the Osorio reference temperature
    fecundity_briere_T0::Float64 = 15.54
    fecundity_briere_Tm::Float64 = 37.39
    fecundity_briere_b::Float64 = 1.35
    fecundity_reference_temperature::Float64 = 27.0

    # fecundity and fertility decrease across gonotrophic cycles
    fecundity_scaling::NTuple{3, Float64} = (1.0, 0.896, 0.590)
    fertility_scaling::NTuple{3, Float64} = (1.0, 0.818, 0.494)

    # cross notation: f = Wolbachia-free, w = Wolbachia-infected
    eggs_r_ff::Float64 = 92.83
    eggs_mu_ff::Float64 = 96.57

    eggs_r_fw::Float64 = 29.03
    eggs_mu_fw::Float64 = 93.91

    eggs_r_wf::Float64 = 10.75
    eggs_mu_wf::Float64 = 81.66

    eggs_r_ww::Float64 = 2.18
    eggs_mu_ww::Float64 = 64.93

    hatch_ff::Float64 = 0.884
    hatch_fw::Float64 = 0.006
    hatch_wf::Float64 = 0.683
    hatch_ww::Float64 = 0.622
end
