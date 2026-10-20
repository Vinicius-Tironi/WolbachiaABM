using Distributed
using Statistics
using Serialization
using CSV

include("main.jl")


const MONTH_DAYS = Int64[
    31,28,31,30,31,30,31,31,30,31,30,31,
    31,28,31,30,31,30,31,31,30,31,30,31,
    31,28,31,30,31,30,31,31,30,31,30,31,
    31,29,31
]

const MONTH_LABELS = vcat(
    ["2017-$(lpad(m, 2, '0'))" for m in 1:12],
    ["2018-$(lpad(m, 2, '0'))" for m in 1:12],
    ["2019-$(lpad(m, 2, '0'))" for m in 1:12],
    ["2020-$(lpad(m, 2, '0'))" for m in 1:3]
)

month_string(x) = first(replace(string(x), "_" => "-"), 7)
anscombe(y::Real, n::Real) = asin(sqrt((y + 3.0 / 8.0) / (n + 3.0 / 4.0)))

# HISTORICAL FIELD RELEASE SCHEDULE
function field_release_schedule(P::MosquitoParameters)

    P.weekly_release || return Dict{Int64,Int64}()

    if P.zone == 1

        schedules = (
            (32, 52, 3500),
            (731, 22, 3500)
        )

    elseif P.zone == 2

        schedules = (
            (152, 35, 22500),
            (731, 22, 22500)
        )

    elseif P.zone == 3

        schedules = (
            (305, 39, 20600),
            (790, 22, 20600)
        )

    else

        error("zone must be 1, 2, or 3")

    end

    release_schedule = Dict{Int64,Int64}()

    for (start_day, number_releases, release_size) in schedules

        for release_number in 0:(number_releases - 1)

            calendar_day =
                start_day +
                release_number * P.weekly_release_interval

            release_schedule[calendar_day] = release_size

        end
    end

    return release_schedule
end


function load_observations(path::String, zone::Int64)
    rows = [(month = month_string(row.year_month), total = Int64(round(row.aeg_total)), wmel = Int64(round(row.aeg_wmel))) for row in CSV.File(path) if Int64(row.zone) == zone && row.aeg_total > 0]
    sort!(rows, by = row -> row.month)

    return (
        month = [row.month for row in rows],
        total = [row.total for row in rows],
        wmel = [row.wmel for row in rows]
    )
end


function monthly_prevalence(free_ctr::Vector{Int64}, wolb_ctr::Vector{Int64}, burn_in_days::Int64, calendar_days::Int64)
    first_index = burn_in_days + 2
    last_index = burn_in_days + calendar_days + 1

    free = free_ctr[first_index:last_index]
    wolb = wolb_ctr[first_index:last_index]

    prevalence = fill(NaN, length(MONTH_DAYS))
    first_day = 1

    for m in eachindex(MONTH_DAYS)
        last_day = first_day + MONTH_DAYS[m] - 1
        monthly_free = sum(free[first_day:last_day])
        monthly_wolb = sum(wolb[first_day:last_day])
        monthly_total = monthly_free + monthly_wolb

        if monthly_total > 0
            prevalence[m] = monthly_wolb / monthly_total
        end

        first_day = last_day + 1
    end

    return prevalence
end


function calculate_distance(simulated::Vector{Float64}, observed_transformed::Vector{Float64}, observed_totals::Vector{Int64}, weights::Vector{Float64})
    any(isnan, simulated) && return Inf

    simulated_transformed = [anscombe(observed_totals[m] * simulated[m], observed_totals[m]) for m in eachindex(simulated)]

    return sum(weights .* (observed_transformed .- simulated_transformed).^2) / sum(weights)
end


function save_burnin(result, P::MosquitoParameters, seed::Int64, burnin_dir::String)
    result.burnin_state === nothing && error("Burn-in snapshot was not generated")

    zone_dir = joinpath(burnin_dir, "zone$(P.zone)")
    mkpath(zone_dir)

    path = joinpath(zone_dir, "K$(Int(round(P.K))).jls")

    burnin_tmean = load_tmean_period(P.burnin_environmental_data_path, P.zone, P.burn_in_days; from_end = true)

    state = (
        format_version = "absolute_clock_v1",
        zone = P.zone,
        K = P.K,
        burn_in_days = P.burn_in_days,
        absolute_day = P.burn_in_days,
        seed = seed,
        parameters = P,
        burnin_tmean = burnin_tmean,
        mosq_fem = deepcopy(result.burnin_state.mosq_fem),
        mosq_male = deepcopy(result.burnin_state.mosq_male),
        immature_free = deepcopy(result.burnin_state.immature_free),
        immature_wolb = deepcopy(result.burnin_state.immature_wolb)
    )

    serialize(path, state)

    return path
end


function evaluate_K(
    K::Float64,
    zone::Int64,
    seed::Int64,
    delta::Float64,
    wolb_lifespan_multiplier::Float64,
    burn_in_days::Int64,
    calendar_days::Int64,
    initial_mosq_free::Int64,
    observations,
    observed_transformed::Vector{Float64},
    weights::Vector{Float64},
    month_indices::Vector{Int64},
    burnin_environmental_data_path::String,
    environmental_data_path::String,
    burnin_dir::String; save_selected_burnin::Bool = false
)

    P = MosquitoParameters(
        initial_mosq_free = initial_mosq_free,
        zone = zone,
        burn_in_days = burn_in_days,
        calendar_days = calendar_days,
        weekly_release = true,
        natural_female_probability = 0.5,
        release_female_probability = 3.0 / 7.0,
        delta = delta,
        wolb_lifespan_multiplier = wolb_lifespan_multiplier,
        K = K,
        burnin_environmental_data_path = burnin_environmental_data_path,
        environmental_data_path = environmental_data_path
    )

    release_schedule = field_release_schedule(P)

    result = main(
        seed,
        P;
        release_schedule = release_schedule
    )
    burnin_path = save_selected_burnin ? save_burnin(result, P, seed, burnin_dir) : "not_saved"

    monthly_simulated = monthly_prevalence(result.free_ctr, result.wolb_ctr, burn_in_days, calendar_days)
    simulated = [monthly_simulated[index] for index in month_indices]
    Q = calculate_distance(simulated, observed_transformed, observations.total, weights)

    burnin_index = burn_in_days + 1

    return (
        Q = Q,
        simulated = simulated,
        monthly_simulated = monthly_simulated,
        burnin_path = burnin_path,
        burnin_free = result.free_ctr[burnin_index],
        burnin_females = result.female_ctr[burnin_index],
        burnin_males = result.male_ctr[burnin_index],
        burnin_immatures = total_immature_population(result.burnin_state.immature_free),
        final_free = result.free_ctr[end],
        final_wolb = result.wolb_ctr[end],
        number_release_events = length(release_schedule),
        total_released = sum(values(release_schedule)),
    )
end



# Replicated M1 K calibration. One SLURM job for the full K grid; Julia workers parallelize seeds.
# Uses the original historical release design, Q objective, and observation weights.
# Existing calibration, burn-in snapshots and other output directories are untouched.

if myid() == 1
    const MODEL_ROOT = dirname(@__DIR__)
    const PROJECT_ROOT = MODEL_ROOT
    const OUTPUT_DIR = joinpath(MODEL_ROOT, "outputs", "SMD_Calibration")
    const K_VALUES = collect(24_500:100:25_500)
    const R = 500
    const FIRST_SEED = 20_001

    function aggregate_results!()
        rows = NamedTuple[]
        for K in K_VALUES
            path = joinpath(OUTPUT_DIR, "K$(K)_replicates.csv")
            isfile(path) || error("Missing results: $path")
            data = collect(CSV.File(path))
            length(data) == R || error("Expected $R replicates for K=$K, found $(length(data))")
            seeds = Int.(getproperty.(data, :seed))
            sort(seeds) == collect(FIRST_SEED:(FIRST_SEED+R-1)) || error("Incomplete seed set for K=$K")
            q = Float64.(getproperty.(data, :Q))
            all(isfinite, q) || error("Nonfinite Q for K=$K; examine $path")
            push!(rows, (K=K, R=R, mean_Q=mean(q), sd_Q=std(q), se_Q=std(q)/sqrt(R), min_Q=minimum(q), max_Q=maximum(q)))
        end
        sort!(rows, by=x -> x.K)
        CSV.write(joinpath(OUTPUT_DIR, "K_grid_R500.csv"), rows)
        ranked = sort(copy(rows), by=x -> x.mean_Q)
        best = first(ranked)
        open(joinpath(OUTPUT_DIR, "summary_R500.txt"), "w") do io
            println(io, "Model: M1")
            println(io, "K grid: 24500:100:25500")
            println(io, "Replicates per K: $R")
            println(io, "Seeds: $FIRST_SEED:$(FIRST_SEED+R-1)")
            println(io, "Selected K (minimum mean Q): $(best.K)")
            println(io, "Mean Q: $(best.mean_Q)")
            println(io, "SE of Q: $(best.se_Q)")
            println(io, "Boundary optimum: $(best.K == first(K_VALUES) || best.K == last(K_VALUES))")
            println(io, "No burn-in snapshots were created during the K-grid evaluations.")
        end
        println("Best K for M1 (R=500): $(best.K) | mean Q=$(best.mean_Q) | SE=$(best.se_Q)")
        return best.K
    end

    # Aggregate-only mode is retained for postprocessing if needed.
    if length(ARGS) == 1 && ARGS[1] == "aggregate"
        aggregate_results!()
    else
        isempty(ARGS) || error("Usage: julia -p <workers> $(basename(@__FILE__)) [aggregate]")

        monitoring_path = get(ENV, "MONITORING_PATH", joinpath(MODEL_ROOT, "data", "monthly_anders.csv"))
        burnin_environmental_data_path = joinpath(MODEL_ROOT, "data", "ambiental16.csv")
        environmental_data_path = joinpath(MODEL_ROOT, "data", "ambiental.csv")
        for path in (monitoring_path, burnin_environmental_data_path, environmental_data_path)
            isfile(path) || error("Required input missing: $path")
        end

        mkpath(OUTPUT_DIR)
        observations = load_observations(monitoring_path, 3)
        observed_csv = joinpath(OUTPUT_DIR, "observed_prevalence.csv")
        if !isfile(observed_csv)
            observed_rows = [
                (month=observations.month[i], observed_total=observations.total[i],
                 observed_wmel=observations.wmel[i],
                 observed_prevalence=observations.wmel[i] / observations.total[i])
                for i in eachindex(observations.month)
            ]
            temp_observed = observed_csv * ".partial.$(getpid())"
            CSV.write(temp_observed, observed_rows)
            mv(temp_observed, observed_csv; force=false)
        end
        month_index = Dict(month => i for (i, month) in enumerate(MONTH_LABELS))
        all(haskey(month_index, month) for month in observations.month) || error("Monitoring months outside model calendar")
        month_indices = Int64[month_index[month] for month in observations.month]
        observed_transformed = [anscombe(observations.wmel[m], observations.total[m]) for m in eachindex(observations.month)]
        weights = Float64.(observations.total) .+ 0.5

        # Load functions on workers; the myid()==1 guard prevents re-entering the driver.
        for pid in workers()
            remotecall_wait(Base.include, pid, Main, abspath(@__FILE__))
        end

        println("M1 | $(length(K_VALUES)) K candidates | R=$R | workers=$(nworkers())")
        for K_integer in K_VALUES
            output_csv = joinpath(OUTPUT_DIR, "K$(K_integer)_replicates.csv")
            trajectory_csv = joinpath(OUTPUT_DIR, "K$(K_integer)_monthly_trajectories.csv")
            if isfile(output_csv) || isfile(trajectory_csv)
                # Resume only when both files are complete; never overwrite partial results.
                isfile(output_csv) && isfile(trajectory_csv) ||
                    error("Incomplete output pair for K=$K_integer; inspect files before resuming")
                data = collect(CSV.File(output_csv))
                seeds = sort(Int.(getproperty.(data, :seed)))
                q = Float64.(getproperty.(data, :Q))
                length(data) == R && seeds == collect(FIRST_SEED:(FIRST_SEED+R-1)) && all(isfinite, q) ||
                    error("Existing replicate results invalid: $output_csv")
                trajectory_rows = collect(CSV.File(trajectory_csv))
                length(trajectory_rows) == R * length(MONTH_LABELS) ||
                    error("Existing monthly trajectories incomplete: $trajectory_csv")
                all(Int(row.K) == K_integer for row in trajectory_rows) ||
                    error("Incorrect K in $trajectory_csv")
                counts = Dict(seed => 0 for seed in FIRST_SEED:(FIRST_SEED+R-1))
                for row in trajectory_rows
                    haskey(counts, Int(row.seed)) || error("Unexpected seed in $trajectory_csv")
                    counts[Int(row.seed)] += 1
                end
                all(count == length(MONTH_LABELS) for count in values(counts)) ||
                    error("Incomplete seed trajectories in $trajectory_csv")
                println("Already complete K=$K_integer | saved=$output_csv and $trajectory_csv")
                continue
            end

            K = Float64(K_integer)
            println("Starting K=$K_integer | R=$R")
            results = pmap(collect(FIRST_SEED:(FIRST_SEED+R-1))) do seed
                evaluation = evaluate_K(
                    K, 3, Int64(seed), 0.95, 1.0, Int64(364), Int64(1186), Int64(100_000),
                    observations, observed_transformed, weights, month_indices,
                    burnin_environmental_data_path, environmental_data_path, ""
                )
                replicate = (model="m1", K=Int(K), seed=Int(seed), Q=evaluation.Q,
                             burnin_free=evaluation.burnin_free, burnin_females=evaluation.burnin_females,
                             burnin_males=evaluation.burnin_males, burnin_immatures=evaluation.burnin_immatures)
                return (replicate=replicate, monthly_simulated=evaluation.monthly_simulated)
            end
            sort!(results, by=x -> x.replicate.seed)
            rows = [x.replicate for x in results]
            all(isfinite(x.Q) for x in rows) || error("Nonfinite Q at K=$K_integer; not saving")
            all(length(x.monthly_simulated) == length(MONTH_LABELS) for x in results) ||
                error("Unexpected monthly trajectory length for K=$K_integer")
            trajectory_rows = [
                (K=K_integer, seed=x.replicate.seed, month=MONTH_LABELS[m],
                 simulated_prevalence=x.monthly_simulated[m])
                for x in results for m in eachindex(MONTH_LABELS)
            ]
            temp_replicates = output_csv * ".partial.$(getpid())"
            temp_trajectories = trajectory_csv * ".partial.$(getpid())"
            CSV.write(temp_replicates, rows)
            CSV.write(temp_trajectories, trajectory_rows)
            mv(temp_trajectories, trajectory_csv; force=false)
            mv(temp_replicates, output_csv; force=false)
            q = [x.Q for x in rows]
            println("Completed K=$K_integer | mean Q=$(mean(q)) | sd Q=$(std(q))")
        end

        # Combine per-K monthly trajectories into a single plotting CSV.
        combined_csv = joinpath(OUTPUT_DIR, "monthly_trajectories_R500.csv")
        if !isfile(combined_csv)
            temp_combined = combined_csv * ".partial.$(getpid())"
            open(temp_combined, "w") do io
                for (i, K) in enumerate(K_VALUES)
                    path = joinpath(OUTPUT_DIR, "K$(K)_monthly_trajectories.csv")
                    open(path, "r") do input
                        for (j, line) in enumerate(eachline(input))
                            if i == 1 || j > 1
                                println(io, line)
                            end
                        end
                    end
                end
            end
            mv(temp_combined, combined_csv; force=false)
        end
        println("Monthly trajectories: $combined_csv")
        selected_K = aggregate_results!()
        selected_burnin = evaluate_K(
            Float64(selected_K), 3, Int64(FIRST_SEED), 0.95, 1.0,
            Int64(364), Int64(1186), Int64(100_000),
            observations, observed_transformed, weights, month_indices,
            burnin_environmental_data_path, environmental_data_path,
            joinpath(OUTPUT_DIR, "BurnIn"); save_selected_burnin = true
        )
        println("Selected burn-in snapshot: $(selected_burnin.burnin_path)")
    end
end
