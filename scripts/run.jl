using Plots
include("main.jl")

use_burnin = false
simulationnumber = 1
weekly_release = true
zone = 3
K = 24600.0
burn_in_days = 364
calendar_days = 1095
initial_mosq_free = 10_000
burnin_dir = joinpath(dirname(@__DIR__), "data", "BurnIn")
burnin_environmental_data_path = joinpath(dirname(@__DIR__), "data", "ambiental16.csv")
environmental_data_path = joinpath(dirname(@__DIR__), "data", "ambiental.csv")

function field_release_schedule(P::MosquitoParameters)
    P.weekly_release || return Dict{Int64,Int64}()

    if P.zone == 1
        schedules = ((32, 52, 3500), (731, 22, 3500))
    elseif P.zone == 2
        schedules = ((152, 35, 22500), (731, 22, 22500))
    elseif P.zone == 3
        schedules = ((305, 39, 20600), (790, 22, 20600))
    else
        error("zone must be 1, 2, or 3")
    end

    release_schedule = Dict{Int64,Int64}()

    for (start_day, number_releases, release_size) in schedules
        for release_number in 0:(number_releases - 1)
            calendar_day = start_day + release_number * P.weekly_release_interval
            release_schedule[calendar_day] = release_size
        end
    end

    return release_schedule
end

if use_burnin
    burnin, burnin_path = load_burnin(burnin_dir, zone, K)
    P = posterior_parameters(burnin, calendar_days, weekly_release, environmental_data_path)
    release_schedule = field_release_schedule(P)
    result = main_from_burnin(simulationnumber, P, burnin; release_schedule = release_schedule)
else
    P = MosquitoParameters(
        initial_mosq_free = initial_mosq_free,
        zone = zone,
        burn_in_days = burn_in_days,
        calendar_days = calendar_days,
        weekly_release = weekly_release,
        K = K,
        burnin_environmental_data_path = burnin_environmental_data_path,
        environmental_data_path = environmental_data_path
    )
    release_schedule = field_release_schedule(P)
    result = main(simulationnumber, P; release_schedule = release_schedule)
end

println("Zone = $zone")
println("K = $K")
println("Calendar days = $calendar_days")
println("Burn-in = $use_burnin")
if use_burnin
    println("Burn-in file = $burnin_path")
    println("Burn-in seed = $(burnin.seed)")
    println("Burn-in days = $(burnin.burn_in_days)")
else
    println("Burn-in days = $burn_in_days")
end
println("Release events = $(length(release_schedule))")
println("Total released mosquitoes = $(sum(values(release_schedule)))")

days = use_burnin ? (0:calendar_days) : (-burn_in_days:calendar_days)
not_yet_eligible_female_ctr = result.not_yet_eligible_female_ctr
eligible_unmated_female_ctr = result.eligible_unmated_female_ctr
mated_female_ctr = result.mated_female_ctr
unmated_female_ctr = not_yet_eligible_female_ctr .+ eligible_unmated_female_ctr
unmated_fraction = [result.female_ctr[i] > 0 ? unmated_female_ctr[i] / result.female_ctr[i] : 0.0 for i in eachindex(result.female_ctr)]

fig_population = plot(
    days,
    result.free_ctr,
    label = "Wolbachia-free",
    xlabel = "Time (days)",
    ylabel = "Adult mosquito population",
    grid = false
)
plot!(fig_population, days, result.wolb_ctr, label = "Wolbachia-infected")

fig_sex_wolb = plot(
    days,
    result.female_free_ctr,
    label = "Female Wolbachia-free",
    xlabel = "Time (days)",
    ylabel = "Adult mosquito population",
    grid = false
)
plot!(fig_sex_wolb, days, result.female_wolb_ctr, label = "Female Wolbachia-infected")
plot!(fig_sex_wolb, days, result.male_free_ctr, label = "Male Wolbachia-free")
plot!(fig_sex_wolb, days, result.male_wolb_ctr, label = "Male Wolbachia-infected")

fig_mating = plot(
    days,
    result.female_ctr,
    label = "Total females",
    xlabel = "Time (days)",
    ylabel = "Female mosquito population",
    grid = false
)
plot!(fig_mating, days, mated_female_ctr, label = "Mated females")
plot!(fig_mating, days, unmated_female_ctr, label = "Unmated females")

fig_unmated_status = plot(
    days,
    not_yet_eligible_female_ctr,
    label = "Not yet eligible",
    xlabel = "Time (days)",
    ylabel = "Unmated female population",
    grid = false
)
plot!(fig_unmated_status, days, eligible_unmated_female_ctr, label = "Eligible unmated")

fig_unmated_fraction = plot(
    days,
    unmated_fraction,
    label = "Unmated fraction",
    xlabel = "Time (days)",
    ylabel = "Fraction of females unmated",
    ylim = (0, 1),
    grid = false
)

if !use_burnin
    for fig in (fig_population, fig_sex_wolb, fig_mating, fig_unmated_status, fig_unmated_fraction)
        vline!(fig, [0], label = "End of burn-in", linestyle = :dash)
    end
end

for fig in (fig_population, fig_sex_wolb, fig_mating, fig_unmated_status, fig_unmated_fraction)
    display(fig)
end

if use_burnin
    free_calendar = result.free_ctr[2:end]
    wolb_calendar = result.wolb_ctr[2:end]
else
    calendar_indices = (burn_in_days + 2):(burn_in_days + calendar_days + 1)
    free_calendar = result.free_ctr[calendar_indices]
    wolb_calendar = result.wolb_ctr[calendar_indices]
end

aeg_total = free_calendar .+ wolb_calendar
wmel_pc = [aeg_total[i] > 0 ? 100.0 * wolb_calendar[i] / aeg_total[i] : 0.0 for i in eachindex(aeg_total)]

fig_wmel = plot(
    1:calendar_days,
    wmel_pc,
    label = "wMel percentage",
    xlabel = "Calendar day",
    ylabel = "wMel-positive adult mosquitoes (%)",
    ylim = (0, 100),
    grid = false
)
display(fig_wmel)

month_days = Int64[
    31,28,31,30,31,30,31,31,30,31,30,31,
    31,28,31,30,31,30,31,31,30,31,30,31,
    31,28,31,30,31,30,31,31,30,31,30,31,
#    31,29,31
]
month_labels = ["$(lpad(mod(i - 1, 12) + 1, 2, '0'))/$(2017 + (i - 1) ÷ 12)" for i in eachindex(month_days)]

monthly_wmel_pc = Float64[]
first_day = 1

for number_days in month_days
    last_day = first_day + number_days - 1
    monthly_wolb = sum(wolb_calendar[first_day:last_day])
    monthly_total = sum(aeg_total[first_day:last_day])
    push!(monthly_wmel_pc, monthly_total > 0 ? 100.0 * monthly_wolb / monthly_total : 0.0)
    first_day = last_day + 1
end

months = eachindex(monthly_wmel_pc)
tick_positions = 1:3:length(month_labels)

fig_monthly_wmel = plot(
    months,
    monthly_wmel_pc,
    label = "Monthly wMel percentage",
    xlabel = "Month",
    ylabel = "wMel-positive adult mosquitoes (%)",
    xticks = (tick_positions, month_labels[tick_positions]),
    xrotation = 45,
    ylim = (0, 100),
    grid = false
)
display(fig_monthly_wmel)