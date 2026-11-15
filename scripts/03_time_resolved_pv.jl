# Time-resolved population vector around movement onset.
#
#   julia --project=. scripts/03_time_resolved_pv.jl [mat_path] [output_dir]
#
# Rates are counted in a 100 ms window sliding in 20 ms steps from 400 ms
# before to 600 ms after movement onset. Tuning comes from TET rates of the
# training folds; population vectors are computed from held-out trials.

using DREAMReach
using DataFrames
using Random
using Statistics
import Plots

const R2_MIN = 0.7
const N_FOLDS = 5
const SEED = 20261116
const TIMES = -0.4:0.02:0.6
const WINDOW = 0.1

function main(args=ARGS)
    root = normpath(joinpath(@__DIR__, ".."))
    mat_path = length(args) >= 1 ? abspath(args[1]) : joinpath(root, "data", "Flint_2012_e1.mat")
    out = length(args) >= 2 ? abspath(args[2]) : joinpath(root, "outputs", "time_resolved")
    mkpath(out)

    session = load_session(mat_path)
    units = candidate_sorted_units(session)
    trials = decoding_trials(session)
    tet_rates = epoch_rates(session, trials, units; epoch=:tet)
    windowed = windowed_rates(session, trials, units; times=TIMES, window=WINDOW)

    # Same folds as 02_population_vector.jl (same seed, same trial order).
    folds = stratified_folds(trials.direction, N_FOLDS, MersenneTwister(SEED))
    pv = time_resolved_pv(windowed, tet_rates, trials, folds; times=TIMES, r2_min=R2_MIN)
    speed = hand_speed_by_direction(session, trials; times=TIMES)

    write_csv(joinpath(out, "pv_time_resolved.csv"), pv)
    write_csv(joinpath(out, "hand_speed.csv"), speed)
    Plots.default(fontfamily="Helvetica", framestyle=:box, grid=false)
    Plots.savefig(plot_time_resolved_pv(pv, speed), joinpath(out, "pv_time_course.png"))
    animate_time_resolved_pv(pv, speed; path=joinpath(out, "pv_time_resolved.gif"))

    times = collect(TIMES)
    pv_curve = [mean(pv.length[pv.time .== t]) for t in times]
    speed_curve = [mean(speed.speed[speed.time .== t]) for t in times]
    pv_half = time_to_fraction(times, pv_curve)
    speed_half = time_to_fraction(times, speed_curve)
    error_at(t) = median(abs.(pv.error_target[isapprox.(pv.time, t; atol=1e-9)]))

    summary = """
    Time-resolved population vector, $(basename(mat_path))
    window $(round(Int, 1000WINDOW)) ms, step $(round(Int, 1000step(TIMES))) ms, $(round(Int, 1000first(TIMES))) to $(round(Int, 1000last(TIMES))) ms from movement onset

    PV length reaches half of its range at   $(round(Int, 1000pv_half)) ms
    hand speed reaches half of its range at  $(round(Int, 1000speed_half)) ms
    PV leads hand speed by                    $(round(Int, 1000(speed_half - pv_half))) ms

    median |PV error| vs target at -300 ms:  $(round(error_at(-0.3); digits=1))°
                                  at -100 ms: $(round(error_at(-0.1); digits=1))°
                                  at    0 ms: $(round(error_at(0.0); digits=1))°
                                  at +200 ms: $(round(error_at(0.2); digits=1))°
    """
    write(joinpath(out, "summary.txt"), summary)
    print(summary)
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()
