# Single-unit directional tuning (Georgopoulos et al., 1982) for one unit:
# RT/MT/TET firing rates by direction, cosine fit, and movement-aligned
# raster/PSTH panels.
#
#   julia --project=. scripts/01_single_unit_tuning.jl [mat_path] [unit] [output_dir]

using DREAMReach
using DataFrames
using Plots

function main(args=ARGS)
    root = normpath(joinpath(@__DIR__, ".."))
    mat_path = length(args) >= 1 ? abspath(args[1]) : joinpath(root, "data", "Flint_2012_e1.mat")
    output_dir = length(args) >= 3 ? abspath(args[3]) : joinpath(root, "outputs", "single_unit")

    session = load_session(mat_path)
    candidates = candidate_sorted_units(session)
    isempty(candidates) && error("No candidate sorted units were found")
    unit = length(args) >= 2 ? parse(Int, args[2]) : first(candidates)
    mkpath(output_dir)

    # Georgopoulos et al. specified a velocity threshold crossing but did not
    # report its numeric value. These two parameters are therefore explicit.
    velocity_threshold = 0.05
    min_consecutive = 5
    rates = firing_rate_table(
        session,
        unit;
        velocity_threshold=velocity_threshold,
        min_consecutive=min_consecutive,
        velocity_mode=:target_projection,
    )

    summaries = DataFrame()
    for epoch in (:rt, :mt, :tet)
        append!(summaries, direction_firing_rates(rates; epoch=epoch), cols=:union)
    end
    tet = direction_firing_rates(rates; epoch=:tet)
    fit = fit_tuning_curve(tet)

    prefix = "unit_$(unit)"
    write_csv(joinpath(output_dir, "$(prefix)_trial_rates.csv"), rates)
    write_csv(joinpath(output_dir, "$(prefix)_direction_summary.csv"), summaries)

    tet_plot = plot_tuning_curve(tet; unit=unit, epoch=:tet)
    savefig(tet_plot, joinpath(output_dir, "$(prefix)_tuning_TET.png"))

    epoch_plot = plot_epoch_tuning(rates; unit=unit)
    savefig(epoch_plot, joinpath(output_dir, "$(prefix)_tuning_RT_MT_TET.png"))

    # Use the observed direction with the highest TET mean for the six-trial
    # paper-style example. Trials are the first six usable trials in time order.
    preferred_row = argmax(tet.mean_rate)
    raster_direction = tet.direction[preferred_row]
    peri = peri_event_data(
        session,
        unit;
        direction=raster_direction,
        n_trials=6,
        bin_width=0.02,
        velocity_threshold=velocity_threshold,
        min_consecutive=min_consecutive,
        velocity_mode=:target_projection,
    )
    peri_plot = plot_perievent(peri)
    savefig(peri_plot, joinpath(output_dir, "$(prefix)_figure2_raster_psth.png"))

    direction_grid = plot_directional_perievent_grid(
        session,
        unit;
        n_trials=6,
        bin_width=0.02,
        velocity_threshold=velocity_threshold,
        min_consecutive=min_consecutive,
        velocity_mode=:target_projection,
    )
    savefig(direction_grid, joinpath(output_dir, "$(prefix)_perievent_8_directions.png"))

    inc_text = peri.initial_change === missing ? "not detected" :
        "$(peri.initial_change.kind) at $(round(peri.initial_change.time * 1000; digits=1)) ms"
    println("Analysis complete")
    println("  MAT file: $mat_path")
    println("  unit: $unit")
    println("  usable trials: $(nrow(rates))")
    println("  velocity onset rule: target projection ≥ $velocity_threshold for $min_consecutive samples")
    println("  TET preferred direction: $(round(fit.preferred_direction; digits=1))°")
    println("  TET R²: $(round(fit.r_squared; digits=3))")
    println("  Figure 2 direction: $(raster_direction)°")
    println("  Inc: $inc_text")
    println("  eight-direction peri-event grid: $(prefix)_perievent_8_directions.png")
    println("  outputs: $output_dir")
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()
