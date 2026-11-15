# Figures for tuning curves, peri-event histograms, and population vectors.

const BLUE = "#2F6BFF"
const BLUE_DARK = "#183B70"
const GOLD = "#D99A20"
const INK = "#20252B"
const GRID = "#D9DEE7"
const LIGHT_BLUE = "#DCE8FF"
const LIGHT_GOLD = "#F7E6B7"
const MUTED = "#68717D"

"""Plot observed direction means ± SEM with the fitted sinusoidal curve."""
function plot_tuning_curve(
    summary::DataFrame;
    unit=nothing,
    epoch::Symbol=first(summary.epoch),
    size=(900, 560),
)
    fit = fit_tuning_curve(summary)
    curve_angles = collect(0.0:1.0:360.0)
    curve_rates = _predict_tuning(fit, curve_angles)
    trial_count = sum(summary.n_trials)
    unit_text = unit === nothing ? "" : " · unit $unit"
    title = "8-direction tuning curve, $(uppercase(String(epoch)))$unit_text"
    subtitle = "mean ± SEM · n=$trial_count trials · PD=$(round(fit.preferred_direction; digits=1))° · R²=$(round(fit.r_squared; digits=3))"

    plot = Plots.plot(
        curve_angles,
        curve_rates;
        color=BLUE,
        linewidth=3,
        label="sinusoidal fit",
        xlabel="Movement direction (degrees)",
        ylabel="Firing rate (spikes/s)",
        title="$title\n$subtitle",
        xlims=(0, 360),
        xticks=0:45:360,
        background_color=:white,
        foreground_color_subplot=INK,
        gridcolor=GRID,
        gridalpha=0.55,
        legend=:topright,
        size=size,
    )
    Plots.scatter!(
        plot,
        summary.direction,
        summary.mean_rate;
        yerror=summary.sem,
        color=GOLD,
        markerstrokecolor=INK,
        markerstrokewidth=1,
        markersize=7,
        label="direction mean ± SEM",
    )
    Plots.vline!(
        plot,
        [fit.preferred_direction];
        color=BLUE_DARK,
        linestyle=:dash,
        linewidth=1.8,
        label="preferred direction",
    )
    Plots.hline!(plot, [0.0]; color=MUTED, linewidth=1, label=false)
    return plot
end

"""Plot RT, MT, and TET tuning curves as linked vertical study panels."""
function plot_epoch_tuning(rate_table::DataFrame; unit=nothing, size=(900, 1450))
    panels = [
        plot_tuning_curve(
            direction_firing_rates(rate_table; epoch=epoch);
            unit=unit,
            epoch=epoch,
            size=(850, 450),
        )
        for epoch in (:rt, :mt, :tet)
    ]
    return Plots.plot(panels...; layout=(3, 1), size=size, link=:y, margin=5Plots.mm)
end

"""Peri-event histogram and spike raster in the style of Georgopoulos et al. (1982), Fig. 2."""
function plot_perievent(data; size=(1150, 760))
    widths = diff(data.edges)
    colors = [value >= 0 ? BLUE : LIGHT_BLUE for value in data.deviations]
    data_min = min(0.0, minimum(data.deviations))
    data_max = max(0.0, maximum(data.deviations))
    vertical_padding = max(0.14 * (data_max - data_min), 1.0)
    top = Plots.bar(
        data.centers .* 1000,
        data.deviations;
        bar_width=widths .* 1000,
        color=colors,
        linecolor=BLUE_DARK,
        linewidth=0.5,
        label=false,
        ylabel="Δ firing rate\n(spikes/s)",
        title="Movement-aligned activity: unit $(data.unit), $(data.direction)°\n20 ms bins · deviation from pre-target control ($(round(data.control_rate; digits=1)) spikes/s) · first 6 usable trials",
        background_color=:white,
        foreground_color_subplot=INK,
        gridcolor=GRID,
        gridalpha=0.5,
        xlims=(data.window[1] * 1000, data.window[2] * 1000),
        ylims=(data_min - 0.2 * vertical_padding, data_max + vertical_padding),
        legend=false,
    )
    Plots.hline!(top, [0.0]; color=INK, linewidth=1, label=false)
    Plots.vline!(top, [0.0]; color=INK, linestyle=:dot, linewidth=1.5, label=false)
    Plots.vspan!(
        top,
        [(data.target_mean - data.target_sd) * 1000, (data.target_mean + data.target_sd) * 1000];
        color=LIGHT_GOLD,
        alpha=0.45,
        label=false,
    )
    Plots.vline!(top, [data.target_mean * 1000]; color=GOLD, linestyle=:dash, linewidth=1.8, label=false)
    if data.initial_change !== missing
        Plots.vline!(
            top,
            [data.initial_change.time * 1000];
            color=BLUE_DARK,
            linestyle=:dashdot,
            linewidth=2,
            label=false,
        )
    end
    label_y = data_max + 0.55 * vertical_padding
    Plots.annotate!(top, 8, label_y, Plots.text("M", 9, INK, :left))
    Plots.annotate!(top, data.target_mean * 1000 + 8, label_y, Plots.text("T ± SD", 9, GOLD, :left))
    if data.initial_change !== missing
        Plots.annotate!(
            top,
            data.initial_change.time * 1000 + 8,
            label_y,
            Plots.text("Inc", 9, BLUE_DARK, :left),
        )
    end

    bottom = Plots.plot(
        ;
        xlabel="Time from movement onset (ms)",
        ylabel="Trial",
        yticks=(1:length(data.trials), string.(data.trials)),
        ylims=(0.4, length(data.trials) + 0.6),
        xlims=(data.window[1] * 1000, data.window[2] * 1000),
        background_color=:white,
        foreground_color_subplot=INK,
        grid=false,
        legend=false,
    )
    for (row, spikes) in enumerate(data.spike_rows)
        for spike in spikes
            x = spike * 1000
            Plots.plot!(bottom, [x, x], [row - 0.28, row + 0.28]; color=INK, linewidth=1.2, label=false)
        end
        target_x = data.target_times[row] * 1000
        Plots.plot!(bottom, [target_x, target_x], [row - 0.42, row + 0.42]; color=GOLD, linewidth=3, label=false)
    end
    Plots.vline!(bottom, [0.0]; color=INK, linestyle=:dot, linewidth=1.5, label=false)
    Plots.vline!(bottom, [data.target_mean * 1000]; color=GOLD, linestyle=:dash, linewidth=1.8, label=false)

    return Plots.plot(
        top,
        bottom;
        layout=Plots.grid(2, 1; heights=[0.48, 0.52]),
        link=:x,
        size=size,
        margin=5Plots.mm,
    )
end

function _direction_panel(data, common_ylims; show_x=false, show_y=false)
    widths = diff(data.edges)
    colors = [value >= 0 ? BLUE : LIGHT_BLUE for value in data.deviations]
    change_text = if data.initial_change === missing
        "Inc: none"
    else
        sign = data.initial_change.kind === :increase ? "↑" : "↓"
        "Inc $sign $(round(Int, data.initial_change.time * 1000)) ms"
    end
    panel = Plots.bar(
        data.centers .* 1000,
        data.deviations;
        bar_width=widths .* 1000,
        color=colors,
        linecolor=BLUE_DARK,
        linewidth=0.35,
        label=false,
        title="$(data.direction)° · $change_text",
        titlefontsize=10,
        xlabel=show_x ? "Time from M (ms)" : "",
        ylabel=show_y ? "Δ spikes/s" : "",
        xlims=(data.window[1] * 1000, data.window[2] * 1000),
        ylims=common_ylims,
        xticks=show_x ? ([-500, 0, 500], ["−500", "0", "500"]) : false,
        yticks=show_y ? :auto : false,
        background_color=:white,
        foreground_color_subplot=INK,
        gridcolor=GRID,
        gridalpha=0.35,
        legend=false,
    )
    Plots.hline!(panel, [0.0]; color=MUTED, linewidth=0.8, label=false)
    Plots.vline!(panel, [0.0]; color=INK, linestyle=:dot, linewidth=1.1, label=false)
    Plots.vline!(panel, [data.target_mean * 1000]; color=GOLD, linestyle=:dash, linewidth=1.2, label=false)
    if data.initial_change !== missing
        Plots.vline!(
            panel,
            [data.initial_change.time * 1000];
            color=BLUE_DARK,
            linestyle=:dashdot,
            linewidth=1.3,
            label=false,
        )
    end
    return panel
end

function _direction_compass()
    panel = Plots.plot(
        ;
        xlims=(-1.25, 1.25),
        ylims=(-1.25, 1.25),
        aspect_ratio=:equal,
        axis=false,
        grid=false,
        legend=false,
        background_color=:white,
        title="Spatial layout",
        titlefontsize=10,
    )
    for direction in 0:45:315
        x = 0.76 * cosd(direction)
        y = 0.76 * sind(direction)
        Plots.plot!(
            panel,
            [0.0, x],
            [0.0, y];
            color=direction in (0, 90, 180, 270) ? BLUE_DARK : MUTED,
            linewidth=1.8,
            arrow=true,
            label=false,
        )
        Plots.annotate!(
            panel,
            1.02 * cosd(direction),
            1.02 * sind(direction),
            Plots.text("$(direction)°", 8, INK),
        )
    end
    Plots.scatter!(panel, [0.0], [0.0]; color=GOLD, markerstrokecolor=INK, markersize=5, label=false)
    Plots.annotate!(panel, 0.0, -1.16, Plots.text("M = 0 ms · T = gold · Inc = blue", 7, MUTED))
    return panel
end

"""
    plot_directional_perievent_grid(session, unit; n_trials=6, kwargs...)

Arrange the eight 20 ms peri-event histograms according to target direction.
All panels use the same x/y ranges, the same number of chronologically first
usable trials, and the Georgopoulos initial-change rule.
"""
function plot_directional_perievent_grid(
    session::ReachSession,
    unit::Integer;
    n_trials::Integer=6,
    size=(1500, 1180),
    kwargs...,
)
    directions = collect(0:45:315)
    data = Dict(
        direction => peri_event_data(
            session,
            unit;
            direction=direction,
            n_trials=n_trials,
            kwargs...,
        )
        for direction in directions
    )
    all_deviations = reduce(vcat, [data[direction].deviations for direction in directions])
    data_min = min(0.0, minimum(all_deviations))
    data_max = max(0.0, maximum(all_deviations))
    padding = max(0.08 * (data_max - data_min), 1.0)
    common_ylims = (data_min - padding, data_max + padding)

    # Row-major placement preserves the physical target geometry.
    placement = [135, 90, 45, 180, :center, 0, 225, 270, 315]
    panels = Any[]
    for (position, direction) in enumerate(placement)
        if direction === :center
            push!(panels, _direction_compass())
        else
            row = cld(position, 3)
            column = mod1(position, 3)
            push!(
                panels,
                _direction_panel(
                    data[direction],
                    common_ylims;
                    show_x=row == 3,
                    show_y=column == 1,
                ),
            )
        end
    end
    return Plots.plot(
        panels...;
        layout=(3, 3),
        size=size,
        margin=3Plots.mm,
        plot_title="Peri-event histograms by direction: unit $unit\n20 ms bins · first $n_trials usable trials/direction · shared axes",
        plot_titlefontsize=15,
    )
end

"""R² histogram and preferred-direction histogram of tuned units."""
function plot_tuning_summary(tuning::DataFrame; r2_min=0.7, epoch=:tet)
    tuned = filter(row -> isfinite(row.r_squared) && row.r_squared >= r2_min, tuning)
    r2 = Plots.histogram(tuning.r_squared; bins=0:0.05:1, color="#9DB4D6", linecolor=:white,
                         legend=false, xlabel="R² of cosine fit", ylabel="units",
                         title="Directional tuning ($(uppercase(string(epoch))))")
    Plots.vline!(r2, [r2_min]; color="#C0392B", linestyle=:dash)
    pd = Plots.histogram(tuned.preferred_direction; bins=0:22.5:360, color=BLUE,
                         linecolor=:white, legend=false, xticks=0:45:360,
                         xlabel="preferred direction (°)", ylabel="units",
                         title="PD of tuned units (R² ≥ $r2_min, n=$(nrow(tuned)))")
    return Plots.plot(r2, pd; layout=(1, 2), size=(1000, 380), margin=5Plots.mm)
end

function _arrow!(plot, angle, len; kwargs...)
    Plots.plot!(plot, [0, len * cosd(angle)], [0, len * sind(angle)];
                arrow=Plots.arrow(:closed, :head, 0.4, 0.3), label="", kwargs...)
end

"""Mean actual movement (grey) and direction-averaged PV (blue) for each target."""
function plot_population_vectors(averaged::DataFrame, trials::DataFrame)
    plot = Plots.plot(; aspect_ratio=:equal, xlims=(-1.35, 1.35), ylims=(-1.35, 1.35),
                      axis=false, ticks=false, grid=false,
                      title="Direction-averaged population vectors (held-out trials)")
    for d in sort(unique(trials.direction))
        movement = circular_mean(trials.hand_angle[trials.direction .== d])
        decoded = circular_mean(averaged.pv_angle[averaged.direction .== d])
        _arrow!(plot, movement, 1.0; color="#B0B7C3", linewidth=5)
        _arrow!(plot, decoded, 0.8; color=BLUE, linewidth=2.5)
        Plots.annotate!(plot, 1.22 * cosd(d), 1.22 * sind(d), Plots.text("$(d)°", 9, MUTED))
    end
    Plots.plot!(plot, [NaN], [NaN]; color="#B0B7C3", linewidth=5, label="actual movement (mean)")
    Plots.plot!(plot, [NaN], [NaN]; color=BLUE, linewidth=2.5, label="population vector")
    return Plots.plot!(plot; legend=:outerbottom, legendcolumns=2, size=(600, 680))
end

"""Single-trial PV against actual movement, and the observed error against the shuffle null."""
function plot_single_trial_pv(single::DataFrame, null_errors::AbstractVector)
    scatter = Plots.scatter(single.hand_angle, single.pv_angle; legend=false,
                            markersize=4, alpha=0.6, color=BLUE, markerstrokewidth=0,
                            xlims=(0, 360), ylims=(0, 360), xticks=0:90:360, yticks=0:90:360,
                            aspect_ratio=:equal, xlabel="actual movement direction (°)",
                            ylabel="PV direction (°)", title="Single-trial PV (cross-validated)")
    Plots.plot!(scatter, [0, 360], [0, 360]; color="#B0B7C3", linestyle=:dash)
    observed = median(abs.(single.error_hand))
    null = Plots.histogram(null_errors; bins=30, color=GRID, linecolor=:white, legend=false,
                           xlabel="median |error| vs movement (°)", ylabel="shuffles",
                           title="Observed vs label-shuffle null")
    Plots.vline!(null, [observed]; color="#C0392B", linewidth=2)
    return Plots.plot(scatter, null; layout=(1, 2), size=(1050, 460), margin=5Plots.mm)
end
