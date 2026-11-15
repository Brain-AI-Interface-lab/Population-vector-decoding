# Movement-aligned raster/PSTH data and the Georgopoulos onset rule.

"""
Find the first three-bin sequence satisfying the Georgopoulos et al. rule.

All three bins must change on the same side of the control mean, and at least
two must differ from control by at least `fraction` (50% by default).
"""
function detect_initial_change(
    rates::AbstractVector{<:Real},
    control_mean::Real;
    times=nothing,
    search_start=(-Inf),
    consecutive::Integer=3,
    min_exceeding::Integer=2,
    fraction::Real=0.5,
)
    control_mean > 0 || return missing
    consecutive >= min_exceeding > 0 || throw(ArgumentError("Require consecutive ≥ min_exceeding > 0"))
    length(rates) >= consecutive || return missing
    centers = times === nothing ? collect(eachindex(rates)) : times
    length(centers) == length(rates) || throw(DimensionMismatch("times and rates must have equal length"))
    for index in 1:(length(rates)-consecutive+1)
        centers[index] < search_start && continue
        window = @view rates[index:(index+consecutive-1)]
        increasing = all(value -> value > control_mean, window)
        decreasing = all(value -> value < control_mean, window)
        large_increase = count(value -> value >= (1 + fraction) * control_mean, window)
        large_decrease = count(value -> value <= (1 - fraction) * control_mean, window)
        if (increasing && large_increase >= min_exceeding) ||
           (decreasing && large_decrease >= min_exceeding)
            return (
                index=index,
                time=centers[index],
                kind=increasing ? :increase : :decrease,
            )
        end
    end
    return missing
end

function _histogram_counts(spike_rows, edges)
    counts = zeros(Int, length(edges) - 1)
    for spikes in spike_rows, spike in spikes
        index = searchsortedlast(edges, spike)
        if 1 <= index < length(edges) && edges[index] <= spike < edges[index+1]
            counts[index] += 1
        end
    end
    return counts
end

"""Prepare six-trial movement-aligned raster/PSTH data for one direction."""
function peri_event_data(
    session::ReachSession,
    unit::Integer;
    direction::Integer,
    n_trials::Integer=6,
    window::Tuple{<:Real,<:Real}=(-0.8, 0.6),
    bin_width::Real=0.02,
    control_duration::Real=0.5,
    trials=good_trials(session),
    kwargs...,
)
    n_trials > 0 || throw(ArgumentError("n_trials must be positive"))
    bin_width > 0 || throw(ArgumentError("bin_width must be positive"))
    window[1] < window[2] || throw(ArgumentError("window start must precede end"))
    events = movement_event_table(session; trials=trials, kwargs...)
    chosen = events[events.direction .== mod(direction, 360), :]
    nrow(chosen) >= n_trials || error("Direction $direction° has only $(nrow(chosen)) usable trials")
    chosen = chosen[1:n_trials, :]

    spike_rows = Vector{Vector{Float64}}()
    target_times = Float64[]
    control_spikes = 0
    control_seconds = 0.0
    for event in eachrow(chosen)
        absolute = spike_times(session, event.trial, unit; relative_to=:absolute)
        relative = absolute .- event.movement_onset
        push!(spike_rows, filter(time -> window[1] <= time <= window[2], relative))
        push!(target_times, event.target_onset - event.movement_onset)
        control_start = max(first(session.times[event.trial]), event.target_onset - control_duration)
        control_spikes += _spike_count(absolute, control_start, event.target_onset)
        control_seconds += event.target_onset - control_start
    end

    edges = collect(Float64(window[1]):Float64(bin_width):Float64(window[2]))
    last(edges) < window[2] - eps(Float64) && push!(edges, Float64(window[2]))
    counts = _histogram_counts(spike_rows, edges)
    widths = diff(edges)
    rates = counts ./ (n_trials .* widths)
    control_rate = control_spikes / control_seconds
    centers = edges[1:(end-1)] .+ widths ./ 2
    deviations = rates .- control_rate
    target_mean = mean(target_times)
    target_sd = length(target_times) > 1 ? std(target_times) : 0.0
    change = detect_initial_change(
        rates,
        control_rate;
        times=centers,
        search_start=target_mean,
    )

    return (
        unit=unit,
        direction=mod(direction, 360),
        trials=Int.(chosen.trial),
        spike_rows=spike_rows,
        target_times=target_times,
        target_mean=target_mean,
        target_sd=target_sd,
        edges=edges,
        centers=centers,
        counts=counts,
        rates=rates,
        control_rate=control_rate,
        deviations=deviations,
        initial_change=change,
        bin_width=Float64(bin_width),
        window=(Float64(window[1]), Float64(window[2])),
    )
end
