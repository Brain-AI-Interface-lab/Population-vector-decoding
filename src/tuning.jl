# Epoch firing rates, per-direction means, and the cosine tuning fit.

_spike_count(spikes, start_time, end_time) = count(time -> start_time <= time < end_time, spikes)

"""
    firing_rate_table(session, unit; control_duration=0.5, kwargs...)

Return trial-level firing rates for the control, reaction-time (RT),
movement-time (MT), and total-experimental-time (TET) epochs.
"""
function firing_rate_table(
    session::ReachSession,
    unit::Integer;
    trials=good_trials(session),
    control_duration::Real=0.5,
    kwargs...,
)
    control_duration > 0 || throw(ArgumentError("control_duration must be positive"))
    events = movement_event_table(session; trials=trials, kwargs...)
    rows = NamedTuple[]
    for event in eachrow(events)
        spikes = spike_times(session, event.trial, unit; relative_to=:absolute)
        control_start = max(first(session.times[event.trial]), event.target_onset - control_duration)
        control_length = event.target_onset - control_start
        rt_count = _spike_count(spikes, event.target_onset, event.movement_onset)
        mt_count = _spike_count(spikes, event.movement_onset, event.movement_end)
        tet_count = _spike_count(spikes, event.target_onset, event.movement_end)
        control_count = _spike_count(spikes, control_start, event.target_onset)
        push!(rows, (
            trial=event.trial,
            direction=event.direction,
            target_onset=event.target_onset,
            acceleration_onset=event.acceleration_onset,
            movement_onset=event.movement_onset,
            movement_end=event.movement_end,
            reaction_time=event.reaction_time,
            movement_time=event.movement_time,
            total_experimental_time=event.total_experimental_time,
            control_count=control_count,
            rt_count=rt_count,
            mt_count=mt_count,
            tet_count=tet_count,
            control_rate=control_count / control_length,
            rt_rate=rt_count / event.reaction_time,
            mt_rate=mt_count / event.movement_time,
            tet_rate=tet_count / event.total_experimental_time,
        ))
    end
    return DataFrame(rows)
end

function _epoch_column(epoch::Symbol)
    epoch in (:control, :rt, :mt, :tet) ||
        throw(ArgumentError("epoch must be :control, :rt, :mt, or :tet"))
    return Symbol(epoch, "_rate")
end

"""Summarize mean firing rate, SD, and SEM for every reach direction."""
function direction_firing_rates(rate_table::DataFrame; epoch::Symbol=:tet)
    column = _epoch_column(epoch)
    rows = NamedTuple[]
    for direction in sort(unique(rate_table.direction))
        values = Float64.(rate_table[rate_table.direction .== direction, column])
        n = length(values)
        spread = n > 1 ? std(values) : NaN
        push!(rows, (
            direction=Int(direction),
            n_trials=n,
            mean_rate=mean(values),
            sd=spread,
            sem=spread / sqrt(n),
            epoch=epoch,
        ))
    end
    return DataFrame(rows)
end


"""Fit `y = b0 + b1*sin(theta) + b2*cos(theta)` to direction means."""
function fit_tuning_curve(summary::DataFrame)
    nrow(summary) >= 3 || throw(ArgumentError("At least three directions are required"))
    angles = Float64.(summary.direction)
    rates = Float64.(summary.mean_rate)
    design = hcat(ones(length(angles)), sind.(angles), cosd.(angles))
    b0, b1, b2 = design \ rates
    fitted = design * [b0, b1, b2]
    residual_sum = sum(abs2, rates .- fitted)
    total_sum = sum(abs2, rates .- mean(rates))
    r_squared = total_sum > 0 ? 1 - residual_sum / total_sum : NaN
    amplitude = hypot(b1, b2)
    preferred_direction = mod(rad2deg(atan(b1, b2)), 360.0)
    modulation_index = b0 > 0 ? amplitude / b0 : NaN
    return (
        b0=b0,
        b1=b1,
        b2=b2,
        amplitude=amplitude,
        preferred_direction=preferred_direction,
        r_squared=r_squared,
        modulation_index=modulation_index,
    )
end


function _predict_tuning(fit, angles)
    return fit.b0 .+ fit.b1 .* sind.(angles) .+ fit.b2 .* cosd.(angles)
end
