# Target onset, movement onset, and movement end for each trial.

function _first_sustained(values, first_index::Integer, run_length::Integer, predicate)
    run_length > 0 || throw(ArgumentError("run_length must be positive"))
    last_start = length(values) - run_length + 1
    first_index > last_start && return nothing
    for index in first_index:last_start
        all(predicate, @view values[index:(index+run_length-1)]) && return index
    end
    return nothing
end

function _target_onset_index(
    session::ReachSession,
    trial::Integer;
    dims::NTuple{2,Int}=DEFAULT_DIMS,
    tolerance::Real=1e-8,
)
    trial = _check_trial(session, trial)
    target = trial_target(session, trial; dims=dims, tolerance=tolerance)
    (target === missing || target.outer === missing) && return nothing
    positions = session.target_pos[trial]
    count = min(length(session.times[trial]), size(positions, 1))
    for index in 1:count
        point = (positions[index, dims[1]], positions[index, dims[2]])
        if all(isfinite, point) && hypot(point[1] - target.center[1], point[2] - target.center[2]) > tolerance
            return index
        end
    end
    return nothing
end

"""Return the first appearance time of the outer target, or `missing`."""
function target_onset(
    session::ReachSession,
    trial::Integer;
    dims::NTuple{2,Int}=DEFAULT_DIMS,
)
    index = _target_onset_index(session, trial; dims=dims)
    return index === nothing ? missing : session.times[trial][index]
end

function _planar_magnitude(matrix::AbstractMatrix, dims::NTuple{2,Int})
    _check_dims(dims, matrix)
    return hypot.(matrix[:, dims[1]], matrix[:, dims[2]])
end

function _target_projection(
    matrix::AbstractMatrix,
    center::NTuple{2,<:Real},
    outer::NTuple{2,<:Real},
    dims::NTuple{2,Int},
)
    _check_dims(dims, matrix)
    dx = outer[1] - center[1]
    dy = outer[2] - center[2]
    norm = hypot(dx, dy)
    norm > 0 || error("Center and outer target are identical")
    ux, uy = dx / norm, dy / norm
    return matrix[:, dims[1]] .* ux .+ matrix[:, dims[2]] .* uy
end

function _acceleration_candidate(
    session::ReachSession,
    trial::Integer,
    target_time::Real,
    target_index::Integer,
    movement_index::Integer;
    dims::NTuple{2,Int}=DEFAULT_DIMS,
    baseline_duration::Real=0.3,
    threshold_sd::Real=3.0,
    min_consecutive::Integer=2,
)
    time = session.times[trial]
    acceleration = _planar_magnitude(session.hand_acc[trial], dims)
    count = min(length(time), length(acceleration))
    baseline = [
        acceleration[index] for index in 1:count
                                if target_time - baseline_duration <= time[index] < target_time
    ]
    isempty(baseline) && return (onset=missing, threshold=NaN)
    center = median(baseline)
    robust_sd = 1.4826 * median(abs.(baseline .- center))
    robust_sd = robust_sd > 0 ? robust_sd : std(baseline)
    threshold = center + threshold_sd * max(robust_sd, eps(Float64))
    stop_index = min(movement_index, count)
    candidate = _first_sustained(
        @view(acceleration[1:stop_index]),
        target_index,
        min_consecutive,
        value -> value >= threshold,
    )
    onset = candidate === nothing ? missing : time[candidate]
    return (onset=onset, threshold=threshold)
end

"""
    trial_events(session, trial; velocity_threshold=0.05, min_consecutive=5,
                 velocity_mode=:target_projection)

Detect target appearance, movement onset, and movement end for one center-out
trial. Georgopoulos et al. reported a velocity threshold crossing but not its
numeric threshold. This implementation therefore exposes the operational
threshold. `:target_projection` requires velocity toward the target;
`:speed` uses planar speed. The acceleration onset is diagnostic only.
"""
function trial_events(
    session::ReachSession,
    trial::Integer;
    velocity_threshold::Real=0.05,
    min_consecutive::Integer=5,
    velocity_mode::Symbol=:target_projection,
    dims::NTuple{2,Int}=DEFAULT_DIMS,
    acceleration_baseline::Real=0.3,
    acceleration_threshold_sd::Real=3.0,
)
    velocity_threshold > 0 || throw(ArgumentError("velocity_threshold must be positive"))
    trial = _check_trial(session, trial)
    target_index = _target_onset_index(session, trial; dims=dims)
    target_index === nothing && return missing
    target = trial_target(session, trial; dims=dims)
    time = session.times[trial]
    velocity = session.hand_vel[trial]
    count = min(length(time), size(velocity, 1))

    trace = if velocity_mode === :target_projection
        _target_projection(velocity, target.center, target.outer, dims)
    elseif velocity_mode === :speed
        _planar_magnitude(velocity, dims)
    else
        throw(ArgumentError("velocity_mode must be :target_projection or :speed"))
    end
    trace = @view trace[1:count]

    movement_index = _first_sustained(
        trace,
        target_index,
        min_consecutive,
        value -> value >= velocity_threshold,
    )
    movement_index === nothing && return missing

    end_index = _first_sustained(
        trace,
        movement_index + min_consecutive,
        min_consecutive,
        value -> value < velocity_threshold,
    )
    end_index === nothing && (end_index = count)

    target_time = time[target_index]
    movement_time = time[movement_index]
    end_time = time[end_index]
    acceleration = _acceleration_candidate(
        session,
        trial,
        target_time,
        target_index,
        movement_index;
        dims=dims,
        baseline_duration=acceleration_baseline,
        threshold_sd=acceleration_threshold_sd,
    )
    angle = mod(round(Int, target.angle / 45) * 45, 360)

    return (
        trial=trial,
        direction=angle,
        target_onset=target_time,
        acceleration_onset=acceleration.onset,
        movement_onset=movement_time,
        movement_end=end_time,
        reaction_time=movement_time - target_time,
        movement_time=end_time - movement_time,
        total_experimental_time=end_time - target_time,
        velocity_threshold=Float64(velocity_threshold),
        acceleration_threshold=acceleration.threshold,
        velocity_mode=velocity_mode,
    )
end

"""Create one row per trial with target, acceleration, and velocity events."""
function movement_event_table(
    session::ReachSession;
    trials=good_trials(session),
    kwargs...,
)
    rows = NamedTuple[]
    for trial in trials
        event = trial_events(session, trial; kwargs...)
        event === missing || push!(rows, event)
    end
    return DataFrame(rows)
end
