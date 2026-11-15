# Population-vector decoding (Georgopoulos et al., 1986) with cross-validation.
#
# Each tuned unit votes along its preferred direction (PD), weighted by how far
# its rate is above or below baseline:
#
#     w_i = (r_i - b0_i) / amplitude_i
#     PV  = sum_i w_i * [cos(PD_i), sin(PD_i)]
#
# Tuning is always fitted on per-direction mean rates of the training trials.

wrap180(angle) = mod(angle + 180, 360) - 180
circular_mean(angles) = mod(atand(sum(sind, angles), sum(cosd, angles)), 360)

"""
    decoding_trials(session; velocity_threshold=0.05, min_consecutive=5)

Event table for good trials, plus `hand_angle`: the direction of the actual
hand displacement from movement onset to movement end.
"""
function decoding_trials(session::ReachSession; velocity_threshold=0.05, min_consecutive=5)
    trials = movement_event_table(
        session;
        velocity_threshold=velocity_threshold,
        min_consecutive=min_consecutive,
        velocity_mode=:target_projection,
    )
    trials.hand_angle = [hand_displacement_angle(session, row) for row in eachrow(trials)]
    return trials
end

function hand_displacement_angle(session::ReachSession, event)
    time = session.times[event.trial]
    position = session.hand_pos[event.trial]
    start = findfirst(>=(event.movement_onset), time)
    stop = findfirst(>=(event.movement_end), time)
    dx = position[stop, 1] - position[start, 1]
    dy = position[stop, 2] - position[start, 2]
    return mod(atand(dy, dx), 360)
end

"""
    epoch_rates(session, trials, units; epoch=:tet)

Firing-rate matrix (trials × units) for one epoch: `:rt`, `:mt`, or `:tet`.
"""
function epoch_rates(session::ReachSession, trials::DataFrame, units; epoch::Symbol=:tet)
    starts, stops = epoch_bounds(trials, epoch)
    rates = zeros(nrow(trials), length(units))
    for (row, trial) in enumerate(trials.trial), (col, unit) in enumerate(units)
        spikes = spike_times(session, trial, unit; relative_to=:absolute)
        rates[row, col] = _spike_count(spikes, starts[row], stops[row]) / (stops[row] - starts[row])
    end
    return rates
end

function epoch_bounds(trials::DataFrame, epoch::Symbol)
    epoch === :rt && return trials.target_onset, trials.movement_onset
    epoch === :mt && return trials.movement_onset, trials.movement_end
    epoch === :tet && return trials.target_onset, trials.movement_end
    throw(ArgumentError("epoch must be :rt, :mt, or :tet"))
end

"""Cosine fit for every unit (column), using per-direction mean rates."""
function fit_all_units(rates::AbstractMatrix, directions::AbstractVector)
    direction_values = sort(unique(directions))
    return map(axes(rates, 2)) do unit
        means = [mean(rates[directions .== d, unit]) for d in direction_values]
        fit_tuning_curve(DataFrame(direction=direction_values, mean_rate=means))
    end
end

is_tuned(fit, r2_min) = isfinite(fit.r_squared) && fit.r_squared >= r2_min && fit.amplitude > 0

"""Direction (degrees) and length of the population vector for one rate vector."""
function population_vector(rates::AbstractVector, fits, units_used)
    x = 0.0
    y = 0.0
    for i in units_used
        fit = fits[i]
        weight = (rates[i] - fit.b0) / fit.amplitude
        x += weight * cosd(fit.preferred_direction)
        y += weight * sind(fit.preferred_direction)
    end
    return (angle=mod(atand(y, x), 360), length=hypot(x, y))
end

"""Assign trials to `k` folds so every direction is spread evenly across folds."""
function stratified_folds(directions::AbstractVector, k::Integer, rng::AbstractRNG)
    folds = zeros(Int, length(directions))
    for d in unique(directions)
        members = shuffle(rng, findall(==(d), directions))
        for (j, trial) in enumerate(members)
            folds[trial] = mod1(j, k)
        end
    end
    return folds
end

"""
    cross_validated_pv(rates, trials, folds; r2_min=0.7, shuffle_labels=nothing)

Fit tuning on the training folds and decode the held-out fold, for every fold.
Returns two tables:

- `single`: one population vector per held-out trial
- `averaged`: one population vector per (fold, direction), computed from
  held-out rates averaged within that direction (the classic PV)

Pass an RNG as `shuffle_labels` to permute training directions (null model).
"""
function cross_validated_pv(rates, trials, folds; r2_min=0.7, shuffle_labels=nothing)
    directions = trials.direction
    single = DataFrame(trial=Int[], fold=Int[], direction=Int[], hand_angle=Float64[],
                       pv_angle=Float64[], pv_length=Float64[], n_units=Int[])
    averaged = DataFrame(fold=Int[], direction=Int[], hand_angle=Float64[],
                         pv_angle=Float64[], pv_length=Float64[], n_units=Int[])

    for fold in sort(unique(folds))
        train = folds .!= fold
        test = findall(folds .== fold)

        train_directions = directions[train]
        if shuffle_labels !== nothing
            train_directions = shuffle(shuffle_labels, train_directions)
        end
        fits = fit_all_units(rates[train, :], train_directions)
        units_used = findall(fit -> is_tuned(fit, r2_min), fits)
        isempty(units_used) && continue

        for row in test
            pv = population_vector(view(rates, row, :), fits, units_used)
            push!(single, (trials.trial[row], fold, directions[row], trials.hand_angle[row],
                           pv.angle, pv.length, length(units_used)))
        end

        for d in sort(unique(directions[test]))
            rows = test[directions[test] .== d]
            pv = population_vector(vec(mean(rates[rows, :]; dims=1)), fits, units_used)
            push!(averaged, (fold, d, circular_mean(trials.hand_angle[rows]),
                             pv.angle, pv.length, length(units_used)))
        end
    end

    for table in (single, averaged)
        table.error_target = wrap180.(table.pv_angle .- table.direction)
        table.error_hand = wrap180.(table.pv_angle .- table.hand_angle)
    end
    return (single=single, averaged=averaged)
end

"""Percentage of decoded directions within ±`limit` degrees (8-way correct for 22.5°)."""
percent_within(errors; limit=22.5) = 100 * mean(abs.(errors) .<= limit)

# ---------------------------------------------------------------------------
# Time-resolved population vector
#
# Rates are counted in a sliding window aligned to movement onset. Tuning is
# still fitted on the training trials' epoch rates (one fit per fold), and the
# population vector is computed from held-out trials at every time step.

"""
    windowed_rates(session, trials, units; times, window=0.1)

Rates in a sliding window of `window` seconds centred on each of `times`
(seconds relative to movement onset). Returns an array (trials × units × times).
"""
function windowed_rates(session::ReachSession, trials::DataFrame, units;
                        times=-0.4:0.02:0.6, window::Real=0.1)
    rates = zeros(nrow(trials), length(units), length(times))
    for (row, event) in enumerate(eachrow(trials)), (col, unit) in enumerate(units)
        spikes = spike_times(session, event.trial, unit; relative_to=:absolute) .- event.movement_onset
        for (k, t) in enumerate(times)
            rates[row, col, k] = _spike_count(spikes, t - window / 2, t + window / 2) / window
        end
    end
    return rates
end

"""Mean planar hand speed per direction at each of `times` (relative to movement onset)."""
function hand_speed_by_direction(session::ReachSession, trials::DataFrame; times=-0.4:0.02:0.6)
    rows = NamedTuple[]
    for d in sort(unique(trials.direction)), (k, t) in enumerate(times)
        speeds = Float64[]
        for event in eachrow(trials[trials.direction .== d, :])
            time = session.times[event.trial]
            i = findfirst(>=(event.movement_onset + t), time)
            i === nothing && continue
            velocity = session.hand_vel[event.trial]
            push!(speeds, hypot(velocity[i, 1], velocity[i, 2]))
        end
        push!(rows, (direction=d, time=t, speed=mean(speeds)))
    end
    return DataFrame(rows)
end

"""
    time_resolved_pv(windowed, epoch_rates, trials, folds; times, r2_min=0.7)

Direction-averaged population vector at every time step, cross-validated.
For each fold, tuning is fitted on the training trials' `epoch_rates`; the
held-out trials' windowed rates are averaged per direction and decoded.
Vectors are averaged over folds. Returns one row per (direction, time).
"""
function time_resolved_pv(windowed::AbstractArray{<:Real,3}, epoch_rates::AbstractMatrix,
                          trials::DataFrame, folds; times=-0.4:0.02:0.6, r2_min=0.7)
    directions = trials.direction
    direction_values = sort(unique(directions))
    x = zeros(length(direction_values), length(times))
    y = zeros(length(direction_values), length(times))
    n_folds = 0

    for fold in sort(unique(folds))
        train = folds .!= fold
        fits = fit_all_units(epoch_rates[train, :], directions[train])
        units_used = findall(fit -> is_tuned(fit, r2_min), fits)
        isempty(units_used) && continue
        n_folds += 1
        for (j, d) in enumerate(direction_values)
            rows = findall((folds .== fold) .& (directions .== d))
            isempty(rows) && continue
            for k in eachindex(times)
                pv = population_vector(vec(mean(windowed[rows, :, k]; dims=1)), fits, units_used)
                x[j, k] += pv.length * cosd(pv.angle)
                y[j, k] += pv.length * sind(pv.angle)
            end
        end
    end

    rows = NamedTuple[]
    for (j, d) in enumerate(direction_values), (k, t) in enumerate(times)
        px = x[j, k] / n_folds
        py = y[j, k] / n_folds
        push!(rows, (direction=d, time=t, x=px, y=py,
                     angle=mod(atand(py, px), 360), length=hypot(px, py)))
    end
    table = DataFrame(rows)
    table.error_target = wrap180.(table.angle .- table.direction)
    return table
end

"""First time at which `values` reaches `fraction` of its maximum."""
function time_to_fraction(times, values; fraction=0.5)
    threshold = minimum(values) + fraction * (maximum(values) - minimum(values))
    return times[findfirst(>=(threshold), values)]
end
