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
