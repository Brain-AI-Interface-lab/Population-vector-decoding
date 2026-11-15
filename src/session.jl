# Loading a DREAM MAT file and accessing trials, targets, and spikes.

const DEFAULT_DIMS = (1, 2)

"""Data needed for trial inspection, spike access, and 2-D reach playback."""
struct ReachSession
    source::String
    times::Vector{Vector{Float64}}
    hand_pos::Vector{Matrix{Float64}}
    hand_vel::Vector{Matrix{Float64}}
    hand_acc::Vector{Matrix{Float64}}
    target_pos::Vector{Matrix{Float64}}
    conditions::Vector{String}
    neurons::Vector{Any}
    neuron_mapping::Matrix{Float64}
end

"""Read a named field from MAT.jl's `MatlabStructArray`."""
function _matfield(struct_array, name::AbstractString)
    index = findfirst(==(name), struct_array.names)
    index === nothing && error("MATLAB struct field not found: $name")
    return struct_array.values[index]
end

function _float_matrices(values)
    return [Matrix{Float64}(value) for value in vec(values)]
end

"""
    load_session(path) -> ReachSession

Load one DREAM Flint MAT file. The planar workspace in this dataset is stored in
columns 1 and 2; column 3 is constant zero.
"""
function load_session(path::AbstractString)
    isfile(path) || error("MAT file not found: $path")
    data = matread(path)
    haskey(data, "Subject") || error("Expected root field `Subject` in $path")

    subject = data["Subject"]
    haskey(subject, "Trial") || error("Expected `Subject.Trial` in $path")
    trials = subject["Trial"]

    times = [Float64.(vec(value)) for value in vec(_matfield(trials, "Time"))]
    hand_pos = _float_matrices(_matfield(trials, "HandPos"))
    hand_vel = _float_matrices(_matfield(trials, "HandVel"))
    hand_acc = _float_matrices(_matfield(trials, "HandAcc"))
    target_pos = _float_matrices(_matfield(trials, "TargetPos"))
    conditions = String.(vec(_matfield(trials, "Condition")))
    neurons = Any[entry for entry in vec(_matfield(trials, "Neuron"))]
    mapping = Matrix{Float64}(subject["Special"]["NeuronMapping"])

    counts = length.((times, hand_pos, hand_vel, hand_acc, target_pos, conditions, neurons))
    all(==(first(counts)), counts) || error("Trial fields have inconsistent lengths: $counts")

    return ReachSession(
        abspath(path),
        times,
        hand_pos,
        hand_vel,
        hand_acc,
        target_pos,
        conditions,
        neurons,
        mapping,
    )
end

ntrials(session::ReachSession) = length(session.times)
good_trials(session::ReachSession) = findall(==("good"), session.conditions)

function _check_trial(session::ReachSession, trial::Integer)
    checkbounds(session.times, trial)
    return Int(trial)
end

function _check_dims(dims::NTuple{2, Int}, matrix::AbstractMatrix)
    all(dim -> dim in axes(matrix, 2), dims) ||
        throw(ArgumentError("dims=$dims are not valid for a matrix with $(size(matrix, 2)) columns"))
    return dims
end

function trial_duration(session::ReachSession, trial::Integer)
    trial = _check_trial(session, trial)
    time = session.times[trial]
    isempty(time) && return 0.0
    return last(time) - first(time)
end

"""
    trial_target(session, trial; dims=(1, 2), tolerance=1e-8)

Return the initial center, farthest outer target, and target angle for one trial.
`outer` and `angle` are `missing` when the trial never presents an outer target.
Rows containing `NaN` are inactive target intervals and are ignored.
"""
function trial_target(
    session::ReachSession,
    trial::Integer;
    dims::NTuple{2, Int}=DEFAULT_DIMS,
    tolerance::Real=1e-8,
)
    trial = _check_trial(session, trial)
    targets = session.target_pos[trial]
    _check_dims(dims, targets)

    rows = [
        row for row in axes(targets, 1)
        if isfinite(targets[row, dims[1]]) && isfinite(targets[row, dims[2]])
    ]
    isempty(rows) && return missing

    first_row = first(rows)
    center = (
        Float64(targets[first_row, dims[1]]),
        Float64(targets[first_row, dims[2]]),
    )
    distances2 = [
        (targets[row, dims[1]] - center[1])^2 +
        (targets[row, dims[2]] - center[2])^2
        for row in rows
    ]
    farthest = argmax(distances2)
    distances2[farthest] <= tolerance^2 &&
        return (center=center, outer=missing, angle=missing)

    outer_row = rows[farthest]
    outer = (
        Float64(targets[outer_row, dims[1]]),
        Float64(targets[outer_row, dims[2]]),
    )
    angle = mod(rad2deg(atan(outer[2] - center[2], outer[1] - center[1])), 360.0)
    return (center=center, outer=outer, angle=angle)
end

function target_direction(
    session::ReachSession,
    trial::Integer;
    dims::NTuple{2, Int}=DEFAULT_DIMS,
)
    target = trial_target(session, trial; dims=dims)
    return target === missing ? missing : target.angle
end

direction_label(::Missing; step::Real=45) = "no outer target"
function direction_label(angle::Real; step::Real=45)
    step > 0 || throw(ArgumentError("step must be positive"))
    rounded = mod(round(Int, angle / step) * step, 360)
    return "$(Int(round(rounded)))°"
end

function _spike_vector(value)
    value isa Number && return [Float64(value)]
    value isa AbstractArray || error("Unsupported spike field type: $(typeof(value))")
    return Float64.(vec(value))
end

"""
    spike_times(session, trial, unit; relative_to=:trial_start)

Return spike timestamps for one trial and one unit item. Use
`relative_to=:absolute` to keep the session timestamps stored in the MAT file.
"""
function spike_times(
    session::ReachSession,
    trial::Integer,
    unit::Integer;
    relative_to::Symbol=:trial_start,
)
    trial = _check_trial(session, trial)
    spike_entries = vec(_matfield(session.neurons[trial], "Spike"))
    checkbounds(spike_entries, unit)
    spikes = _spike_vector(spike_entries[unit])

    if relative_to === :absolute
        return spikes
    elseif relative_to === :trial_start
        return spikes .- first(session.times[trial])
    else
        throw(ArgumentError("relative_to must be :trial_start or :absolute"))
    end
end

"""
Return candidate sorted-unit rows using mapping codes 1-3.

The code meaning follows the conservative DREAM-format interpretation used in
this project; it is not explicitly documented inside the MAT file itself.
"""
function candidate_sorted_units(session::ReachSession)
    return findall(code -> code in 1:3, session.neuron_mapping[:, 2])
end

"""Basic counts and ranges that describe one session."""
function session_summary(session::ReachSession)
    condition_counts = Dict{String, Int}()
    for condition in session.conditions
        condition_counts[condition] = get(condition_counts, condition, 0) + 1
    end

    sample_counts = length.(session.times)
    durations = [trial_duration(session, trial) for trial in 1:ntrials(session)]
    time_steps = [
        median(diff(time)) for time in session.times
        if length(time) > 1
    ]
    direction_counts = Dict{Int, Int}()
    no_outer_targets = 0
    for trial in 1:ntrials(session)
        angle = target_direction(session, trial)
        if angle === missing
            no_outer_targets += 1
        else
            direction = mod(round(Int, angle / 45) * 45, 360)
            direction_counts[direction] = get(direction_counts, direction, 0) + 1
        end
    end

    return (
        source=session.source,
        trials=ntrials(session),
        condition_counts=condition_counts,
        samples_per_trial=extrema(sample_counts),
        duration_seconds=extrema(durations),
        kinematics_hz=1 / median(time_steps),
        unit_items=size(session.neuron_mapping, 1),
        mapping_codes=sort(unique(session.neuron_mapping[:, 2])),
        direction_counts=sort(collect(direction_counts)),
        no_outer_targets=no_outer_targets,
    )
end

"""Write a DataFrame as CSV without adding a CSV package dependency."""
function write_csv(path::AbstractString, table::DataFrame)
    format(value) = value === missing ? "" :
        occursin(r"[,\"\n]", string(value)) ? "\"" * replace(string(value), '"' => "\"\"") * "\"" :
        string(value)
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, join(names(table), ","))
        for row in eachrow(table)
            println(io, join(format.(values(row)), ","))
        end
    end
    return path
end
