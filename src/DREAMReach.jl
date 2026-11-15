"""
Directional tuning and population-vector decoding for the Flint et al. (2012)
center-out reaching data (DREAM format).

Files, in the order the analysis uses them:

- `session.jl`            load the MAT file; trials, targets, spikes
- `events.jl`             target onset, movement onset and end
- `tuning.jl`             epoch firing rates and the cosine tuning fit
- `perievent.jl`          movement-aligned raster/PSTH data
- `population_vector.jl`  cross-validated population-vector decoding
- `plotting.jl`           figures
"""
module DREAMReach

using DataFrames
using LinearAlgebra
using MAT
using Random
using Statistics
import Plots

include("session.jl")
include("events.jl")
include("tuning.jl")
include("perievent.jl")
include("population_vector.jl")
include("plotting.jl")

# data access
export ReachSession, load_session, ntrials, good_trials, session_summary,
       trial_duration, trial_target, target_direction, direction_label,
       spike_times, candidate_sorted_units, write_csv

# single-unit analysis
export target_onset, trial_events, movement_event_table,
       firing_rate_table, direction_firing_rates, fit_tuning_curve,
       detect_initial_change, peri_event_data

# population vector
export decoding_trials, epoch_rates, fit_all_units, is_tuned, population_vector,
       stratified_folds, cross_validated_pv, percent_within, circular_mean, wrap180,
       windowed_rates, hand_speed_by_direction, time_resolved_pv, time_to_fraction

# figures
export plot_tuning_curve, plot_epoch_tuning, plot_perievent, plot_directional_perievent_grid,
       plot_tuning_summary, plot_population_vectors, plot_single_trial_pv,
       animate_population_vectors, animate_time_resolved_pv, plot_time_resolved_pv

end
