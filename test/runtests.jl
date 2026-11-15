using Test
using DREAMReach
using DataFrames
using Statistics
using Random

mat_path = normpath(joinpath(@__DIR__, "..", "data", "Flint_2012_e1.mat"))
session = load_session(mat_path)
summary = session_summary(session)

@testset "Flint session structure" begin
    @test ntrials(session) == 194
    @test length(good_trials(session)) == 175
    @test summary.condition_counts == Dict("good" => 175, "bad" => 19)
    @test summary.samples_per_trial == (220, 514)
    @test summary.unit_items == 196
    @test isapprox(summary.kinematics_hz, 100.0; atol=0.1)
    @test !isempty(candidate_sorted_units(session))
end

@testset "Target direction" begin
    @test isapprox(target_direction(session, 1), 225.0; atol=0.2)
    @test target_direction(session, 3) === missing
end

@testset "Spike access" begin
    unit = first(candidate_sorted_units(session))
    relative = spike_times(session, 1, unit)
    absolute = spike_times(session, 1, unit; relative_to=:absolute)
    @test length(relative) == length(absolute)
    @test relative ≈ absolute .- first(session.times[1])
end

@testset "Movement events" begin
    event = trial_events(session, 1)
    @test event !== missing
    @test event.target_onset < event.movement_onset < event.movement_end
    @test event.direction == 225
    @test event.velocity_threshold == 0.05

    events = movement_event_table(session)
    @test nrow(events) == 175
    @test all(events.reaction_time .> 0)
    @test all(events.movement_time .> 0)
end

@testset "Directional tuning" begin
    unit = first(candidate_sorted_units(session))
    rates = firing_rate_table(session, unit)
    summary = direction_firing_rates(rates; epoch=:tet)
    fit = fit_tuning_curve(summary)

    @test nrow(rates) == 175
    @test summary.direction == collect(0:45:315)
    @test all(summary.n_trials .>= 14)
    @test isapprox(fit.b0, mean(summary.mean_rate); atol=1e-10)
    @test 0 <= fit.preferred_direction < 360
    @test fit.r_squared > 0.7
end

@testset "Georgopoulos perievent rule" begin
    rates = [10.0, 11.0, 16.0, 17.0, 12.0]
    change = detect_initial_change(rates, 10.0)
    @test change.index == 2
    @test change.kind == :increase

    unit = first(candidate_sorted_units(session))
    peri = peri_event_data(session, unit; direction=45)
    @test length(peri.trials) == 6
    @test isapprox(peri.bin_width, 0.02)
    @test length(peri.counts) == length(peri.centers)
end

@testset "Population vector" begin
    @test wrap180(350.0) == -10.0
    @test isapprox(circular_mean([350.0, 10.0]), 0.0; atol=1e-9) ||
          isapprox(circular_mean([350.0, 10.0]), 360.0; atol=1e-9)

    # Two ideal cosine units with PDs at 0° and 90° decode 45° exactly.
    fits = [(b0=10.0, amplitude=5.0, preferred_direction=0.0),
            (b0=10.0, amplitude=5.0, preferred_direction=90.0)]
    rates = [10.0 + 5cosd(45), 10.0 + 5cosd(45)]
    @test isapprox(population_vector(rates, fits, 1:2).angle, 45.0; atol=1e-9)

    folds = stratified_folds(repeat(0:45:315, 5), 5, MersenneTwister(1))
    @test all(count(==(k), folds) == 8 for k in 1:5)
end
