# Cross-validated population-vector decoding of reach direction
# (Georgopoulos et al., 1986).
#
#   julia --project=. scripts/02_population_vector.jl [mat_path] [epoch] [output_dir]
#
# epoch: tet (target onset -> movement end, default), mt, or rt

using DREAMReach
using DataFrames
using Random
using Statistics
import Plots

const R2_MIN = 0.7      # units must fit a cosine this well to vote
const N_FOLDS = 5
const N_SHUFFLES = 200
const SEED = 20261116

function main(args=ARGS)
    root = normpath(joinpath(@__DIR__, ".."))
    mat_path = length(args) >= 1 ? abspath(args[1]) : joinpath(root, "data", "Flint_2012_e1.mat")
    epoch = length(args) >= 2 ? Symbol(args[2]) : :tet
    out = length(args) >= 3 ? abspath(args[3]) : joinpath(root, "outputs", "population_vector")
    mkpath(out)
    rng = MersenneTwister(SEED)

    # 1. Trials, units, and firing rates
    session = load_session(mat_path)
    units = candidate_sorted_units(session)
    trials = decoding_trials(session)
    rates = epoch_rates(session, trials, units; epoch=epoch)

    # 2. Tuning of every unit on all trials (descriptive only)
    fits = fit_all_units(rates, trials.direction)
    tuning = DataFrame(
        unit=units,
        channel=Int.(session.neuron_mapping[units, 1]),
        sort_code=Int.(session.neuron_mapping[units, 2]),
        mean_rate=vec(mean(rates; dims=1)),
        b0=[f.b0 for f in fits],
        amplitude=[f.amplitude for f in fits],
        preferred_direction=[f.preferred_direction for f in fits],
        r_squared=[f.r_squared for f in fits],
    )
    n_tuned = count(f -> is_tuned(f, R2_MIN), fits)

    # 3. Cross-validated decoding
    folds = stratified_folds(trials.direction, N_FOLDS, rng)
    result = cross_validated_pv(rates, trials, folds; r2_min=R2_MIN)
    single, averaged = result.single, result.averaged

    # 4. Null model: shuffle direction labels in the training trials
    null_median_error = Float64[]
    null_percent_correct = Float64[]
    for _ in 1:N_SHUFFLES
        shuffled = cross_validated_pv(rates, trials, folds; r2_min=R2_MIN, shuffle_labels=rng).single
        nrow(shuffled) == 0 && continue
        push!(null_median_error, median(abs.(shuffled.error_hand)))
        push!(null_percent_correct, percent_within(shuffled.error_target))
    end
    observed = median(abs.(single.error_hand))
    p_value = (1 + count(<=(observed), null_median_error)) / (1 + length(null_median_error))

    # 5. Save tables and figures
    tag = string(epoch)
    write_csv(joinpath(out, "unit_tuning_$tag.csv"), tuning)
    write_csv(joinpath(out, "pv_single_trial_$tag.csv"), single)
    write_csv(joinpath(out, "pv_direction_mean_$tag.csv"), averaged)
    Plots.default(fontfamily="Helvetica", framestyle=:box, grid=false)
    Plots.savefig(plot_tuning_summary(tuning; r2_min=R2_MIN, epoch=epoch),
                  joinpath(out, "tuning_summary_$tag.png"))
    Plots.savefig(plot_population_vectors(averaged, trials), joinpath(out, "pv_direction_mean_$tag.png"))
    Plots.savefig(plot_single_trial_pv(single, null_median_error), joinpath(out, "pv_single_trial_$tag.png"))
    animate_population_vectors(averaged, trials, rates, fits, findall(f -> is_tuned(f, R2_MIN), fits);
                               path=joinpath(out, "pv_direction_mean_$tag.gif"))

    deg(x) = "$(round(x; digits=1))°"
    summary = """
    Population-vector decoding, $(basename(mat_path)), epoch $(uppercase(tag))
    trials: $(nrow(trials))   candidate units: $(length(units))   tuned (R² ≥ $R2_MIN): $n_tuned
    cross-validation: $N_FOLDS folds stratified by direction; units chosen on training folds
                      (mean $(round(mean(single.n_units); digits=1)) units per fold)

    Direction-averaged PV ($(nrow(averaged)) fold × direction means)
      median |error| vs target:          $(deg(median(abs.(averaged.error_target))))
      median |error| vs actual movement: $(deg(median(abs.(averaged.error_hand))))

    Single-trial PV ($(nrow(single)) trials)
      median |error| vs target:          $(deg(median(abs.(single.error_target))))
      median |error| vs actual movement: $(deg(median(abs.(single.error_hand))))
      within ±22.5° of target:           $(round(percent_within(single.error_target); digits=1))% (chance 12.5%)

    Label-shuffle null ($(length(null_median_error)) shuffles)
      median |error| vs movement:        $(deg(median(null_median_error))) (5th percentile $(deg(quantile(null_median_error, 0.05))))
      within ±22.5° of target:           $(round(mean(null_percent_correct); digits=1))%
      p (null error ≤ observed):         $(round(p_value; digits=4))
    """
    write(joinpath(out, "summary_$tag.txt"), summary)
    print(summary)
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()
