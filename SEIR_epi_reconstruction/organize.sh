#!/bin/bash
# =============================================================================
#  organize.sh  — create folder structure and move data/result files
#
#  Run once from the SEIR_epi_reconstruction directory:
#      bash organize.sh
#
#  What this does:
#    1. Creates results/ and plots/comparison/ and plots/windows/ dirs
#    2. Moves CSV result files (method summaries, daily, etc.) -> results/
#    3. Moves generated PNG plots -> plots/comparison/ or plots/windows/
#    4. Does NOT touch scripts (.R, .py) or input data (test_actual_*.csv,
#       test_incidence_raw.*, bilstm_timing.csv, data_construction/, model/)
#
#  After running, update PROJECT_DIR references only if you move scripts.
#  Currently all scripts keep PROJECT_DIR at the root of SEIR_epi_reconstruction.
# =============================================================================

set -e
cd "$(dirname "$0")"

echo "=== Creating directory structure ==="
mkdir -p results
mkdir -p plots/comparison
mkdir -p plots/windows
mkdir -p archive_old_results

echo "=== Moving method result CSVs -> results/ ==="
for f in \
    abc_seir_summary.csv abc_seir_daily.csv \
    abcsmc_seir_summary.csv abcsmc_seir_daily.csv \
    abcsmc_seir_trace.csv abcsmc_seir_particles.csv \
    nm_seir_summary.csv nm_seir_daily.csv \
    de_seir_summary.csv de_seir_daily.csv \
    comparison_5method_metrics.csv \
    comparison_5method_summary.csv \
    comparison_5method_cost.csv \
    comparison_4method_metrics.csv \
    comparison_4method_metrics_with_oracle.csv \
    comparison_4method_summary.csv \
    comparison_4method_cost.csv \
    comparison_oracle_reference.csv \
    comparison_bayesian_coverage.csv \
    validation_parameter_recovery.csv \
    validation_interval_coverage.csv \
    validation_posterior_geometry.csv \
    validation_ridge_ratio.csv \
    mae_by_r0_window_summary.csv \
    part2_random_pipeline_summary.csv; do
  [ -f "$f" ] && mv "$f" results/ && echo "  moved $f" || true
done

echo "=== Moving comparison PNGs -> plots/comparison/ ==="
for f in \
    comparison_4method_all_metrics.png \
    comparison_4method_accuracy_vs_cost.png \
    comparison_4method_excess_smape.png \
    comparison_4method_incidence_curves.png \
    comparison_4method_mean_bar.png \
    comparison_4method_param_scatter.png \
    comparison_4method_smape.png \
    comparison_5method_all_metrics.png \
    comparison_5method_accuracy_vs_cost.png \
    comparison_5method_incidence_curves.png \
    comparison_5method_mean_bar.png \
    comparison_5method_param_scatter.png \
    comparison_5method_smape.png \
    comparison_5method_timing.png \
    validation_ape_beta_R0.png \
    validation_ape_by_quantity.png \
    validation_ridge_plane.png \
    validation_scatter_all.png \
    mae_by_r0_window.png \
    mae_by_r0_window_combined.png \
    mae_by_r0_window_params.png \
    part1_seir_incidence_curves.png \
    part1_seir_per_day_mae.png \
    part2_seir_random_pipeline.png \
    part3A_seir_param_error_vs_window.png \
    part3B_seir_curve_mae_vs_window.png \
    part3C_seir_example_windows.png \
    part3D_seir_peak_prediction.png; do
  [ -f "$f" ] && mv "$f" plots/comparison/ && echo "  moved $f" || true
done

echo "=== Moving window comparison results -> results/ and plots/windows/ ==="
[ -f window_all_methods_results.csv  ] && mv window_all_methods_results.csv  results/ && echo "  moved window_all_methods_results.csv"  || true
[ -f window_all_methods_summary.csv  ] && mv window_all_methods_summary.csv  results/ && echo "  moved window_all_methods_summary.csv"  || true

echo "=== Archiving old/stale files ==="
for f in \
    Rplots.pdf \
    optimize_seir_submit.R.bak \
    "test_actual_parameters.csv.bak_20260803_174325" \
    "test_incidence_raw.npy.bak_20260803_174325" \
    "test_actual_parameters (2).csv" \
    "test_bilstm_predictions (2).csv" \
    abc_collect.log abcsmc_collect.log opt_collect.log smc_collect.log \
    compare.log validate.log build_test_set.log; do
  [ -f "$f" ] && mv "$f" archive_old_results/ && echo "  archived $f" || true
done

echo ""
echo "=== Directory structure after organizing ==="
ls -1
echo ""
echo "results/:"
ls -1 results/ 2>/dev/null | head -30 || true
echo ""
echo "plots/comparison/:"
ls -1 plots/comparison/ 2>/dev/null | head -30 || true
echo ""
echo "plots/windows/:"
ls -1 plots/windows/ 2>/dev/null | head -30 || true

echo ""
echo "Done. Scripts remain at the root of SEIR_epi_reconstruction/."
echo "Input data files (test_actual_parameters.csv, test_incidence_raw.*,"
echo "  bilstm_timing.csv, data_construction/, model/) are unchanged."
