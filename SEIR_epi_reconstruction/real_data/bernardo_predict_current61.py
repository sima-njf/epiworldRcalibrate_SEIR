"""
bernardo_predict_current61.py
==============================
BiLSTM prediction for a single window: the CURRENT last-61-day Utah COVID
window (utah_covid_data as packaged today, not the 2020-21 wave1.csv used
elsewhere in this project). Exists because seir_utah_calibration.R fits that
current window, and the pre-computed bernardo_real_covid_predictions.csv
was built from a completely different period (2020-03-18 to 2021-03-17) --
comparing the two would silently mix eras.

Reuses the model-loading / predict_window code from bernardo_predict_real.py
and the same scale (n_bilstm) from seir_scale_config.json, so this is exactly
the same model call, just on the current window instead of the historical one.

Prerequisites:
    Rscript -e 'data("utah_covid_data", package="epiworldRcalibrate"); ...'
    (see calling code in seir_utah_calibration.R / the CSV export step)

Output:
    real_data/bernardo_real_covid_current61_predictions.csv
"""

import os, sys, json
import numpy as np
import pandas as pd
import joblib
import warnings

warnings.filterwarnings("ignore")

REAL_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, REAL_DIR)
from bernardo_predict_real import load_model, predict_window  # noqa: E402

INC_CSV  = os.path.join(REAL_DIR, "utah_covid_current61.csv")
OUT_CSV  = os.path.join(REAL_DIR, "bernardo_real_covid_current61_predictions.csv")

# Same known/assumed values as seir_utah_calibration.R -- keep in sync.
RECOV_RATE = 1 / 7      # INFECTIOUS_PERIOD_DAYS
INCUB_DAYS = 5          # INCUBATION_DAYS


def main():
    with open(os.path.join(REAL_DIR, "seir_scale_config.json")) as f:
        scale_cfg = json.load(f)
    n_bilstm = float(scale_cfg["n_bilstm"])

    inc_df = pd.read_csv(INC_CSV)
    inc = inc_df["daily_cases"].values.astype(np.float32)
    n_days = len(inc)
    print(f"Window: {inc_df['date'].iloc[0]} to {inc_df['date'].iloc[-1]} "
          f"({n_days} days), total cases={inc.sum():.0f}, mean={inc.mean():.1f}")

    print("Loading Bernardo BiLSTM model...")
    model = load_model()
    scaler_add = joblib.load(os.path.join(REAL_DIR, "..", "bernardo_model",
                                           "scaler_additional (5).pkl"))
    scaler_tgt = joblib.load(os.path.join(REAL_DIR, "..", "bernardo_model",
                                           "scaler_targets (5).pkl"))

    # This window's raw counts (mean ~33/day) are already small relative to
    # n_bilstm (8000), unlike the true Utah population case -- but scale
    # anyway for consistency with how every other BiLSTM call in this
    # project is made (population feature must be n_bilstm exactly, and
    # scaling incidence by the same factor keeps the (population, incidence
    # magnitude) relationship the network saw in training).
    n_pop_true = 3.34e6  # UTAH_POPULATION, matches seir_utah_calibration.R
    scale_factor = min(1.0, n_bilstm / n_pop_true)
    inc_scaled = inc * scale_factor

    beta, R0, elapsed = predict_window(
        model, scaler_add, scaler_tgt, inc_scaled, n_bilstm, RECOV_RATE, INCUB_DAYS
    )
    print(f"current_061d  beta={beta:.4f}  R0={R0:.3f}  ({elapsed*1000:.1f} ms)")

    pd.DataFrame([{
        "window": "current_061d", "regime": "current", "win_len": n_days, "t_start": 0,
        "beta_pred": beta, "R0_pred": R0, "time_sec": elapsed,
        "n_bilstm_used": n_bilstm,
    }]).to_csv(OUT_CSV, index=False)
    print(f"-> saved {OUT_CSV}")


if __name__ == "__main__":
    main()
