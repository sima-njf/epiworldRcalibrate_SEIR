"""
bernardo_predict_measles1861.py
================================
BiLSTM prediction for the EpiEstim Measles1861 smoke-test dataset (48-day
full window only, per instructions -- no 60/90/180/365-day sub-windows).

Reuses the model-loading / predict_window code from bernardo_predict_real.py,
same scale (n_bilstm) from seir_scale_config.json.

Output row is tagged "early_048d" (t_start=0, win_len=48) to match
calibrate_real_5method.R's win_tag = sprintf("early_%03dd", win_len) lookup
convention exactly, so it plugs into process_dataset() unmodified.

Output:
    real_data/bernardo_real_measles1861_predictions.csv
"""
import os, sys, json
import numpy as np
import pandas as pd

REAL_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, REAL_DIR)
from bernardo_predict_real import load_model, predict_window  # noqa: E402

INC_CSV  = os.path.join(REAL_DIR, "measles1861_epiestim.csv")
META_CSV = os.path.join(REAL_DIR, "measles1861_epiestim_meta.csv")
OUT_CSV  = os.path.join(REAL_DIR, "bernardo_real_measles1861_predictions.csv")

RECOV_RATE = 1 / 8      # matches prepare_measles1861_epiestim.R
INCUB_DAYS = 10


def main():
    with open(os.path.join(REAL_DIR, "seir_scale_config.json")) as f:
        scale_cfg = json.load(f)
    n_bilstm = float(scale_cfg["n_bilstm"])

    inc_df = pd.read_csv(INC_CSV)
    meta = pd.read_csv(META_CSV).iloc[0]
    inc = inc_df["daily_cases"].values.astype(np.float32)
    n_pop_true = float(meta["n_pop"])
    n_days = len(inc)
    print(f"Window: {n_days} days, total cases={inc.sum():.0f}, mean={inc.mean():.2f}, "
          f"n_pop_true={n_pop_true:.0f}")

    print("Loading Bernardo BiLSTM model...")
    model = load_model()
    import joblib
    scaler_add = joblib.load(os.path.join(REAL_DIR, "..", "bernardo_model", "scaler_additional (5).pkl"))
    scaler_tgt = joblib.load(os.path.join(REAL_DIR, "..", "bernardo_model", "scaler_targets (5).pkl"))

    # n_pop_true (187) is already below n_bilstm (8000) -- same situation as
    # the existing individual-level Hagelloch measles predictions, so no
    # scaling is applied (scale_factor clamps to 1.0), consistent with how
    # bernardo_predict_real.py's run_dataset() already handles this dataset.
    scale_factor = min(1.0, n_bilstm / n_pop_true) if n_pop_true > 0 else 1.0
    n_pop_model = n_pop_true * scale_factor
    inc_model = inc * scale_factor
    print(f"Scaling for BiLSTM: n_pop {n_pop_true:,.0f} -> {n_pop_model:,.0f} (factor={scale_factor:.6g})")

    beta, R0, elapsed = predict_window(
        model, scaler_add, scaler_tgt, inc_model, n_pop_model, RECOV_RATE, INCUB_DAYS
    )
    print(f"early_{n_days:03d}d  beta={beta:.4f}  R0={R0:.3f}  ({elapsed*1000:.1f} ms)")

    pd.DataFrame([{
        "window": f"early_{n_days:03d}d", "regime": "early", "win_len": n_days, "t_start": 0,
        "beta_pred": beta, "R0_pred": R0, "time_sec": elapsed,
        "n_bilstm_used": n_bilstm,
    }]).to_csv(OUT_CSV, index=False)
    print(f"-> saved {OUT_CSV}")


if __name__ == "__main__":
    main()
