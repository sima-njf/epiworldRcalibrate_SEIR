"""
generate_tuned_predictions.py
==============================
Generates test-set predictions using the hypertuned BiLSTM model (Optuna best trial).
Output: test_bilstm_predictions_tuned.csv  (same schema as test_bilstm_predictions.csv)
        bilstm_timing.csv                  (per-sim mean inference time in seconds)

Usage:
    python generate_tuned_predictions.py

Requires: torch, joblib, numpy, pandas, scikit-learn
"""

import os, sys, time
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
from torch.nn.utils.rnn import pack_padded_sequence
import joblib
import warnings

warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=FutureWarning)

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
BASE_DIR   = os.path.expanduser(
    "~/sima/epiworldRcalibrate_SEIR/SEIR_epi_reconstruction"
)
MODEL_DIR  = os.path.join(BASE_DIR, "model")

MODEL_PATH     = os.path.join(MODEL_DIR, "model_bilstm_tuned.pt")
SCALER_ADD     = os.path.join(MODEL_DIR, "scaler_additional.pkl")
SCALER_TGT     = os.path.join(MODEL_DIR, "scaler_targets.pkl")

TEST_INC_NPY   = os.path.join(BASE_DIR, "test_incidence_raw (2).npy")
TEST_PARAMS    = os.path.join(BASE_DIR, "test_actual_parameters (2).csv")
OUT_CSV        = os.path.join(BASE_DIR, "test_bilstm_predictions_tuned.csv")
TIMING_CSV     = os.path.join(BASE_DIR, "bilstm_timing.csv")

T_MAX   = 365
MIN_LEN = 15
DEVICE  = torch.device("cpu")

# ---------------------------------------------------------------------------
# Model architecture — must match hypertuned BiLSTM exactly
# ---------------------------------------------------------------------------
class BiLSTMRegressor(nn.Module):
    def __init__(self):
        super().__init__()
        self.lstm = nn.LSTM(input_size=1, hidden_size=160, num_layers=1,
                            batch_first=True, bidirectional=True, dropout=0.0)
        self.fc1       = nn.Linear(2 * 160 + 6, 64)
        self.fc2       = nn.Linear(64, 2)
        self.act       = nn.GELU()
        self.head_drop = nn.Dropout(0.18242000617087442)

    def forward(self, x, add_inputs, lengths):
        packed = pack_padded_sequence(x, lengths.cpu(),
                                     batch_first=True, enforce_sorted=False)
        _, (h, _) = self.lstm(packed)
        rep = torch.cat([h[-2], h[-1]], dim=-1)
        z   = self.head_drop(self.act(self.fc1(torch.cat([rep, add_inputs], dim=-1))))
        return torch.sigmoid(self.fc2(z))


def load_model():
    model = BiLSTMRegressor()
    state = torch.load(MODEL_PATH, map_location=DEVICE, weights_only=True)
    model.load_state_dict(state)
    model.to(DEVICE).eval()
    return model

# ---------------------------------------------------------------------------
# RevIN + predict (returns [beta, R0] and elapsed seconds)
# ---------------------------------------------------------------------------
def predict_window(model, scaler_add, scaler_tgt, seq, n, recov, incub):
    x_raw   = np.asarray(seq, dtype=np.float32)
    win_len = len(x_raw)

    wm, ws = float(x_raw.mean()), float(x_raw.std())
    if ws < 1e-6:
        x_norm   = np.zeros_like(x_raw)
        log_mean = float(np.log1p(wm))
        log_std  = 0.0
    else:
        x_norm   = (x_raw - wm) / ws
        log_mean = float(np.log1p(wm))
        log_std  = float(np.log1p(ws))

    add_sc  = scaler_add.transform(np.array([[n, recov, incub]], dtype=np.float32))[0]
    add_row = np.array([
        add_sc[0], add_sc[1], add_sc[2],
        win_len / T_MAX,
        log_mean / 10.0,
        log_std  / 10.0,
    ], dtype=np.float32)

    x_t   = torch.tensor(x_norm,  dtype=torch.float32).unsqueeze(0).unsqueeze(-1)
    add_t = torch.tensor(add_row, dtype=torch.float32).unsqueeze(0)
    len_t = torch.tensor([win_len], dtype=torch.long)

    t0 = time.perf_counter()
    with torch.no_grad():
        pred_sc = model(x_t, add_t, len_t).numpy()
    elapsed = time.perf_counter() - t0

    result = scaler_tgt.inverse_transform(pred_sc)[0]  # [beta, R0]
    return result, elapsed


# ---------------------------------------------------------------------------
# Window definitions (18 windows: 6 lengths x 3 regimes)
# ---------------------------------------------------------------------------
LENGTHS = [15, 30, 60, 90, 180, 365]

def make_windows(t_max=365):
    wins = {}
    for L in LENGTHS:
        if L > t_max:
            continue
        wins[f"early_{L:03d}d"] = {"start": 0,              "len": L, "regime": "early"}
        mid_s = max(0, (t_max - L) // 2)
        wins[f"mid_{L:03d}d"]   = {"start": mid_s,          "len": L, "regime": "mid"}
        late_s = max(0, t_max - L)
        wins[f"late_{L:03d}d"]  = {"start": late_s,         "len": L, "regime": "late"}
    # Special window: first 30 days excluding seed day (days 2-31, 0-indexed start=1)
    wins["noseed_030d"] = {"start": 1, "len": 30, "regime": "noseed"}
    return wins


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    for p in [MODEL_PATH, SCALER_ADD, SCALER_TGT, TEST_INC_NPY, TEST_PARAMS]:
        if not os.path.exists(p):
            sys.exit(f"ERROR: missing file: {p}")

    print("Loading model and scalers...")
    model      = load_model()
    scaler_add = joblib.load(SCALER_ADD)
    scaler_tgt = joblib.load(SCALER_TGT)

    print("Loading test data...")
    inc_raw = np.load(TEST_INC_NPY)         # shape: (n_test, 365)
    params  = pd.read_csv(TEST_PARAMS)

    windows = make_windows()
    n_sims  = len(params)
    print(f"Test sims: {n_sims}  |  Windows: {len(windows)}")

    rows        = []
    timing_rows = []  # per-sim: mean inference time across all windows

    for i, row in params.iterrows():
        sid    = int(row["sim_idx"])
        n_pop  = float(row["n"])
        recov  = float(row["recov"])
        incub  = float(row["incub"])
        beta_t = float(row["beta"])
        R0_t   = float(row["R0"])

        if (i + 1) % 500 == 0:
            print(f"  [{i+1}/{n_sims}]")

        inc_full = inc_raw[i].astype(np.float32)   # 365 days

        sim_times = []

        for win_tag, w in windows.items():
            t0 = w["start"]          # 0-indexed
            t1 = t0 + w["len"]
            seq = inc_full[t0:t1]

            if len(seq) < MIN_LEN:
                continue

            try:
                (beta_pred, R0_pred), elapsed = predict_window(
                    model, scaler_add, scaler_tgt, seq, n_pop, recov, incub
                )
            except Exception as e:
                print(f"  WARN sim {sid} {win_tag}: {e}")
                continue

            sim_times.append(elapsed)

            rows.append({
                "sim_idx":    sid,
                "window":     win_tag,
                "regime":     w["regime"],
                "win_len":    w["len"],
                "t_start":    t0,
                "beta_true":  beta_t,
                "R0_true":    R0_t,
                "beta_pred":  float(beta_pred),
                "R0_pred":    float(R0_pred),
                "time_sec":   elapsed,
            })

        if sim_times:
            timing_rows.append({
                "sim_idx":         sid,
                "mean_time_sec":   float(np.mean(sim_times)),
                "median_time_sec": float(np.median(sim_times)),
                "n_windows":       len(sim_times),
            })

    out = pd.DataFrame(rows)
    out.to_csv(OUT_CSV, index=False)
    print(f"\nSaved {len(out)} rows to:\n  {OUT_CSV}")

    timing_df = pd.DataFrame(timing_rows)
    timing_df.to_csv(TIMING_CSV, index=False)
    mean_t = timing_df["mean_time_sec"].mean()
    print(f"Saved timing to:\n  {TIMING_CSV}")
    print(f"Mean per-window inference time: {mean_t*1000:.2f} ms")

    # Quick sanity check
    for win_tag in ["early_015d", "early_090d", "late_365d"]:
        sub = out[out["window"] == win_tag]
        mae_b = (sub["beta_pred"] - sub["beta_true"]).abs().mean()
        mae_r = (sub["R0_pred"]   - sub["R0_true"]).abs().mean()
        print(f"  {win_tag:15s}  MAE beta={mae_b:.4f}  MAE R0={mae_r:.4f}")


if __name__ == "__main__":
    main()
