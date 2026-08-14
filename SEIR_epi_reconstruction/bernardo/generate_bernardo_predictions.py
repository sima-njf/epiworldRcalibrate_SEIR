"""
generate_bernardo_predictions.py
=================================
Generates test-set predictions using the Bernardo BiLSTM model.

Architecture (from model_config.json):
  hidden=64, num_layers=3, bidirectional, additional_dim=10
  10 additional features:
    [n, recov, incub]  <- scaled by scaler_additional (5).pkl
    [win_len/T_MAX, log_mean/10, log_std/10,
     shape_logslope, shape_curvature, shape_argmax_frac, shape_late_over_early]

Output:
  bernardo_bilstm_predictions.csv
  bernardo_bilstm_timing.csv

Usage:
    python generate_bernardo_predictions.py
"""

import os, sys, time, json
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
BASE_DIR    = os.path.expanduser(
    "~/sima/epiworldRcalibrate_SEIR/SEIR_epi_reconstruction"
)
MODEL_DIR   = os.path.join(BASE_DIR, "bernardo_model")
BERNARDO_DIR = os.path.join(BASE_DIR, "bernardo")

MODEL_PATH  = os.path.join(MODEL_DIR, "model_bilstm_tuned.pt")
SCALER_ADD  = os.path.join(MODEL_DIR, "scaler_additional (5).pkl")
SCALER_TGT  = os.path.join(MODEL_DIR, "scaler_targets (5).pkl")
MODEL_CFG   = os.path.join(MODEL_DIR, "model_config.json")

TEST_INC_NPY = os.path.join(BASE_DIR, "test_incidence_raw (2).npy")
TEST_PARAMS  = os.path.join(BASE_DIR, "test_actual_parameters (2).csv")
OUT_CSV      = os.path.join(BERNARDO_DIR, "bernardo_bilstm_predictions.csv")
TIMING_CSV   = os.path.join(BERNARDO_DIR, "bernardo_bilstm_timing.csv")

T_MAX   = 365
MIN_LEN = 15
DEVICE  = torch.device("cpu")

# ---------------------------------------------------------------------------
# Model architecture (matches model_config.json)
# ---------------------------------------------------------------------------
class BernardoBiLSTM(nn.Module):
    def __init__(self, hidden, num_layers, lstm_dropout, head_dropout,
                 additional_dim, fc_dim, output_dim=2):
        super().__init__()
        drop = lstm_dropout if num_layers > 1 else 0.0
        self.lstm = nn.LSTM(input_size=1, hidden_size=hidden,
                            num_layers=num_layers, batch_first=True,
                            bidirectional=True, dropout=drop)
        self.fc1       = nn.Linear(2 * hidden + additional_dim, fc_dim)
        self.fc2       = nn.Linear(fc_dim, output_dim)
        self.act       = nn.GELU()
        self.head_drop = nn.Dropout(head_dropout)

    def forward(self, x, add_inputs, lengths):
        packed = pack_padded_sequence(x, lengths.cpu(),
                                      batch_first=True, enforce_sorted=False)
        _, (h, _) = self.lstm(packed)
        rep = torch.cat([h[-2], h[-1]], dim=-1)   # last fwd + last bwd
        z   = self.head_drop(self.act(
                  self.fc1(torch.cat([rep, add_inputs], dim=-1))))
        return torch.sigmoid(self.fc2(z))


def load_model():
    with open(MODEL_CFG) as f:
        cfg = json.load(f)["architecture"]

    state = torch.load(MODEL_PATH, map_location=DEVICE, weights_only=True)

    # Infer fc_dim from saved weights
    fc_dim = state["fc1.weight"].shape[0]

    model = BernardoBiLSTM(
        hidden         = cfg["hidden"],
        num_layers     = cfg["num_layers"],
        lstm_dropout   = cfg["lstm_dropout"],
        head_dropout   = cfg["head_dropout"],
        additional_dim = cfg["additional_dim"],
        fc_dim         = fc_dim,
        output_dim     = cfg["output_dim"],
    )
    model.load_state_dict(state)
    model.to(DEVICE).eval()
    print(f"  hidden={cfg['hidden']}, layers={cfg['num_layers']}, "
          f"additional_dim={cfg['additional_dim']}, fc_dim={fc_dim}")
    return model


# ---------------------------------------------------------------------------
# Shape features — matches DeepIMC/generate_tuned_predictions.py exactly
# ---------------------------------------------------------------------------
def shape_features(x_raw):
    x_raw = np.asarray(x_raw, dtype=np.float64)
    n     = len(x_raw)
    lc    = np.log1p(x_raw)
    t     = np.arange(n, dtype=np.float64)

    slope = float(np.polyfit(t, lc, 1)[0]) if n > 1 else 0.0

    half = max(1, n // 2)
    if half > 1 and n - half > 1:
        s1        = float(np.polyfit(t[:half], lc[:half], 1)[0])
        s2        = float(np.polyfit(t[half:], lc[half:], 1)[0])
        curvature = s2 - s1
    else:
        curvature = 0.0

    argmax_frac = float(np.argmax(x_raw)) / max(1, n - 1)

    q               = max(1, n // 4)
    first_q         = float(x_raw[:q].mean())
    last_q          = float(x_raw[-q:].mean())
    late_over_early = (last_q + 1.0) / (first_q + 1.0)

    return [float(np.tanh(slope)),
            float(np.tanh(curvature)),
            argmax_frac,
            float(np.tanh(np.log1p(late_over_early) / 3))]


# ---------------------------------------------------------------------------
# Predict one window
# ---------------------------------------------------------------------------
def predict_window(model, scaler_add, scaler_tgt, seq, n_pop, recov, incub):
    x_raw   = np.asarray(seq, dtype=np.float32)
    win_len = len(x_raw)

    # RevIN normalisation
    wm, ws = float(x_raw.mean()), float(x_raw.std())
    if ws < 1e-6:
        x_norm   = np.zeros_like(x_raw)
        log_mean, log_std = float(np.log1p(wm)), 0.0
    else:
        x_norm   = (x_raw - wm) / ws
        log_mean = float(np.log1p(wm))
        log_std  = float(np.log1p(ws))

    # Scaled population features (n, recov, incub)
    add_sc = scaler_add.transform(
        np.array([[n_pop, recov, incub]], dtype=np.float32))[0]

    # Full 10-dim additional vector (shape features match DeepIMC exactly)
    add_row = np.array([
        add_sc[0], add_sc[1], add_sc[2],   # scaled n, recov, incub
        win_len / T_MAX,                     # win_len/T_MAX
        log_mean / 10.0,                     # log_mean/10
        log_std  / 10.0,                     # log_std/10
    ] + shape_features(x_raw), dtype=np.float32)

    x_t   = torch.tensor(x_norm, dtype=torch.float32).unsqueeze(0).unsqueeze(-1)
    add_t = torch.tensor(add_row, dtype=torch.float32).unsqueeze(0)
    len_t = torch.tensor([win_len], dtype=torch.long)

    t0 = time.perf_counter()
    with torch.no_grad():
        pred_sc = model(x_t, add_t, len_t).numpy()
    elapsed = time.perf_counter() - t0

    result = scaler_tgt.inverse_transform(pred_sc)[0]   # [beta, R0]
    return result, elapsed


# ---------------------------------------------------------------------------
# Window definitions (18 windows: 6 lengths x 3 regimes)
# ---------------------------------------------------------------------------
LENGTHS = [15, 30, 60, 90, 180, 365]

def make_windows(t_max=365):
    wins = {}
    # Standard windows: 6 lengths x 3 regimes
    for L in LENGTHS:
        if L > t_max:
            continue
        wins[f"early_{L:03d}d"] = {"start": 0,               "len": L, "regime": "early"}
        mid_s = max(0, (t_max - L) // 2)
        wins[f"mid_{L:03d}d"]   = {"start": mid_s,           "len": L, "regime": "mid"}
        late_s = max(0, t_max - L)
        wins[f"late_{L:03d}d"]  = {"start": late_s,          "len": L, "regime": "late"}
    # Special window: first 30 days excluding seed day (days 2-31, 0-indexed start=1)
    wins["noseed_030d"] = {"start": 1, "len": 30, "regime": "noseed"}
    return wins


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    for p in [MODEL_PATH, SCALER_ADD, SCALER_TGT, MODEL_CFG,
              TEST_INC_NPY, TEST_PARAMS]:
        if not os.path.exists(p):
            sys.exit(f"ERROR: missing file: {p}")

    os.makedirs(BERNARDO_DIR, exist_ok=True)

    print("Loading Bernardo model...")
    model      = load_model()
    scaler_add = joblib.load(SCALER_ADD)
    scaler_tgt = joblib.load(SCALER_TGT)

    print("Loading test data...")
    inc_raw = np.load(TEST_INC_NPY)
    params  = pd.read_csv(TEST_PARAMS)

    windows = make_windows()
    n_sims  = len(params)
    print(f"Test sims: {n_sims}  |  Windows: {len(windows)}")

    rows        = []
    timing_rows = []

    for i, row in params.iterrows():
        sid    = int(row["sim_idx"])
        n_pop  = float(row["n"])
        recov  = float(row["recov"])
        incub  = float(row["incub"])
        beta_t = float(row["beta"])
        R0_t   = float(row["R0"])

        if (i + 1) % 500 == 0:
            print(f"  [{i+1}/{n_sims}]")

        inc_full  = inc_raw[i].astype(np.float32)
        sim_times = []

        for win_tag, w in windows.items():
            t0  = w["start"]
            t1  = t0 + w["len"]
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
                "sim_idx":   sid,
                "window":    win_tag,
                "regime":    w["regime"],
                "win_len":   w["len"],
                "t_start":   t0,
                "beta_true": beta_t,
                "R0_true":   R0_t,
                "beta_pred": float(beta_pred),
                "R0_pred":   float(R0_pred),
                "time_sec":  elapsed,
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

    # Sanity check
    for win_tag in ["early_015d", "early_090d", "late_365d"]:
        sub = out[out["window"] == win_tag]
        if len(sub) == 0:
            continue
        mae_b = (sub["beta_pred"] - sub["beta_true"]).abs().mean()
        mae_r = (sub["R0_pred"]   - sub["R0_true"]).abs().mean()
        print(f"  {win_tag:15s}  MAE beta={mae_b:.4f}  MAE R0={mae_r:.4f}")


if __name__ == "__main__":
    main()
