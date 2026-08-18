"""
bernardo_predict_real.py
========================
Runs the Bernardo BiLSTM on real incidence data (no ground-truth parameters).
Tries every combination of window length × regime that has at least MIN_LEN
informative days.

Prerequisites (run first):
    Rscript real_data/prepare_utah_covid.R
    Rscript real_data/prepare_measles.R

Output:
    real_data/bernardo_real_covid_predictions.csv
    real_data/bernardo_real_measles_predictions.csv

Usage:
    cd SEIR_epi_reconstruction
    python real_data/bernardo_predict_real.py
"""

import os, sys, time, json
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
from torch.nn.utils.rnn import pack_padded_sequence
import joblib
import warnings

warnings.filterwarnings("ignore")

BASE_DIR  = os.path.expanduser(
    "~/sima/epiworldRcalibrate_SEIR/SEIR_epi_reconstruction"
)
MODEL_DIR = os.path.join(BASE_DIR, "bernardo_model")
REAL_DIR  = os.path.join(BASE_DIR, "real_data")

MODEL_PATH = os.path.join(MODEL_DIR, "model_bilstm_tuned.pt")
SCALER_ADD = os.path.join(MODEL_DIR, "scaler_additional (5).pkl")
SCALER_TGT = os.path.join(MODEL_DIR, "scaler_targets (5).pkl")
MODEL_CFG  = os.path.join(MODEL_DIR, "model_config.json")

T_MAX   = 365
MIN_LEN = 15
MIN_SUM = 5       # skip windows with fewer total cases (uninformative)
DEVICE  = torch.device("cpu")

# The BiLSTM was trained on simulated populations with n in [5000, 10000]
# (see test_actual_parameters (2).csv). Real population sizes (e.g. Utah's
# 3.27M) are ~300-650x outside that range: scaler_additional.transform()
# maps them to huge out-of-[0,1] values, which saturates the sigmoid output
# head to a fixed corner (beta_pred/R0_pred pinned at the scaler's max,
# identical for every window). To keep the model in-distribution we scale
# the population (and incidence, proportionally, since the "additional"
# features include log(mean)/log(std) of the window in absolute case
# counts) down to a population inside the trained range before predicting.
#
# This scale is read from seir_scale_config.json -- the single shared source
# of truth also read by calibrate_real_5method.R (n_classical there). Do NOT
# hard-code a value here: two scripts each inventing their own population
# assumption with no cross-check is exactly the bug pattern found 2026-08-17
# on the SIR side of this project. The scale actually used is stamped into
# every output row (n_bilstm_used column) so a stale/edited config is
# detectable downstream instead of silently producing wrong predictions.
with open(os.path.join(REAL_DIR, "seir_scale_config.json")) as _f:
    _scale_cfg = json.load(_f)
N_SCALED = float(_scale_cfg["n_bilstm"])


# ── Model (identical architecture to generate_bernardo_predictions.py) ─────────

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
        rep = torch.cat([h[-2], h[-1]], dim=-1)
        z   = self.head_drop(self.act(
                  self.fc1(torch.cat([rep, add_inputs], dim=-1))))
        return torch.sigmoid(self.fc2(z))


def load_model():
    with open(MODEL_CFG) as f:
        cfg = json.load(f)["architecture"]
    state  = torch.load(MODEL_PATH, map_location=DEVICE, weights_only=True)
    fc_dim = state["fc1.weight"].shape[0]
    model  = BernardoBiLSTM(
        hidden         = cfg["hidden"],
        num_layers     = cfg["num_layers"],
        lstm_dropout   = cfg["lstm_dropout"],
        head_dropout   = cfg["head_dropout"],
        additional_dim = cfg["additional_dim"],
        fc_dim         = fc_dim,
    )
    model.load_state_dict(state)
    model.to(DEVICE).eval()
    print(f"  hidden={cfg['hidden']}, layers={cfg['num_layers']}, "
          f"additional_dim={cfg['additional_dim']}, fc_dim={fc_dim}")
    return model


def shape_features(x_raw):
    x_raw = np.asarray(x_raw, dtype=np.float64); n = len(x_raw)
    lc    = np.log1p(x_raw); t = np.arange(n, dtype=np.float64)
    slope = float(np.polyfit(t, lc, 1)[0]) if n > 1 else 0.0
    half  = max(1, n // 2)
    if half > 1 and n - half > 1:
        s1 = float(np.polyfit(t[:half], lc[:half], 1)[0])
        s2 = float(np.polyfit(t[half:], lc[half:], 1)[0])
        curvature = s2 - s1
    else:
        curvature = 0.0
    argmax_frac     = float(np.argmax(x_raw)) / max(1, n - 1)
    q               = max(1, n // 4)
    late_over_early = (float(x_raw[-q:].mean()) + 1.0) / (float(x_raw[:q].mean()) + 1.0)
    return [float(np.tanh(slope)), float(np.tanh(curvature)),
            argmax_frac, float(np.tanh(np.log1p(late_over_early) / 3))]


def predict_window(model, scaler_add, scaler_tgt, seq, n_pop, recov, incub):
    x_raw   = np.asarray(seq, dtype=np.float32)
    win_len = len(x_raw)
    wm, ws  = float(x_raw.mean()), float(x_raw.std())
    if ws < 1e-6:
        x_norm   = np.zeros_like(x_raw)
        log_mean, log_std = float(np.log1p(wm)), 0.0
    else:
        x_norm   = (x_raw - wm) / ws
        log_mean = float(np.log1p(wm))
        log_std  = float(np.log1p(ws))

    add_sc  = scaler_add.transform(
        np.array([[n_pop, recov, incub]], dtype=np.float32))[0]
    add_row = np.array([
        add_sc[0], add_sc[1], add_sc[2],
        win_len / T_MAX,
        log_mean / 10.0,
        log_std  / 10.0,
    ] + shape_features(x_raw), dtype=np.float32)

    x_t   = torch.tensor(x_norm, dtype=torch.float32).unsqueeze(0).unsqueeze(-1)
    add_t = torch.tensor(add_row, dtype=torch.float32).unsqueeze(0)
    len_t = torch.tensor([win_len], dtype=torch.long)

    t0 = time.perf_counter()
    with torch.no_grad():
        pred_sc = model(x_t, add_t, len_t).numpy()
    elapsed = time.perf_counter() - t0

    result = scaler_tgt.inverse_transform(pred_sc)[0]
    return float(result[0]), float(result[1]), elapsed   # beta, R0, sec


def run_dataset(model, scaler_add, scaler_tgt,
                inc, n_pop, recov, incub, label, out_csv):
    """Run all windows on a single incidence array and save predictions."""
    n_days  = len(inc)
    lengths = [L for L in [15, 30, 60, 90, 180, 365] if L <= n_days]

    # Scale population + incidence into the model's trained range (see
    # N_SCALED comment above). beta = crate * ptran is an intensive
    # (population-independent) rate, so predictions need no re-scaling;
    # only the raw case counts and the population feature are scaled down.
    scale_factor = min(1.0, N_SCALED / n_pop) if n_pop > 0 else 1.0
    n_pop_model  = n_pop * scale_factor
    inc_model    = inc * scale_factor
    print(f"  Scaling for BiLSTM: n_pop {n_pop:,.0f} -> {n_pop_model:,.0f} "
          f"(factor={scale_factor:.6g})")

    rows = []
    for L in lengths:
        for regime, t0 in [
            ("early", 0),
            ("mid",   max(0, (n_days - L) // 2)),
            ("late",  max(0, n_days - L)),
        ]:
            seq = inc_model[t0 : t0 + L]
            if len(seq) < MIN_LEN or float(np.sum(seq)) < MIN_SUM:
                continue
            try:
                beta, R0, elapsed = predict_window(
                    model, scaler_add, scaler_tgt,
                    seq, n_pop_model, recov, incub
                )
                rows.append({
                    "window":    f"{regime}_{L:03d}d",
                    "regime":    regime,
                    "win_len":   L,
                    "t_start":   t0,
                    "beta_pred": beta,
                    "R0_pred":   R0,
                    "time_sec":  elapsed,
                    "n_bilstm_used": N_SCALED,
                })
                print(f"  {regime}_{L:03d}d  beta={beta:.4f}  R0={R0:.3f}"
                      f"  ({elapsed*1000:.1f} ms)")
            except Exception as e:
                print(f"  WARN {regime}_{L:03d}d: {e}")

    # Full-outbreak window when shorter than 365 d
    if n_days < 365 and n_days >= MIN_LEN and float(np.sum(inc)) >= MIN_SUM:
        try:
            beta, R0, elapsed = predict_window(
                model, scaler_add, scaler_tgt,
                inc_model, n_pop_model, recov, incub
            )
            rows.append({
                "window":    f"full_{n_days:03d}d",
                "regime":    "full",
                "win_len":   n_days,
                "t_start":   0,
                "beta_pred": beta,
                "R0_pred":   R0,
                "time_sec":  elapsed,
            })
            print(f"  full_{n_days:03d}d  beta={beta:.4f}  R0={R0:.3f}")
        except Exception as e:
            print(f"  WARN full window: {e}")

    df = pd.DataFrame(rows)
    df.to_csv(out_csv, index=False)
    print(f"\n  -> {len(df)} predictions saved to:\n     {out_csv}\n")
    return df


def main():
    os.makedirs(REAL_DIR, exist_ok=True)

    for p in [MODEL_PATH, SCALER_ADD, SCALER_TGT, MODEL_CFG]:
        if not os.path.exists(p):
            sys.exit(f"ERROR: missing file: {p}")

    print("Loading Bernardo BiLSTM model...")
    model      = load_model()
    scaler_add = joblib.load(SCALER_ADD)
    scaler_tgt = joblib.load(SCALER_TGT)

    # ── Utah COVID ─────────────────────────────────────────────────────────────
    covid_inc  = os.path.join(REAL_DIR, "utah_covid_wave1.csv")
    covid_meta = os.path.join(REAL_DIR, "utah_covid_meta.csv")

    if os.path.exists(covid_inc) and os.path.exists(covid_meta):
        print("\n=== Utah COVID-19 (wave 1, 365 days) ===")
        inc_df   = pd.read_csv(covid_inc)
        meta_df  = pd.read_csv(covid_meta).iloc[0]
        inc      = inc_df["daily_cases"].values.astype(np.float32)
        run_dataset(
            model, scaler_add, scaler_tgt,
            inc,
            n_pop  = float(meta_df["n_pop"]),
            recov  = float(meta_df["recov_rate"]),
            incub  = float(meta_df["incub_days"]),
            label  = "covid",
            out_csv = os.path.join(REAL_DIR, "bernardo_real_covid_predictions.csv"),
        )
    else:
        print(f"\nSkipping COVID — run prepare_utah_covid.R first.")

    # ── Hagelloch Measles ──────────────────────────────────────────────────────
    meas_inc  = os.path.join(REAL_DIR, "measles_hagelloch_incidence.csv")
    meas_meta = os.path.join(REAL_DIR, "measles_hagelloch_meta.csv")

    if os.path.exists(meas_inc) and os.path.exists(meas_meta):
        print("\n=== Hagelloch 1861 Measles ===")
        inc_df   = pd.read_csv(meas_inc)
        meta_df  = pd.read_csv(meas_meta).iloc[0]
        inc      = inc_df["daily_cases"].values.astype(np.float32)
        run_dataset(
            model, scaler_add, scaler_tgt,
            inc,
            n_pop  = float(meta_df["n_pop"]),
            recov  = float(meta_df["recov_rate"]),
            incub  = float(meta_df["incub_days"]),
            label  = "measles",
            out_csv = os.path.join(REAL_DIR, "bernardo_real_measles_predictions.csv"),
        )
    else:
        print(f"\nSkipping measles — run prepare_measles.R first.")


if __name__ == "__main__":
    main()
