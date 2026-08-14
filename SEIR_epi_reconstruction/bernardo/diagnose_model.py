"""
diagnose_model.py
=================
Sanity-checks the Bernardo BiLSTM model and scalers:
  1. Prints scaler dimensions and parameter ranges
  2. Runs predictions on 5 test cases from the test set
  3. Reports whether predictions are sensible

Usage:
    cd SEIR_epi_reconstruction
    python bernardo/diagnose_model.py
"""

import os, sys, json
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
from torch.nn.utils.rnn import pack_padded_sequence
import joblib

BASE_DIR   = os.path.expanduser(
    "~/sima/epiworldRcalibrate_SEIR/SEIR_epi_reconstruction"
)
MODEL_DIR  = os.path.join(BASE_DIR, "bernardo_model")
MODEL_PATH = os.path.join(MODEL_DIR, "model_bilstm_tuned.pt")
SCALER_ADD = os.path.join(MODEL_DIR, "scaler_additional (5).pkl")
SCALER_TGT = os.path.join(MODEL_DIR, "scaler_targets (5).pkl")
MODEL_CFG  = os.path.join(MODEL_DIR, "model_config.json")
TEST_INC   = os.path.join(BASE_DIR,  "test_incidence_raw (2).npy")
TEST_PARAMS = os.path.join(BASE_DIR, "test_actual_parameters (2).csv")

# ── 1. Scaler inspection ──────────────────────────────────────────────────────
print("=" * 60)
print("SCALER INSPECTION")
print("=" * 60)

sc_add = joblib.load(SCALER_ADD)
sc_tgt = joblib.load(SCALER_TGT)

print(f"\nscaler_additional (5).pkl:")
print(f"  type             : {type(sc_add).__name__}")
print(f"  n_features_in_   : {sc_add.n_features_in_}")
print(f"  feature_names_in_: {getattr(sc_add, 'feature_names_in_', 'not set')}")
print(f"  data_min_        : {sc_add.data_min_}")
print(f"  data_max_        : {sc_add.data_max_}")
print(f"  scale_           : {sc_add.scale_}")

print(f"\nscaler_targets (5).pkl:")
print(f"  type             : {type(sc_tgt).__name__}")
print(f"  n_features_in_   : {sc_tgt.n_features_in_}")
print(f"  feature_names_in_: {getattr(sc_tgt, 'feature_names_in_', 'not set')}")
print(f"  data_min_        : {sc_tgt.data_min_}")
print(f"  data_max_        : {sc_tgt.data_max_}")
print(f"  scale_           : {sc_tgt.scale_}")

# Check what range of targets the scaler covers
with open(MODEL_CFG) as f:
    cfg = json.load(f)

print(f"\nmodel_config.json:")
print(f"  hidden={cfg['architecture']['hidden']}, "
      f"num_layers={cfg['architecture']['num_layers']}, "
      f"additional_dim={cfg['architecture']['additional_dim']}, "
      f"output_dim={cfg['architecture']['output_dim']}")

# ── 2. Load model and check state dict ───────────────────────────────────────
print("\n" + "=" * 60)
print("MODEL INSPECTION")
print("=" * 60)

state = torch.load(MODEL_PATH, map_location="cpu", weights_only=True)
print("\nState dict keys and shapes:")
for k, v in state.items():
    print(f"  {k:40s}  {tuple(v.shape)}")

# Infer fc_dim from weights
fc_dim      = state["fc1.weight"].shape[0]
add_dim_actual = state["fc1.weight"].shape[1] - 2 * cfg["architecture"]["hidden"]
print(f"\n  Inferred fc_dim         : {fc_dim}")
print(f"  Inferred additional_dim : {add_dim_actual}  "
      f"(config says {cfg['architecture']['additional_dim']})")
if add_dim_actual != cfg["architecture"]["additional_dim"]:
    print("  *** MISMATCH — additional_dim in weights != model_config.json ***")

# ── 3. Rebuild model ─────────────────────────────────────────────────────────
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

arch = cfg["architecture"]
model = BernardoBiLSTM(
    hidden         = arch["hidden"],
    num_layers     = arch["num_layers"],
    lstm_dropout   = arch["lstm_dropout"],
    head_dropout   = arch["head_dropout"],
    additional_dim = add_dim_actual,   # use actual dim from weights
    fc_dim         = fc_dim,
)
result = model.load_state_dict(state, strict=True)
model.eval()
print(f"\n  load_state_dict: missing={result.missing_keys}, "
      f"unexpected={result.unexpected_keys}")
print("  (both should be empty [])")

# ── 4. Predictions on test cases ─────────────────────────────────────────────
print("\n" + "=" * 60)
print("TEST PREDICTIONS vs GROUND TRUTH")
print("=" * 60)

if not os.path.exists(TEST_INC) or not os.path.exists(TEST_PARAMS):
    print("Test data not found — skipping prediction check.")
    sys.exit(0)

inc_raw = np.load(TEST_INC)
params  = pd.read_csv(TEST_PARAMS)

T_MAX = 365
np.random.seed(0)
sample_idx = np.random.choice(len(params), size=10, replace=False)

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

print(f"\n{'Sim':>6}  {'beta_true':>9}  {'R0_true':>7}  "
      f"{'beta_pred':>9}  {'R0_pred':>7}  {'beta_err%':>9}  {'R0_err%':>7}")
print("-" * 70)

errors_beta, errors_R0 = [], []

for i in sample_idx:
    row   = params.iloc[i]
    n_pop = float(row["n"])
    recov = float(row["recov"])
    incub = float(row["incub"])
    b_true = float(row["beta"])
    R0_true = float(row["R0"])

    # Full 365-day window
    x_raw   = inc_raw[i].astype(np.float32)
    win_len = len(x_raw)
    wm, ws  = float(x_raw.mean()), float(x_raw.std())
    if ws < 1e-6:
        x_norm   = np.zeros_like(x_raw); log_mean, log_std = float(np.log1p(wm)), 0.0
    else:
        x_norm   = (x_raw - wm) / ws
        log_mean = float(np.log1p(wm)); log_std = float(np.log1p(ws))

    add_sc  = sc_add.transform(np.array([[n_pop, recov, incub]], dtype=np.float32))[0]
    base    = [add_sc[0], add_sc[1], add_sc[2],
               win_len / T_MAX, log_mean / 10.0, log_std / 10.0]
    if add_dim_actual == 10:
        add_row = np.array(base + shape_features(x_raw), dtype=np.float32)
    elif add_dim_actual == 6:
        add_row = np.array(base, dtype=np.float32)
    else:
        print(f"  Unknown additional_dim={add_dim_actual}!")
        sys.exit(1)

    x_t   = torch.tensor(x_norm, dtype=torch.float32).unsqueeze(0).unsqueeze(-1)
    add_t = torch.tensor(add_row, dtype=torch.float32).unsqueeze(0)
    len_t = torch.tensor([win_len], dtype=torch.long)

    with torch.no_grad():
        pred_sc = model(x_t, add_t, len_t).numpy()
    result_vals = sc_tgt.inverse_transform(pred_sc)[0]
    b_pred, R0_pred = float(result_vals[0]), float(result_vals[1])

    err_b  = abs(b_pred - b_true) / (b_true + 1e-9) * 100
    err_R0 = abs(R0_pred - R0_true) / (R0_true + 1e-9) * 100
    errors_beta.append(err_b); errors_R0.append(err_R0)

    print(f"  {int(row['sim_idx']):>6}  {b_true:>9.4f}  {R0_true:>7.3f}  "
          f"{b_pred:>9.4f}  {R0_pred:>7.3f}  {err_b:>9.1f}%  {err_R0:>7.1f}%")

print("-" * 70)
print(f"  {'MEAN':>6}  {'':>9}  {'':>7}  "
      f"{'':>9}  {'':>7}  {np.mean(errors_beta):>9.1f}%  {np.mean(errors_R0):>7.1f}%")

print("\nDone. Check that:")
print("  - load_state_dict shows []  (no missing/unexpected keys)")
print("  - Inferred additional_dim matches config (10)")
print("  - Mean beta error is similar to or better than the original BiLSTM")
