"""
bernardo_predict_oldmodel.py
=============================
Runs the ARCHIVED old BiLSTM (SEIR_epi_reconstruction/model/model_bilstm (1).pt,
architecture + feature pipeline lifted verbatim from
archive_old_model/calibrate_seir.R's embedded Python) on the same real windows
the current bernardo_model/model_bilstm_tuned.pt is evaluated on, so the two
can be compared honestly on identical inputs.

Old model: hidden=160, num_layers=3, additional_dim=6 (no shape features),
targets scaler trained on R0 in [1.0, 5.0] -- i.e. it CANNOT output a
subcritical/decaying epidemic by construction.
"""
import os, json
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import joblib
import warnings
warnings.filterwarnings("ignore")

REAL_DIR  = os.path.dirname(os.path.abspath(__file__))
MODEL_DIR = os.path.join(REAL_DIR, "..", "model")
DEVICE = torch.device("cpu")
T_MAX  = 365

class BiLSTMRegressorOld(nn.Module):
    def __init__(self):
        super().__init__()
        self.lstm = nn.LSTM(input_size=1, hidden_size=160, num_layers=3,
                             batch_first=True, bidirectional=True, dropout=0.5)
        self.fc1 = nn.Linear(2 * 160 + 6, 64)
        self.fc2 = nn.Linear(64, 2)
        self.act = nn.GELU()
        self.head_drop = nn.Dropout(0.5)

    def forward(self, x, add_inputs, lengths):
        from torch.nn.utils.rnn import pack_padded_sequence
        packed = pack_padded_sequence(x, lengths.cpu(), batch_first=True, enforce_sorted=False)
        _, (h, _) = self.lstm(packed)
        rep = torch.cat([h[-2], h[-1]], dim=-1)
        z = self.head_drop(self.act(self.fc1(torch.cat([rep, add_inputs], dim=-1))))
        return torch.sigmoid(self.fc2(z))

def load():
    model = BiLSTMRegressorOld()
    sd = torch.load(os.path.join(MODEL_DIR, "model_bilstm (1).pt"), map_location=DEVICE, weights_only=True)
    model.load_state_dict(sd)
    model.to(DEVICE).eval()
    sa = joblib.load(os.path.join(MODEL_DIR, "scaler_additional (1).pkl"))
    st = joblib.load(os.path.join(MODEL_DIR, "scaler_targets (1).pkl"))
    return model, sa, st

def predict(model, sa, st, seq, n_pop, recov, incub):
    x_raw = np.asarray(seq, dtype=np.float32)
    win_len = len(x_raw)
    wm, ws = float(x_raw.mean()), float(x_raw.std())
    if ws < 1e-6:
        x_norm = np.zeros_like(x_raw); log_mean = float(np.log1p(wm)); log_std = 0.0
    else:
        x_norm = (x_raw - wm) / ws
        log_mean = float(np.log1p(wm)); log_std = float(np.log1p(ws))
    add_sc = sa.transform(np.array([[n_pop, recov, incub]], dtype=np.float32))[0]
    add_row = np.array([add_sc[0], add_sc[1], add_sc[2],
                         win_len / T_MAX, log_mean / 10.0, log_std / 10.0], dtype=np.float32)
    x_t = torch.tensor(x_norm, dtype=torch.float32).unsqueeze(0).unsqueeze(-1)
    add_t = torch.tensor(add_row, dtype=torch.float32).unsqueeze(0)
    len_t = torch.tensor([win_len], dtype=torch.long)
    with torch.no_grad():
        pred_sc = model(x_t, add_t, len_t).numpy()
    beta, R0 = st.inverse_transform(pred_sc)[0].tolist()
    return beta, R0

def main():
    with open(os.path.join(REAL_DIR, "seir_scale_config.json")) as f:
        scale_cfg = json.load(f)
    n_bilstm = float(scale_cfg["n_bilstm"])   # same scaled population the current model uses
    recov = 1/7; incub = 5

    model, sa, st = load()

    print(f"Old model additional-feature trained ranges: n in {list(sa.data_min_[:1])}-{list(sa.data_max_[:1])}, "
          f"recov in {sa.data_min_[1]:.3f}-{sa.data_max_[1]:.3f}, incub in {sa.data_min_[2]:.1f}-{sa.data_max_[2]:.1f}")
    print(f"Old model target trained ranges: beta in [{st.data_min_[0]:.4f}, {st.data_max_[0]:.4f}], "
          f"R0 in [{st.data_min_[1]:.4f}, {st.data_max_[1]:.4f}]  <-- R0 floor is 1.0, cannot predict decay\n")

    for label, csv in [("current61", "utah_covid_current61.csv"),
                        ("wave1_first60", "utah_covid_wave1.csv")]:
        inc_df = pd.read_csv(os.path.join(REAL_DIR, csv))
        inc = inc_df["daily_cases"].values.astype(np.float32)
        if label == "wave1_first60":
            inc = inc[:60]
        n_pop_true = 3.34e6
        scale_factor = min(1.0, n_bilstm / n_pop_true)
        inc_scaled = inc * scale_factor
        beta, R0 = predict(model, sa, st, inc_scaled, n_bilstm, recov, incub)
        print(f"{label}: n={len(inc)}d  beta={beta:.4f}  R0={R0:.4f}")

if __name__ == "__main__":
    main()
