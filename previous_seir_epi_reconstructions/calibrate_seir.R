library(reticulate)
library(data.table)

torch  <- import("torch")
np     <- import("numpy")
joblib <- import("joblib")

# ── 1. Path ───────────────────────────────────────────────────────────
OUTPUT_DIR <- path.expand("~/epiworldRcalibrate_SEIR/model")

cat("Path exists:", dir.exists(OUTPUT_DIR), "\n")
cat("Files:\n"); print(list.files(OUTPUT_DIR))

# ── 2. Load model + scalers ───────────────────────────────────────────
py_run_string(paste0("
import torch
import torch.nn as nn
from torch.nn.utils.rnn import pack_padded_sequence
import joblib, os

OUTPUT_DIR = '", OUTPUT_DIR, "'

# ── Architecture (must match training exactly) ────────────────────────
#   BiLSTM → GELU → fc1(2*hidden+5, 64) → GELU → fc2(64, 3) → sigmoid
#   additional_dim = 5: [n, recov, win_len/T_MAX, log_mean/10, log_std/10]
#   output: [beta, R0, incub]  all in scaled [0,1] space, inverted after

class BiLSTMRegressor(nn.Module):
    def __init__(self, hidden, num_layers, dropout,
                 additional_dim=5, output_dim=3):
        super().__init__()
        self.lstm = nn.LSTM(
            input_size=1,
            hidden_size=hidden,
            num_layers=num_layers,
            batch_first=True,
            bidirectional=True,
            dropout=(dropout if num_layers > 1 else 0.0),
        )
        self.fc1       = nn.Linear(2 * hidden + additional_dim, 64)
        self.fc2       = nn.Linear(64, output_dim)
        self.act       = nn.GELU()
        self.head_drop = nn.Dropout(dropout)

    def forward(self, x, add_inputs, lengths=None):
        if lengths is not None:
            packed = pack_padded_sequence(
                x, lengths.cpu(), batch_first=True, enforce_sorted=False)
            _, (h, _) = self.lstm(packed)
        else:
            _, (h, _) = self.lstm(x)
        rep = torch.cat([h[-2], h[-1]], dim=-1)
        z   = torch.cat([rep, add_inputs], dim=-1)
        z   = self.head_drop(self.act(self.fc1(z)))
        return torch.sigmoid(self.fc2(z))

# ── Load hyperparameters ──────────────────────────────────────────────
hp         = joblib.load(os.path.join(OUTPUT_DIR, 'best_hyperparams.pkl'))
hidden_dim = hp['hidden_dim']
num_layers = hp['num_layers']
dropout    = hp['dropout']
print(f'Hyperparams: hidden={hidden_dim}, layers={num_layers}, dropout={dropout:.4f}')

# ── Build model and load weights ──────────────────────────────────────
model = BiLSTMRegressor(
    hidden     = hidden_dim,
    num_layers = num_layers,
    dropout    = dropout,
)
state = torch.load(
    os.path.join(OUTPUT_DIR, 'model_bilstm_revin_v3.pt'),
    map_location='cpu',
)
model.load_state_dict(state)
model.eval()
print('Model loaded.')

# ── Load scalers ──────────────────────────────────────────────────────
# scaler_additional: fitted on [n, recov]       (2 columns)
# scaler_targets   : fitted on [beta, R0, incub] (3 columns)
scaler_additional = joblib.load(os.path.join(OUTPUT_DIR, 'scaler_additional.pkl'))
scaler_targets    = joblib.load(os.path.join(OUTPUT_DIR, 'scaler_targets.pkl'))
print('Scalers loaded.')
print('Model and scalers ready!')
"))

cat("Model ready!\n")

# ── 3. Inference function ─────────────────────────────────────────────
#
# Parameters:
#   incidence_vec : numeric vector of daily incidence counts (any length)
#   n             : population size
#   recov         : recovery rate gamma  (e.g. 1/7 ≈ 0.143)
#   T_MAX         : epidemic horizon used in training (default 365)
#
# Returns a data.table with columns:
#   beta  = ptran * crate  (overall transmission rate)
#   R0    = basic reproduction number
#   incub = incubation period (days)

calibrate_seir_revin <- function(incidence_vec, n, recov, T_MAX = 365) {

  win_len  <- length(incidence_vec)
  win_mean <- mean(incidence_vec)
  win_std  <- sd(incidence_vec)

  # RevIN: per-window z-score (no global incidence scaler)
  if (win_std < 1e-6) {
    x_norm   <- rep(0.0, win_len)
    log_mean <- log1p(win_mean)
    log_std  <- 0.0
  } else {
    x_norm   <- (incidence_vec - win_mean) / win_std
    log_mean <- log1p(win_mean)
    log_std  <- log1p(win_std)
  }

  py_run_string(paste0("
import numpy as np
import torch

# ── RevIN-normalised sequence ─────────────────────────────────────────
x_norm  = np.array([", paste(x_norm,   collapse = ", "), "], dtype='float32')
X_seq   = torch.tensor(x_norm).unsqueeze(0).unsqueeze(-1)  # (1, T, 1)
lengths = torch.tensor([", win_len, "])

# ── Additional features (5) ───────────────────────────────────────────
#   scaler_additional was fitted on [n, recov] only;
#   the remaining 3 features are computed directly (no scaler needed)
add_raw = np.array([[", n, ", ", recov, "]], dtype='float32')
add_sc  = scaler_additional.transform(add_raw)   # shape (1, 2)

add_features = np.array([[
    float(add_sc[0, 0]),      # n      (MinMax scaled)
    float(add_sc[0, 1]),      # recov  (MinMax scaled)
    ", win_len / T_MAX, ",    # win_len / T_MAX
    ", log_mean / 10.0, ",    # log(1 + mean) / 10
    ", log_std  / 10.0, ",    # log(1 + std)  / 10
]], dtype='float32')

add_tensor = torch.tensor(add_features)

# ── Predict ───────────────────────────────────────────────────────────
with torch.no_grad():
    pred_scaled = model(X_seq, add_tensor, lengths=lengths).numpy()

pred_nat = scaler_targets.inverse_transform(pred_scaled)

beta_hat  = float(pred_nat[0, 0])   # beta  = ptran * crate
R0_hat    = float(pred_nat[0, 1])   # R0
incub_hat = float(pred_nat[0, 2])   # incubation period (days)
  "))

  data.table(
    beta  = py$beta_hat,
    R0    = py$R0_hat,
    incub = py$incub_hat
  )
}

# ── 4. Example usage ──────────────────────────────────────────────────
incidence_vec <- c(
  103, 37, 60, 74, 108, 125, 138, 186, 215, 276, 318, 331, 414, 402, 446,
  454, 405, 401, 373, 334, 285, 241, 219, 156, 140, 108, 93, 82, 73, 78,
  48, 38, 34, 22, 27, 22, 20, 11, 14, 8, 11, 14, 8, 6, 6, 2, 0, 7, 2, 4,
  4, 1, 1, 2, 2, 1, 0, 2, 1, 1, 0
)

recov_rate <- 0.203
n_pop      <- 7087

result <- calibrate_seir_revin(
  incidence_vec = incidence_vec,
  n             = n_pop,
  recov         = recov_rate
)

cat("==============================================\n")
cat("Predicted Parameters (RevIN v3 model):\n")
cat("==============================================\n")
cat(sprintf("  beta  : %.6f  (= ptran * crate)\n", result$beta))
cat(sprintf("  R0    : %.4f\n",  result$R0))
cat(sprintf("  incub : %.4f days\n", result$incub))

# ── 5. SEIR consistency check  (R0 = beta / recov) ───────────────────
R0_calc <- result$beta / recov_rate
R0_diff  <- abs(result$R0 - R0_calc)

cat("\nSEIR Consistency Check (R0 = beta / recov):\n")
cat(sprintf("  R0 from model  : %.4f\n", result$R0))
cat(sprintf("  beta / recov   : %.4f\n", R0_calc))
cat(sprintf("  Difference     : %.6f\n", R0_diff))
if (R0_diff < 0.05) {
  cat("  CHECK PASSED — predictions are SEIR-consistent\n")
} else {
  cat("  CHECK FAILED — consider increasing LAMBDA_PHYS and retraining\n")
}
