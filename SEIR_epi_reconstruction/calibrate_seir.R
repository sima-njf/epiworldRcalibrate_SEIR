# =====================================================================
# BiLSTM Epicurve Reconstruction  (v4 model — no hyperparameter file)
# =====================================================================
# Required files in OUTPUT_DIR:
#   model_bilstm.pt         — trained weights
#   scaler_additional.pkl   — MinMaxScaler fitted on [n, recov, incub]
#   scaler_targets.pkl      — MinMaxScaler fitted on [beta, R0]
#
# The model predicts beta (= ptran * crate) and R0.
# incubation days (incub) are a KNOWN input — you must supply them.
# =====================================================================

library(tidyverse)
library(epiworldR)
library(reticulate)
library(data.table)

torch  <- import("torch")
np     <- import("numpy")
joblib <- import("joblib")

# ── Path to your model folder ─────────────────────────────────────────
OUTPUT_DIR <- path.expand("~/sima/epiworldRcalibrate_SEIR/SEIR_epi_reconstruction/model")

cat("Path exists:", dir.exists(OUTPUT_DIR), "\n")
cat("Files found:\n")
print(list.files(OUTPUT_DIR))

# =====================================================================
# 1.  Load BiLSTM model + scalers
#     Architecture hardcoded from notebook (no best_hyperparams.pkl):
#       LSTM_HIDDEN  = 160
#       LSTM_LAYERS  = 3
#       additional_dim = 6   [n, recov, incub, win_len/T_MAX, log_mean/10, log_std/10]
#       output_dim     = 2   [beta, R0]
# =====================================================================

py_run_string(paste0("
import torch
import torch.nn as nn
from torch.nn.utils.rnn import pack_padded_sequence
import joblib, os

OUTPUT_DIR = '", OUTPUT_DIR, "'

# ── BiLSTM architecture (must match training exactly) ────────────────
class BiLSTMRegressor(nn.Module):
    def __init__(self, hidden, num_layers, dropout,
                 additional_dim=6, output_dim=2):
        super().__init__()
        self.lstm = nn.LSTM(
            input_size    = 1,
            hidden_size   = hidden,
            num_layers    = num_layers,
            batch_first   = True,
            bidirectional = True,
            dropout       = (dropout if num_layers > 1 else 0.0)
        )
        self.fc1       = nn.Linear(2 * hidden + additional_dim, 64)
        self.fc2       = nn.Linear(64, output_dim)
        self.act       = nn.GELU()
        self.head_drop = nn.Dropout(dropout)

    def forward(self, x, pad_mask, add_inputs, lengths=None):
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

# ── Hardcoded hyperparameters (from notebook — no Optuna file needed) ─
HIDDEN_DIM = 160
NUM_LAYERS = 3
DROPOUT    = 0.0   # irrelevant at inference; model.eval() disables dropout

model = BiLSTMRegressor(
    hidden         = HIDDEN_DIM,
    num_layers     = NUM_LAYERS,
    dropout        = DROPOUT,
    additional_dim = 6,
    output_dim     = 2
)

state = torch.load(
    os.path.join(OUTPUT_DIR, 'model_bilstm.pt'),
    map_location = 'cpu'
)
model.load_state_dict(state)
model.eval()
print(f'BiLSTM loaded  |  hidden={HIDDEN_DIM}  layers={NUM_LAYERS}')

# ── Scalers ───────────────────────────────────────────────────────────
# scaler_additional : MinMaxScaler fitted on [n, recov, incub]  (3 cols)
# scaler_targets    : MinMaxScaler fitted on [beta, R0]          (2 cols)
scaler_additional = joblib.load(os.path.join(OUTPUT_DIR, 'scaler_additional.pkl'))
scaler_targets    = joblib.load(os.path.join(OUTPUT_DIR, 'scaler_targets.pkl'))
print('Scalers loaded.')
print('Ready.')
"))

cat("Model ready!\n")

# =====================================================================
# 2.  calibrate_seir_bilstm()
#     Input  : observed incidence vector + known epi parameters
#     Returns: data.table with beta, R0, incub
# =====================================================================

calibrate_seir_bilstm <- function(incidence_vec,
                                  n,      # population size
                                  recov,  # recovery rate
                                  incub,  # incubation days (known)
                                  T_MAX = 365) {

  win_len  <- length(incidence_vec)
  win_mean <- mean(incidence_vec)
  win_std  <- sd(incidence_vec)

  # RevIN: per-window z-score normalisation
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

# ── Sequence (RevIN-normalised) ───────────────────────────────────────
x_norm   = np.array([", paste(x_norm, collapse = ", "), "], dtype='float32')
X_seq    = torch.tensor(x_norm).unsqueeze(0).unsqueeze(-1)  # (1, T, 1)
lengths  = torch.tensor([", win_len, "])
pad_mask = torch.zeros(1, ", win_len, ", dtype=torch.bool)  # no padding

# ── Additional features (6) ───────────────────────────────────────────
#   scaler_additional was fitted on [n, recov, incub]
add_raw = np.array([[", n, ", ", recov, ", ", incub, "]], dtype='float32')
add_sc  = scaler_additional.transform(add_raw)              # (1, 3)

add_features = np.array([[
    float(add_sc[0, 0]),        # n      (MinMax scaled)
    float(add_sc[0, 1]),        # recov  (MinMax scaled)
    float(add_sc[0, 2]),        # incub  (MinMax scaled)
    ", win_len / T_MAX, ",      # win_len / T_MAX
    ", log_mean / 10.0, ",      # log(1 + mean) / 10
    ", log_std  / 10.0, ",      # log(1 + std)  / 10
]], dtype='float32')

add_tensor = torch.tensor(add_features)

# ── Forward pass ──────────────────────────────────────────────────────
with torch.no_grad():
    pred_scaled = model(X_seq, pad_mask, add_tensor, lengths=lengths).numpy()

pred_nat = scaler_targets.inverse_transform(pred_scaled)

beta_hat  = float(pred_nat[0, 0])   # beta = ptran * crate
R0_hat    = float(pred_nat[0, 1])   # R0
"))

  data.table(
    beta  = py$beta_hat,
    R0    = py$R0_hat,
    incub = incub           # pass-through — incub was input, not predicted
  )
}

# =====================================================================
# 3.  simulate_seir()
#     Runs ModelSEIRCONN via epiworldR and returns the median
#     daily E→I incidence trajectory across nsims replicates.
#
#     Beta decomposition:
#       BiLSTM predicts beta = ptran * crate (the product only).
#       ModelSEIRCONN needs them separately.
#       We fix crate = 1.0  →  ptran = beta.
#       Any factorisation that keeps ptran * crate = beta gives
#       the same mean epidemic curve.
# =====================================================================

simulate_seir <- function(ptran, crate, incub, n_days, N, recov, prev,
                          nsims = 100) {
  model_epi <- ModelSEIRCONN(
    name              = "BiLSTM reconstruction",
    n                 = N,
    prevalence        = prev,
    contact_rate      = crate,
    transmission_rate = ptran,
    incubation_days   = incub,
    recovery_rate     = recov
  )

  saver <- make_saver("transition")

  run_multiple(
    model_epi,
    ndays    = n_days,
    nsims    = nsims,
    saver    = saver,
    nthreads = 8
  )

  sim_results <- run_multiple_get_results(model_epi, nthreads = 8)

  # Median daily E→I transitions across all replicates
  sim_results$transition %>%
    filter(from == "Exposed", to == "Infected") %>%
    group_by(sim_num, date) %>%
    summarise(counts = sum(counts), .groups = "drop") %>%
    group_by(date) %>%
    summarise(median_count = median(counts), .groups = "drop") %>%
    arrange(date) %>%
    pull(median_count)
}

# =====================================================================
# 4.  reconstruct_epicurve()
#     Full pipeline: calibrate → simulate → plot → return results
# =====================================================================

reconstruct_epicurve <- function(
    incidence_vec,        # numeric vector — observed daily incidence
    n,                    # integer        — population size
    recov,                # numeric        — recovery rate (e.g. 1/7 = 0.143)
    prevalence,           # numeric        — initial infected fraction (e.g. 1/N)
    incub,                # numeric        — incubation period in days (known input)
    n_days      = NULL,   # integer        — simulation length; defaults to length(incidence_vec)
    T_MAX       = 365,    # integer        — epidemic horizon used during training
    nsims       = 100     # integer        — stochastic replicates for the simulation
) {

  # ── Step 1: Predict beta and R0 with BiLSTM ──────────────────────
  cal <- calibrate_seir_bilstm(
    incidence_vec = incidence_vec,
    n             = n,
    recov         = recov,
    incub         = incub,
    T_MAX         = T_MAX
  )

  beta_hat  <- cal$beta
  R0_hat    <- cal$R0
  incub_hat <- cal$incub

  # ── Step 2: Decompose beta → ptran, crate ────────────────────────
  ptran_hat <- beta_hat     # crate fixed at 1.0  →  ptran = beta / 1
  crate_hat <- 1.0

  # Consistency check: R0 should ≈ beta / recov
  R0_check <- beta_hat / recov

  cat(sprintf("\n── BiLSTM calibration ───────────────────────────────\n"))
  cat(sprintf("  beta   = %.6f\n",  beta_hat))
  cat(sprintf("  R0     = %.4f  (model)  |  %.4f  (beta/recov)  |  diff = %.4f\n",
              R0_hat, R0_check, abs(R0_hat - R0_check)))
  cat(sprintf("  incub  = %.2f days\n", incub_hat))
  cat(sprintf("  → ptran = %.6f  |  crate = %.2f\n", ptran_hat, crate_hat))

  # ── Step 3: Run SEIR simulation ──────────────────────────────────
  sim_days  <- if (is.null(n_days)) length(incidence_vec) else n_days

  simulated <- simulate_seir(
    ptran  = ptran_hat,
    crate  = crate_hat,
    incub  = incub_hat,
    n_days = sim_days,
    N      = n,
    recov  = recov,
    prev   = prevalence,
    nsims  = nsims
  )

  simulated <- simulated[-1]   # drop day 0 (epiworldR always includes it)

  # ── Step 4: Align lengths ─────────────────────────────────────────
  min_len   <- min(length(simulated), length(incidence_vec))
  simulated <- simulated[1:min_len]
  observed  <- incidence_vec[1:min_len]

  # ── Step 5: MAE ──────────────────────────────────────────────────
  mae_per_day <- abs(simulated - observed)
  mae         <- mean(mae_per_day)

  cat(sprintf("  MAE    = %.2f cases/day  (min %.2f  max %.2f)\n\n",
              mae, min(mae_per_day), max(mae_per_day)))

  # ── Step 6: Plot ──────────────────────────────────────────────────
  plot_df <- tibble(
    day       = 1:min_len,
    observed  = observed,
    simulated = simulated
  )

  subtitle_txt <- sprintf(
    "N = %d  |  recov = %.4f  |  incub = %.1f days  |  prevalence = %.6f\nBiLSTM → beta = %.5f  |  R0 = %.2f  |  MAE = %.2f cases/day",
    n, recov, incub_hat, prevalence,
    beta_hat, R0_hat, mae
  )

  p <- ggplot(plot_df, aes(x = day)) +
    geom_line(aes(y = observed,  color = "Observed"),
              linewidth = 1.2) +
    geom_point(aes(y = observed, color = "Observed"),
               size = 2.5, alpha = 0.7) +
    geom_line(aes(y = simulated, color = "Reconstructed (BiLSTM)"),
              linewidth = 1.1, linetype = "dashed") +
    scale_color_manual(
      values = c("Observed"               = "steelblue",
                 "Reconstructed (BiLSTM)" = "tomato")
    ) +
    labs(
      title    = "Epidemic Curve Reconstruction — BiLSTM",
      subtitle = subtitle_txt,
      x        = "Day",
      y        = "Daily incidence (E\u2192I transitions)",
      color    = ""
    ) +
    theme_minimal(base_size = 13) +
    theme(
      legend.position = "bottom",
      plot.subtitle   = element_text(size = 9, family = "mono")
    )

  print(p)

  # ── Return all results invisibly ─────────────────────────────────
  invisible(list(
    beta      = beta_hat,
    R0        = R0_hat,
    incub     = incub_hat,
    ptran     = ptran_hat,
    crate     = crate_hat,
    simulated = simulated,
    observed  = observed,
    mae       = mae,
    plot      = p
  ))
}

# =====================================================================
# 5.  Run
#     Replace the values below with your actual epidemic data.
# =====================================================================

incidence_vec <- c(
  103, 37, 60, 74, 108, 125, 138, 186, 215, 276, 318, 331, 414, 402, 446,
  454, 405, 401, 373, 334, 285, 241, 219, 156, 140, 108,  93,  82,  73,  78,
  48,  38,  34,  22,  27,  22,  20,  11,  14,   8,  11,  14,   8,   6,   6,
  2,   0,   7,   2,   4,   4,   1,   1,   2,   2,   1,   0,   2,   1,   1, 0
)

n_pop      <- 7087
recov_rate <- 0.203
incub_days <- 5.0          # <-- set to your known incubation period
prevalence <- 1 / n_pop   # initial fraction infected; adjust as needed

result <- reconstruct_epicurve(
  incidence_vec = incidence_vec,
  n             = n_pop,
  recov         = recov_rate,
  prevalence    = prevalence,
  incub         = incub_days,
  nsims         = 100
)

cat("══════════════════════════════════════════\n")
cat(sprintf("  beta   : %.6f\n",    result$beta))
cat(sprintf("  R0     : %.4f\n",    result$R0))
cat(sprintf("  incub  : %.2f days\n", result$incub))
cat(sprintf("  MAE    : %.2f cases/day\n", result$mae))
cat("══════════════════════════════════════════\n")
