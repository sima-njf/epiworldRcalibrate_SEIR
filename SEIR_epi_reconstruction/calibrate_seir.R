#' @keywords internal
#' @importFrom stats setNames
"_PACKAGE"

.bilstm_env <- new.env(parent = emptyenv())
.bilstm_env$loaded <- FALSE

# ── Embedded Python ───────────────────────────────────────────────────────────
.python_code <- '
import torch
import torch.nn as nn
from torch.nn.utils.rnn import pack_padded_sequence
import joblib
import numpy as np
import warnings
warnings.filterwarnings("ignore", category=UserWarning, module="sklearn")
warnings.filterwarnings("ignore", category=FutureWarning, module="torch")

_model = None; _scaler_add = None; _scaler_tgt = None
_device = torch.device("cpu"); T_MAX = 365

class BernardoBiLSTM(nn.Module):
    """
    Bernardo production BiLSTM (Optuna tuned, exp_production.py):
      hidden=64, num_layers=3, bidirectional
      additional_dim=10  [n, recov, incub (scaled),
                          win_len/T_MAX, log_mean/10, log_std/10,
                          shape_logslope, shape_curvature,
                          shape_argmax_frac, shape_late_over_early]
      output_dim=2       [sigmoid(beta), sigmoid(R0)]
    """
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
        rep = torch.cat([h[-2], h[-1]], dim=-1)   # last fwd + last bwd hidden
        z   = self.head_drop(self.act(
                  self.fc1(torch.cat([rep, add_inputs], dim=-1))))
        return torch.sigmoid(self.fc2(z))   # [beta, R0]


def _shape_features(x_raw):
    """Matches DeepIMC/generate_tuned_predictions.py exactly."""
    x_raw = np.asarray(x_raw, dtype=np.float64)
    n     = len(x_raw)
    lc    = np.log1p(x_raw)
    t     = np.arange(n, dtype=np.float64)
    slope = float(np.polyfit(t, lc, 1)[0]) if n > 1 else 0.0
    half  = max(1, n // 2)
    if half > 1 and n - half > 1:
        s1        = float(np.polyfit(t[:half], lc[:half], 1)[0])
        s2        = float(np.polyfit(t[half:], lc[half:], 1)[0])
        curvature = s2 - s1
    else:
        curvature = 0.0
    argmax_frac     = float(np.argmax(x_raw)) / max(1, n - 1)
    q               = max(1, n // 4)
    first_q         = float(x_raw[:q].mean())
    last_q          = float(x_raw[-q:].mean())
    late_over_early = (last_q + 1.0) / (first_q + 1.0)
    return [float(np.tanh(slope)), float(np.tanh(curvature)),
            argmax_frac, float(np.tanh(np.log1p(late_over_early) / 3))]


def load_model(model_path, scaler_add_path, scaler_tgt_path):
    import json, os
    global _model, _scaler_add, _scaler_tgt

    _scaler_add = joblib.load(scaler_add_path)
    _scaler_tgt = joblib.load(scaler_tgt_path)

    # Read architecture from model_config.json sitting next to the .pt file
    cfg_path = os.path.join(os.path.dirname(model_path), "model_config.json")
    with open(cfg_path) as f:
        cfg = json.load(f)["architecture"]

    state  = torch.load(model_path, map_location=_device, weights_only=True)
    fc_dim = state["fc1.weight"].shape[0]
    add_dim = state["fc1.weight"].shape[1] - 2 * cfg["hidden"]

    _model = BernardoBiLSTM(
        hidden         = cfg["hidden"],
        num_layers     = cfg["num_layers"],
        lstm_dropout   = cfg["lstm_dropout"],
        head_dropout   = cfg["head_dropout"],
        additional_dim = add_dim,
        fc_dim         = fc_dim,
    )
    _model.load_state_dict(state)
    _model.to(_device).eval()


def predict(seq, n, recov, incub):
    """
    seq   : list or 1-D array, length 15-365, raw daily incidence (NOT normalised)
    n     : population size
    recov : recovery rate per day (e.g. 1/10)
    incub : incubation period in days (known parameter)
    Returns [beta, R0] in natural units.
    """
    x_raw   = np.asarray(seq, dtype=np.float32)
    win_len = len(x_raw)

    # RevIN: per-window z-score normalisation for LSTM input
    wm, ws = float(x_raw.mean()), float(x_raw.std())
    if ws < 1e-6:
        x_norm   = np.zeros_like(x_raw)
        log_mean, log_std = float(np.log1p(wm)), 0.0
    else:
        x_norm   = (x_raw - wm) / ws
        log_mean = float(np.log1p(wm))
        log_std  = float(np.log1p(ws))

    # Scale population features [n, recov, incub]
    add_sc = _scaler_add.transform(
        np.array([[n, recov, incub]], dtype=np.float32))[0]

    # Full 10-dim additional vector (shape features match DeepIMC exactly)
    add_row = np.array([
        add_sc[0], add_sc[1], add_sc[2],   # scaled n, recov, incub
        win_len / T_MAX,                     # fractional window length
        log_mean / 10.0,                     # log-scale level
        log_std  / 10.0,                     # log-scale spread
    ] + _shape_features(x_raw), dtype=np.float32)

    x_t   = torch.tensor(x_norm,   dtype=torch.float32).unsqueeze(0).unsqueeze(-1)
    add_t = torch.tensor(add_row,  dtype=torch.float32).unsqueeze(0)
    len_t = torch.tensor([win_len], dtype=torch.long)

    with torch.no_grad():
        pred_sc = _model(x_t, add_t, len_t).numpy()

    return _scaler_tgt.inverse_transform(pred_sc)[0].tolist()


def cleanup_model():
    global _model, _scaler_add, _scaler_tgt
    _model = None; _scaler_add = None; _scaler_tgt = None
    return True
'

# ── init_bilstm_model ─────────────────────────────────────────────────────────
#' Load the hypertuned SEIR BiLSTM model
#'
#' @param model_dir Path to the folder containing
#'   \code{model_bilstm_tuned.pt}, \code{scaler_additional.pkl},
#'   \code{scaler_targets.pkl}.
#' @param force_reload Reload even if already loaded.
#' @return Invisibly \code{TRUE}.
#' @export
init_bilstm_model <- function(model_dir, force_reload = FALSE) {
  if (.bilstm_env$loaded && !force_reload) {
    message("Bernardo BiLSTM model already loaded.")
    return(invisible(TRUE))
  }

  model_dir <- normalizePath(model_dir, winslash = "/", mustWork = TRUE)

  # Scaler files: prefer versioned "(5)" names, fall back to unversioned
  .pick_scaler <- function(dir, base) {
    versioned <- file.path(dir, paste0(base, " (5).pkl"))
    plain     <- file.path(dir, paste0(base, ".pkl"))
    if (file.exists(versioned)) return(versioned)
    if (file.exists(plain))     return(plain)
    versioned   # return versioned name so the missing-file check reports it
  }

  files <- list(
    model  = file.path(model_dir, "model_bilstm_tuned.pt"),
    s_add  = .pick_scaler(model_dir, "scaler_additional"),
    s_tgt  = .pick_scaler(model_dir, "scaler_targets")
  )
  missing <- names(files)[!vapply(files, file.exists, logical(1))]
  if (length(missing))
    stop("Missing files in model_dir: ", paste(missing, collapse = ", "),
         call. = FALSE)

  reticulate::py_run_string(.python_code)
  reticulate::py$load_model(files$model, files$s_add, files$s_tgt)

  .bilstm_env$loaded <- TRUE
  message("Bernardo BiLSTM loaded (3-layer, 10 additional features, RevIN, 15-365 days).")
  invisible(TRUE)
}

# ── Informativeness check ─────────────────────────────────────────────────────
#' @keywords internal
.last_informative_day <- function(daily_cases, k = 7, min_count = 1) {
  n <- length(daily_cases)
  if (n < k) {
    return(if (mean(daily_cases) >= min_count) n else 0L)
  }
  roll <- vapply(k:n, function(i) mean(daily_cases[(i - k + 1):i]), numeric(1))
  below <- which(roll < min_count)
  if (length(below) == 0) return(n)
  cut_at <- (k - 1) + below[1]
  max(cut_at, 0L)
}

#' @keywords internal
.is_informative <- function(daily_cases, k = 7, min_count = 1,
                            min_len = 15, min_total_cases = 5) {
  n <- length(daily_cases)
  if (n < min_len) return(FALSE)
  if (sum(daily_cases) < min_total_cases) return(FALSE)
  if (.last_informative_day(daily_cases, k = k, min_count = min_count) <= 0) return(FALSE)
  TRUE
}

# ── calibrate_seir ────────────────────────────────────────────────────────────
#' Predict SEIR parameters from a daily incidence window
#'
#' @param daily_cases Numeric vector, length 15-365, raw daily incidence counts.
#'   Do not normalise — RevIN is applied internally.
#' @param population_size Single numeric: population size (n).
#' @param recovery_rate Single numeric: recovery rate (e.g. \code{1/7}).
#' @param incubation_days Single numeric: incubation period in days (known).
#' @param auto_trim If TRUE (default FALSE), automatically trim daily_cases
#'   to its informative prefix before calibrating.
#' @param warn_uninformative If TRUE (default), warn when input is near-zero.
#'
#' @return Named numeric vector: \code{beta}, \code{R0}.
#' @export
calibrate_seir <- function(daily_cases, population_size,
                           recovery_rate, incubation_days,
                           auto_trim = FALSE, warn_uninformative = TRUE) {
  if (!.bilstm_env$loaded)
    stop("Model not loaded. Call init_bilstm_model() first.", call. = FALSE)

  n <- length(daily_cases)
  if (!is.numeric(daily_cases) || n < 15 || n > 365 || any(daily_cases < 0))
    stop("daily_cases must be numeric with 15-365 non-negative values.", call. = FALSE)
  if (!is.numeric(population_size) || population_size <= 0)
    stop("population_size must be a positive number.", call. = FALSE)
  if (!is.numeric(recovery_rate) || recovery_rate <= 0)
    stop("recovery_rate must be a positive number.", call. = FALSE)
  if (!is.numeric(incubation_days) || incubation_days <= 0)
    stop("incubation_days must be a positive number.", call. = FALSE)

  informative <- .is_informative(daily_cases)

  if (auto_trim) {
    cut <- .last_informative_day(daily_cases)
    cut <- max(cut, 15L)
    cut <- min(cut, n)
    daily_cases <- daily_cases[seq_len(cut)]
    informative <- .is_informative(daily_cases)
  }

  if (!informative && warn_uninformative) {
    warning(
      "daily_cases looks like a near-zero / uninformative window. ",
      "Treat this prediction with caution, or pass auto_trim = TRUE.",
      call. = FALSE
    )
  }

  out <- reticulate::py$predict(
    as.numeric(daily_cases),
    as.numeric(population_size),
    as.numeric(recovery_rate),
    as.numeric(incubation_days)
  )
  names(out) <- c("beta", "R0")
  attr(out, "informative") <- informative
  out
}

# ── cleanup_model ─────────────────────────────────────────────────────────────
#' Unload the model from memory
#' @return Invisibly \code{TRUE}.
#' @export
cleanup_model <- function() {
  if (!.bilstm_env$loaded) { message("No model loaded."); return(invisible(TRUE)) }
  try(reticulate::py$cleanup_model(), silent = TRUE)
  .bilstm_env$loaded <- FALSE
  message("Model unloaded.")
  invisible(TRUE)
}

# ── .onAttach ─────────────────────────────────────────────────────────────────
#' @keywords internal
.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    "epiworldRcalibrate (SEIR) loaded — Bernardo BiLSTM (3-layer, 10 features).\n",
    "  1. init_bilstm_model('path/to/bernardo_model')  # needs model_bilstm_tuned.pt\n",
    "  2. calibrate_seir(daily_cases, population_size, recovery_rate, incubation_days)"
  )
}
