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

class BiLSTMRegressor(nn.Module):
    """
    Exact SEIR architecture from SEIR_FourModels_RevIN notebook:
      hidden=160, num_layers=3, dropout=0.5
      additional_dim=6  [n, recov, incub, win_len/T_MAX, log_mean/10, log_std/10]
      output_dim=2      [sigmoid(beta), sigmoid(R0)]
    """
    def __init__(self):
        super().__init__()
        self.lstm = nn.LSTM(input_size=1, hidden_size=160, num_layers=3,
                            batch_first=True, bidirectional=True, dropout=0.5)
        self.fc1       = nn.Linear(2 * 160 + 6, 64)
        self.fc2       = nn.Linear(64, 2)
        self.act       = nn.GELU()
        self.head_drop = nn.Dropout(0.5)

    def forward(self, x, pad_mask, add_inputs, lengths=None):
        if lengths is not None:
            packed = pack_padded_sequence(x, lengths.cpu(),
                                         batch_first=True, enforce_sorted=False)
            _, (h, _) = self.lstm(packed)
        else:
            _, (h, _) = self.lstm(x)
        rep = torch.cat([h[-2], h[-1]], dim=-1)
        z   = self.head_drop(self.act(self.fc1(torch.cat([rep, add_inputs], dim=-1))))
        return torch.sigmoid(self.fc2(z))   # [beta, R0]


def load_model(model_path, scaler_add_path, scaler_tgt_path):
    global _model, _scaler_add, _scaler_tgt
    _scaler_add = joblib.load(scaler_add_path)   # fitted on [n, recov, incub]
    _scaler_tgt = joblib.load(scaler_tgt_path)   # fitted on [beta, R0]
    _model = BiLSTMRegressor()
    _model.load_state_dict(torch.load(model_path, map_location=_device, weights_only=True))
    _model.to(_device).eval()


def predict(seq, n, recov, incub):
    """
    seq   : list or 1-D array, length 15-365, raw daily incidence (NOT normalised)
    n     : population size
    recov : recovery rate (e.g. 1/7)
    incub : incubation period in days (known parameter)
    Returns [beta, R0] in natural units.
    """
    x_raw   = np.asarray(seq, dtype=np.float32)
    win_len = len(x_raw)

    # RevIN: z-score by this window own mean and std
    wm, ws = float(x_raw.mean()), float(x_raw.std())
    if ws < 1e-6:
        x_norm   = np.zeros_like(x_raw)
        log_mean = float(np.log1p(wm))
        log_std  = 0.0
    else:
        x_norm   = (x_raw - wm) / ws
        log_mean = float(np.log1p(wm))
        log_std  = float(np.log1p(ws))

    # Scale [n, recov, incub] with the fitted MinMaxScaler
    add_sc = _scaler_add.transform(
        np.array([[n, recov, incub]], dtype=np.float32))[0]

    # 6-feature additional vector
    add_row = np.array([
        add_sc[0],         # n      (MinMax scaled)
        add_sc[1],         # recov  (MinMax scaled)
        add_sc[2],         # incub  (MinMax scaled) — known parameter
        win_len / T_MAX,   # relative window length
        log_mean / 10.0,   # RevIN level
        log_std  / 10.0,   # RevIN scale
    ], dtype=np.float32)

    x_t    = torch.tensor(x_norm,  dtype=torch.float32).unsqueeze(0).unsqueeze(-1)
    add_t  = torch.tensor(add_row, dtype=torch.float32).unsqueeze(0)
    mask_t = torch.zeros(1, win_len, dtype=torch.bool)
    len_t  = torch.tensor([win_len], dtype=torch.long)

    with torch.no_grad():
        pred_sc = _model(x_t, mask_t, add_t, lengths=len_t).numpy()

    return _scaler_tgt.inverse_transform(pred_sc)[0].tolist()


def cleanup_model():
    global _model, _scaler_add, _scaler_tgt
    _model = None; _scaler_add = None; _scaler_tgt = None
    return True
'

# ── init_bilstm_model ─────────────────────────────────────────────────────────
#' Load the SEIR BiLSTM model
#'
#' @param model_dir Path to the \code{model_output_revin_v3} folder containing
#'   \code{model_bilstm.pt}, \code{scaler_additional.pkl},
#'   \code{scaler_targets.pkl}.
#' @param force_reload Reload even if already loaded.
#' @return Invisibly \code{TRUE}.
#' @export
init_bilstm_model <- function(model_dir, force_reload = FALSE) {
  if (.bilstm_env$loaded && !force_reload) {
    message("SEIR BiLSTM model already loaded.")
    return(invisible(TRUE))
  }

  model_dir <- normalizePath(model_dir, winslash = "/", mustWork = TRUE)
  files <- list(
    model  = file.path(model_dir, "model_bilstm (1).pt"),
    s_add  = file.path(model_dir, "scaler_additional (1).pkl"),
    s_tgt  = file.path(model_dir, "scaler_targets (1).pkl")
  )
  missing <- names(files)[!vapply(files, file.exists, logical(1))]
  if (length(missing))
    stop("Missing files: ", paste(missing, collapse = ", "), call. = FALSE)

  reticulate::py_run_string(.python_code)
  reticulate::py$load_model(files$model, files$s_add, files$s_tgt)

  .bilstm_env$loaded <- TRUE
  message("SEIR BiLSTM model loaded (RevIN, variable-length 15-365 days).")
  invisible(TRUE)
}

# ── Informativeness check (zero-tail fix, R side) ──────────────────────────────
# Mirrors the Python informative_window_utils.py used to fix the training
# pipeline. Rule: 7-day rolling mean must stay >= min_count. This is applied
# here, at the deployment entry point, because this is where real (possibly
# zero-tail) curves actually get passed in by callers.

#' Find the last "informative" day of an incidence curve
#' @keywords internal
.last_informative_day <- function(daily_cases, k = 7, min_count = 1) {
  n <- length(daily_cases)
  if (n < k) {
    return(if (mean(daily_cases) >= min_count) n else 0L)
  }
  roll <- vapply(k:n, function(i) mean(daily_cases[(i - k + 1):i]), numeric(1))
  below <- which(roll < min_count)
  if (length(below) == 0) return(n)
  # roll[j] corresponds to the window ending at day (k - 1 + j)
  cut_at <- (k - 1) + below[1]
  max(cut_at, 0L)
}

#' Is this incidence window informative enough to calibrate from?
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
#'   to its informative prefix before calibrating, using the same rule as
#'   training (7-day rolling mean, threshold 1 case). If the curve is
#'   uninformative even after trimming to the minimum length, the call
#'   proceeds but a warning is issued -- the prediction should not be trusted.
#' @param warn_uninformative If TRUE (default), issue a warning when the
#'   input is a near-zero / uninformative window, regardless of auto_trim.
#'
#' @return Named numeric vector: \code{beta}, \code{R0}. If the input was
#'   flagged as uninformative, this is returned with an additional
#'   \code{"informative"} attribute set to \code{FALSE} so callers can check
#'   \code{attr(out, "informative")} programmatically instead of parsing
#'   warnings.
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

  # ── zero-tail check ──────────────────────────────────────────────────────
  informative <- .is_informative(daily_cases)

  if (auto_trim) {
    cut <- .last_informative_day(daily_cases)
    cut <- max(cut, 15L)                 # never go below the model's min length
    cut <- min(cut, n)
    daily_cases <- daily_cases[seq_len(cut)]
    informative <- .is_informative(daily_cases)
  }

  if (!informative && warn_uninformative) {
    warning(
      "daily_cases looks like a near-zero / uninformative window ",
      "(insufficient recent case counts). The model was trained to be less ",
      "biased on such windows, but a truly zero-information input still ",
      "carries no real signal about R0 -- treat this prediction with caution, ",
      "or pass auto_trim = TRUE to calibrate from the informative prefix only.",
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
    "epiworldRcalibrate (SEIR) loaded.\n",
    "  1. init_bilstm_model('path/to/model_output_revin_v3')\n",
    "  2. calibrate_seir(daily_cases, population_size, recovery_rate, incubation_days)"
  )
}
