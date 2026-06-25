#' @keywords internal
#' @importFrom stats setNames
"_PACKAGE"

.bilstm_env <- new.env(parent = emptyenv())
.bilstm_env$loaded <- FALSE

# ── Embedded Python ───────────────────────────────────────────────────────────
.python_code <- '
import torch, torch.nn as nn, torch.nn.functional as F
from torch.nn.utils.rnn import pack_padded_sequence
import joblib, numpy as np, warnings
warnings.filterwarnings("ignore", category=UserWarning, module="sklearn")

_model = None; _scaler_add = None; _scaler_tgt = None
_device = torch.device("cpu"); T_MAX = 365

class BiLSTMRegressor(nn.Module):
    def __init__(self):
        super().__init__()
        self.lstm = nn.LSTM(input_size=1, hidden_size=160, num_layers=3,
                            batch_first=True, bidirectional=True, dropout=0.5)
        self.fc1       = nn.Linear(2 * 160 + 5, 64)
        self.fc2       = nn.Linear(64, 3)
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
        raw = self.fc2(z)
        return torch.stack([torch.sigmoid(raw[:, 0]),
                            F.softplus(raw[:, 1]),
                            F.softplus(raw[:, 2])], dim=1)

def load_model(model_path, scaler_add_path, scaler_tgt_path):
    global _model, _scaler_add, _scaler_tgt
    _scaler_add = joblib.load(scaler_add_path)
    _scaler_tgt = joblib.load(scaler_tgt_path)
    _model = BiLSTMRegressor()
    _model.load_state_dict(torch.load(model_path, map_location=_device))
    _model.to(_device).eval()

def predict(seq, n, recov):
    x_raw = np.asarray(seq, dtype=np.float32)
    wm, ws = float(x_raw.mean()), float(x_raw.std())
    if ws < 1e-6:
        x_norm, log_mean, log_std = np.zeros_like(x_raw), float(np.log1p(wm)), 0.0
    else:
        x_norm = (x_raw - wm) / ws
        log_mean, log_std = float(np.log1p(wm)), float(np.log1p(ws))
    add_sc  = _scaler_add.transform(np.array([[n, recov]], dtype=np.float32))[0]
    add_row = np.array([add_sc[0], add_sc[1],
                        len(x_raw) / T_MAX, log_mean / 10.0, log_std / 10.0],
                       dtype=np.float32)
    x_t   = torch.tensor(x_norm,  dtype=torch.float32).unsqueeze(0).unsqueeze(-1)
    add_t = torch.tensor(add_row, dtype=torch.float32).unsqueeze(0)
    len_t = torch.tensor([len(x_raw)], dtype=torch.long)
    mask_t = torch.zeros(1, len(x_raw), dtype=torch.bool)
    with torch.no_grad():
        out = _model(x_t, mask_t, add_t, lengths=len_t).numpy()
    return _scaler_tgt.inverse_transform(out)[0].tolist()
'

# ── init_bilstm_model ─────────────────────────────────────────────────────────
#' Load the BiLSTM model
#'
#' @param model_dir Path to the folder containing \code{model_bilstm_sir.pt},
#'   \code{scaler_additional.pkl}, and \code{scaler_targets.pkl}.
#' @param force_reload Reload even if already loaded.
#' @return Invisibly \code{TRUE}.
#' @export
init_bilstm_model <- function(model_dir, force_reload = FALSE) {
  if (.bilstm_env$loaded && !force_reload) {
    message("BiLSTM model already loaded.")
    return(invisible(TRUE))
  }

  model_dir <- normalizePath(model_dir, winslash = "/", mustWork = TRUE)
  files <- list(
    model  = file.path(model_dir, "model_bilstm_sir.pt"),
    s_add  = file.path(model_dir, "scaler_additional.pkl"),
    s_tgt  = file.path(model_dir, "scaler_targets.pkl")
  )
  missing <- names(files)[!vapply(files, file.exists, logical(1))]
  if (length(missing))
    stop("Missing files: ", paste(missing, collapse = ", "), call. = FALSE)

  reticulate::py_run_string(.python_code)
  reticulate::py$load_model(files$model, files$s_add, files$s_tgt)

  .bilstm_env$loaded <- TRUE
  message("BiLSTM model loaded.")
  invisible(TRUE)
}

# ── calibrate_sir ─────────────────────────────────────────────────────────────
#' Predict SIR parameters from a daily incidence window
#'
#' @param daily_cases Numeric vector, length 15-365, raw daily incidence counts.
#'   Do not normalise — RevIN is applied internally.
#' @param population_size Single numeric: population size (n).
#' @param recovery_rate Single numeric: recovery rate (e.g. \code{1/7}).
#'
#' @return Named numeric vector: \code{ptran}, \code{crate}, \code{R0}.
#'   \code{crate} is back-calculated as \code{R0 * recovery_rate / ptran}
#'   to enforce SIR physics exactly.
#' @export
calibrate_sir <- function(daily_cases, population_size, recovery_rate) {
  if (!.bilstm_env$loaded)
    stop("Model not loaded. Call init_bilstm_model() first.", call. = FALSE)

  n <- length(daily_cases)
  if (!is.numeric(daily_cases) || n < 15 || n > 365 || any(daily_cases < 0))
    stop("daily_cases must be a numeric vector with 15-365 non-negative values.",
         call. = FALSE)
  if (!is.numeric(population_size) || population_size <= 0)
    stop("population_size must be a positive number.", call. = FALSE)
  if (!is.numeric(recovery_rate) || recovery_rate <= 0)
    stop("recovery_rate must be a positive number.", call. = FALSE)

  out <- reticulate::py$predict(as.numeric(daily_cases),
                                as.numeric(population_size),
                                as.numeric(recovery_rate))
  names(out) <- c("ptran", "crate", "R0")
  out[["crate"]] <- out[["R0"]] * recovery_rate / out[["ptran"]]
  out
}

# ── .onAttach ─────────────────────────────────────────────────────────────────
#' @keywords internal
.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    "epiworldRcalibrate loaded.\n",
    "  1. init_bilstm_model(\'path/to/model_output_revin_sir\')\n",
    "  2. calibrate_sir(daily_cases, population_size, recovery_rate)"
  )
}
