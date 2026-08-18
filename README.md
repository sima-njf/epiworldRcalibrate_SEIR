# epiworldR calibration experiments

This repository compares approaches for recovering epidemic parameters from
agent-based-model (ABM) incidence curves. The main SEIR experiments compare
ABC, ABC-SMC, Nelder-Mead, differential evolution, and a BiLSTM estimator,
then reconstruct incidence with `epiworldR::ModelSEIRCONN` to evaluate the
estimated parameters.

## Repository structure

- `SEIR_epi_reconstruction/`: current SEIR simulation, calibration, comparison,
  and plotting pipeline. Shared model conventions live in `seir_common.R`.
- `SEIR_epi_reconstruction/data_construction/`: synthetic parameter and
  incidence generation.
- `SEIR_epi_reconstruction/bernardo/` and
  `SEIR_epi_reconstruction/bernardo_model/`: BiLSTM inference, trained-model
  artifacts, and BiLSTM comparison outputs.
- `SEIR_epi_reconstruction/real_data/`: Utah COVID-19 and Hagelloch measles data
  preparation, real-data calibration, diagnostics, and plots.
- `SIR_epi_reconstruction/`: earlier SIR reconstruction experiments.
- `previous_seir_epi_reconstructions/` and
  `SEIR_epi_reconstruction/archive_old_model/`: retained legacy implementations;
  these are not the current pipeline.

Several scripts contain cluster-specific paths and SLURM settings. Review those
settings before running them in a different environment.

## Experiment structure

The synthetic experiments draw SEIR parameters, generate daily Exposed-to-
Infected (E-to-I) transition counts, and hold out simulations for evaluation.
Each calibration method observes the same incidence window and estimates the
effective transmission rate

$$
\beta = c p, \qquad R_0 = \frac{\beta}{\gamma},
$$

where $c$ is the contact rate, $p$ is the per-contact transmission probability,
and $\gamma$ is the daily recovery probability. The estimated parameters are
put back into the same ABM, which is run repeatedly to obtain reconstructed
incidence, uncertainty bands, accuracy metrics, and computational cost.

The real-data workflow follows the same comparison:

1. prepare an observed daily-incidence series and fixed model inputs;
2. obtain a BiLSTM prediction for the selected observation window;
3. calibrate the simulation-based methods on that same window;
4. re-simulate every estimate under common initial conditions; and
5. compare the observed curve with the reconstructed ensembles.

The principal entry point is
`SEIR_epi_reconstruction/real_data/calibrate_real_5method.R`.

## Establishing initial prevalence from incidence

Prevalence is a stock, whereas daily incidence is a flow. Consequently, using
the first observed count as `prevalence * N` is generally not a compatible SEIR
initial condition. It also leaves the real-data comparison misaligned because
`ModelSEIRCONN` places all prevalence seeds in Exposed by default, even though
an observed outbreak normally already has both Exposed and Infected people.

Let

- $y_t$ be observed E-to-I incidence on day $t$;
- $N$ be the modeled population size;
- $L$ be the mean latent duration (`incub_days`);
- $\gamma$ be the daily recovery probability, so the mean infectious duration
  is $D = 1/\gamma$; and
- $k$ be a short opening smoothing window (currently $k=3$).

We estimate the opening incidence rate as

$$
\widehat{\lambda}_0 = \frac{1}{k}\sum_{t=1}^{k} y_t.
$$

Assuming incidence is approximately constant over the short initialization
interval, a stock-flow approximation gives

$$
\widehat E_0 = \widehat{\lambda}_0 L,
\qquad
\widehat I_0 = \widehat{\lambda}_0 D
              = \frac{\widehat{\lambda}_0}{\gamma}.
$$

The initial active prevalence and compartment split are therefore

$$
\widehat\pi_0
  = \frac{\widehat E_0 + \widehat I_0}{N}
  = \frac{\widehat{\lambda}_0(L + 1/\gamma)}{N},
$$

$$
q_E = \frac{L}{L + 1/\gamma},
\qquad
q_I = \frac{1/\gamma}{L + 1/\gamma}.
$$

The simulator receives $N\widehat\pi_0$ active seeds and assigns fractions
$q_E$ and $q_I$ to Exposed and Infected. This makes the expected opening E-to-I
flow approximately

$$
\frac{\widehat E_0}{L} = \widehat{\lambda}_0,
$$

so calibration and reconstruction begin on the same incidence scale. The
implementation is `seir_initial_conditions()` in
`SEIR_epi_reconstruction/seir_common.R`. As defensive guards, the implementation
retains at least half the population as Susceptible and uses one Exposed seed
when the opening incidence estimate is zero.

This is an initialization approximation, not an independent prevalence
measurement. If reliable $E_0$ and $I_0$ estimates are available, they should
replace the approximation directly. If only total active prevalence $\pi_0$ is
known, use $N\pi_0$ as the seed total and either use an independently informed
E/I split or use $q_E,q_I$ above. When incidence is under-reported, delayed, or
changing rapidly, the reporting model and uncertainty in $L$, $\gamma$, and the
initial state should be included in sensitivity analyses.
