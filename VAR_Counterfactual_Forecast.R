#=========================================================================================================================================================
# Nafiu Bashir A. (PhD)                                                                                                                                  #
# Central Bank of Nigeria (CBN)                                                                                                                          #
# Monetary Policy Department                                                                                                                             #
# Email: nafiu13bashir@gmail.com                                                                                                                         #
#=========================================================================================================================================================
#
# POLICY COUNTERFACTUAL ANALYSIS AND FORECAST SCENARIOS
# South African Quarterly Data: 1981Q2 - 2012Q4
#
# Two analyses:
#   (1) COUNTERFACTUAL  -- What would inflation & output have looked like
#                          if there had been NO monetary policy shocks?
#                          Identification: Recursive (Cholesky) SVAR
#                          Ordering: gdp_gap -> inf -> int
#
#   (2) FORECAST SCENARIOS -- Inflation & output forecasts under:
#                             - Baseline (unrestricted VAR)
#                             - Moderate tightening  (+50 bp/Q for 4 quarters)
#                             - Aggressive tightening (+100 bp/Q for 4 quarters)
#=========================================================================================================================================================

rm(list = ls())
graphics.off()

install.packages(c("vars", "mFilter", "ggplot2", "dplyr", "tidyr", "gridExtra"))

#-------------------------------------------------#
#  Libraries                                      #
#-------------------------------------------------#
library(vars)
library(mFilter)
library(ggplot2)
library(dplyr)
library(tidyr)

# gridExtra for multi-panel layouts (install if missing)
if (!requireNamespace("gridExtra", quietly = TRUE)) install.packages("gridExtra")
library(gridExtra)


#=========================================================================================
# 1.  LOAD AND PREPARE DATA
#=========================================================================================

# Set working directory if needed
# setwd("C:/path/to/VARs")

data_sa <- read.csv("data_sa.csv")
nT      <- nrow(data_sa)   # 127 observations

# Time series objects (quarterly, start 1981 Q2)
gdp_raw <- ts(log(data_sa$gdp), start = c(1981, 2), frequency = 4)
inf     <- ts(data_sa$inf,      start = c(1981, 2), frequency = 4)
int     <- ts(data_sa$int,      start = c(1981, 2), frequency = 4)

# HP-filter real GDP to obtain output gap (lambda = 1600 for quarterly data)
gdp_gap <- hpfilter(gdp_raw, freq = 1600)$cycle

# Combine into VAR data matrix
# Ordering matters for Cholesky ID: gdp_gap first (most exogenous),
# inflation second, interest rate last (policy instrument)
dat_var <- cbind(gdp_gap, inf, int)
colnames(dat_var) <- c("gdp_gap", "inf", "int")

cat("=== Data summary ===\n")
print(summary(dat_var))
plot(dat_var, main = "GDP Gap, Inflation and Interest Rate")


#=========================================================================================
# 2.  ESTIMATE TRIVARIATE VAR
#=========================================================================================

info_var <- VARselect(dat_var, lag.max = 8, type = "const")
cat("\n=== Lag selection criteria ===\n")
print(info_var$selection)

p_sel <- as.integer(info_var$selection["SC(n)"])   # BIC/SC lag order
p_sel <- max(p_sel, 2)                              # minimum 2 lags
cat(sprintf("\nUsing p = %d lags\n", p_sel))

bv_est <- VAR(dat_var, p = p_sel, type = "const")
summary(bv_est)

#-- Diagnostics --
cat("\n=== Serial correlation test ===\n")
print(serial.test(bv_est, lags.pt = 12, type = "PT.asymptotic"))

cat("\n=== Stability (eigenvalues, all must be < 1) ===\n")
print(roots(bv_est))


#=========================================================================================
# 3.  SVAR IDENTIFICATION: RECURSIVE (CHOLESKY)
#     Ordering: gdp_gap -> inf -> int
#     The monetary policy (MP) shock is the orthogonalized residual in the
#     interest rate equation AFTER controlling for gdp_gap and inflation.
#=========================================================================================

resid_rf <- residuals(bv_est)           # (T-p) x 3 reduced-form residuals
T_resid  <- nrow(resid_rf)             # number of residuals = nT - p_sel

# Cholesky decomposition of residual covariance matrix
# Sigma = A0 %*% t(A0),  A0 is lower-triangular
Sigma <- crossprod(resid_rf) / T_resid
A0    <- t(chol(Sigma))                 # lower-triangular impact matrix

# Structural shocks: eps_t = A0^{-1} %*% u_t
A0_inv <- solve(A0)
eps    <- t(A0_inv %*% t(resid_rf))    # (T-p) x 3
colnames(eps) <- c("eps_gdp", "eps_inf", "eps_mp")

cat("\n=== Structural shock variances (should be ~1) ===\n")
print(round(apply(eps, 2, var), 4))

# Quick plot: MP shocks over time
eps_ts  <- ts(eps, start = tsp(dat_var)[1] + p_sel / 4, frequency = 4)
par(mfrow = c(1, 1))
plot(eps_ts[, "eps_mp"],
     main = "Structural Monetary Policy Shocks (Cholesky SVAR)",
     ylab = "Std. deviations", col = "steelblue", lwd = 1.5)
abline(h = 0, lty = 2, col = "grey50")


#=========================================================================================
# 4.  HISTORICAL DECOMPOSITION VIA MA REPRESENTATION
#
#     y_t  = mu  +  sum_{s=0}^{inf}  Psi_s  %*%  u_{t-s}
#           = mu  +  sum_{s=0}^{inf}  Theta_s  %*%  eps_{t-s}
#
#     where Psi_s   = reduced-form MA matrices (from irf with ortho=FALSE)
#           Theta_s = Psi_s %*% A0  (structural MA matrices)
#
#     Contribution of MP shock (column 3) to variable i at time t:
#       HD_mp[t, i] = sum_{s=0}^{t-1}  Theta_s[i, 3]  *  eps[t-s, 3]
#=========================================================================================

cat("\nComputing MA representation for historical decomposition ...\n")

# Get reduced-form IRF matrices (Psi_s), horizon 0 to T_resid-1
irf_rf <- irf(bv_est,
              n.ahead  = T_resid - 1,
              boot     = FALSE,
              ortho    = FALSE)    # ortho=FALSE => Psi_s matrices

k         <- ncol(dat_var)
var_names <- colnames(dat_var)

# Build Psi array: Psi[response, impulse, horizon+1]
Psi <- array(0, dim = c(k, k, T_resid))
for (j in 1:k) {
  for (s in 0:(T_resid - 1)) {
    Psi[, j, s + 1] <- irf_rf$irf[[var_names[j]]][s + 1, ]
  }
}

# Structural MA: Theta_s = Psi_s %*% A0
Theta <- array(0, dim = dim(Psi))
for (s in 0:(T_resid - 1)) {
  Theta[, , s + 1] <- Psi[, , s + 1] %*% A0
}

# Historical decomposition: contribution of MP shock (j = 3)
HD_mp <- matrix(0, nrow = T_resid, ncol = k)
colnames(HD_mp) <- var_names

for (t in 1:T_resid) {
  for (s in 0:(t - 1)) {
    shock_t <- t - s                   # index into eps (1-based)
    if (shock_t >= 1 && shock_t <= T_resid) {
      HD_mp[t, ] <- HD_mp[t, ] + Theta[, 3, s + 1] * eps[shock_t, 3]
    }
  }
}

cat("Historical decomposition complete.\n")


#=========================================================================================
# 5.  COUNTERFACTUAL PATHS
#     Counterfactual = Actual path  -  Cumulative MP shock contribution
#     Interpretation: the path that inflation / output would have followed
#     if all monetary policy shocks had been zero throughout history.
#=========================================================================================

# Align actual data with residual time span (drop first p_sel obs)
actual   <- as.matrix(dat_var)[(p_sel + 1):nT, ]
cf_path  <- actual - HD_mp

# Time index (numeric fractional years for ggplot)
time_idx <- as.numeric(time(dat_var))[(p_sel + 1):nT]

# Build long-format data frame for ggplot
cf_df <- data.frame(
  time        = time_idx,
  Actual_inf  = actual[, "inf"],
  CF_inf      = cf_path[, "inf"],
  Actual_gdp  = actual[, "gdp_gap"],
  CF_gdp      = cf_path[, "gdp_gap"],
  MP_shock    = eps[, "eps_mp"],
  HD_inf      = HD_mp[, "inf"],
  HD_gdp      = HD_mp[, "gdp_gap"]
)


#=========================================================================================
# 6.  PLOT COUNTERFACTUAL RESULTS
#=========================================================================================

#-- Panel A: Inflation --
pA <- ggplot(cf_df, aes(x = time)) +
  geom_ribbon(aes(ymin = pmin(Actual_inf, CF_inf),
                  ymax = pmax(Actual_inf, CF_inf)),
              fill = "steelblue", alpha = 0.20) +
  geom_line(aes(y = Actual_inf, colour = "Actual"),
            linewidth = 0.85) +
  geom_line(aes(y = CF_inf,
                colour = "Counterfactual (No MP Shocks)"),
            linewidth = 0.85, linetype = "dashed") +
  scale_colour_manual(
    values = c("Actual" = "black",
               "Counterfactual (No MP Shocks)" = "steelblue")) +
  labs(title    = "Inflation: Actual vs. Counterfactual",
       subtitle = "What would inflation have been with no monetary policy shocks?",
       x = NULL, y = "Inflation (%, q-o-q)", colour = NULL) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom",
        plot.title = element_text(face = "bold"))

#-- Panel B: Output Gap --
pB <- ggplot(cf_df, aes(x = time)) +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "grey50") +
  geom_ribbon(aes(ymin = pmin(Actual_gdp, CF_gdp),
                  ymax = pmax(Actual_gdp, CF_gdp)),
              fill = "tomato", alpha = 0.20) +
  geom_line(aes(y = Actual_gdp, colour = "Actual"),
            linewidth = 0.85) +
  geom_line(aes(y = CF_gdp,
                colour = "Counterfactual (No MP Shocks)"),
            linewidth = 0.85, linetype = "dashed") +
  scale_colour_manual(
    values = c("Actual" = "black",
               "Counterfactual (No MP Shocks)" = "tomato")) +
  labs(title    = "Output Gap: Actual vs. Counterfactual",
       subtitle = "What would the output gap have been with no monetary policy shocks?",
       x = NULL, y = "Output Gap (HP-filtered log GDP)", colour = NULL) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom",
        plot.title = element_text(face = "bold"))

#-- Panel C: MP Shock Contribution to Inflation --
pC <- ggplot(cf_df, aes(x = time)) +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "grey50") +
  geom_col(aes(y = HD_inf), fill = "steelblue", alpha = 0.75) +
  labs(title    = "Contribution of Monetary Policy Shocks to Inflation",
       subtitle = "Historical decomposition: MP shock contribution (pp)",
       x = NULL, y = "Contribution (pp)") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

#-- Panel D: Structural MP Shocks --
pD <- ggplot(cf_df, aes(x = time)) +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "grey50") +
  geom_col(aes(y = MP_shock), fill = "darkblue", alpha = 0.65) +
  labs(title    = "Structural Monetary Policy Shocks",
       subtitle = "Cholesky-identified SVAR | Positive = contractionary surprise",
       x = NULL, y = "Shock (std. units)") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

# Display all counterfactual panels
grid.arrange(pA, pB, ncol = 1,
             top = "SVAR Counterfactual Analysis -- South Africa (1981Q4-2012Q4)")
grid.arrange(pC, pD, ncol = 1,
             top = "MP Shock Identification and Historical Decomposition")

# Save
ggsave("CF_inflation.png",    pA, width = 11, height = 5, dpi = 150)
ggsave("CF_output_gap.png",   pB, width = 11, height = 5, dpi = 150)
ggsave("CF_HD_inflation.png", pC, width = 11, height = 4, dpi = 150)
ggsave("CF_MP_shocks.png",    pD, width = 11, height = 4, dpi = 150)

cat("\n=== Counterfactual Summary (last 8 quarters) ===\n")
cf_summary <- data.frame(
  Time             = round(tail(time_idx, 8), 2),
  Actual_Inf       = round(tail(actual[, "inf"], 8), 3),
  CF_Inf_NoMP      = round(tail(cf_path[, "inf"], 8), 3),
  MP_Contribution  = round(tail(HD_mp[, "inf"], 8), 3),
  Actual_GDPgap    = round(tail(actual[, "gdp_gap"], 8), 4),
  CF_GDPgap_NoMP   = round(tail(cf_path[, "gdp_gap"], 8), 4)
)
print(cf_summary)


#=========================================================================================
# 7.  FORECAST SCENARIOS: SUSTAINED MONETARY POLICY TIGHTENING
#
#     Approach: Conditional (scenario) forecast
#       - Fix the interest rate path to reflect a sustained hike
#       - Propagate through the VAR to obtain conditional paths for
#         inflation and output
#
#     Scenarios:
#       S0 - Baseline          : unrestricted VAR point forecast
#       S1 - Moderate tightening  : +50 bp per quarter for 4 quarters,
#                                   then held at the terminal rate
#       S2 - Aggressive tightening: +100 bp per quarter for 4 quarters,
#                                   then held at the terminal rate
#=========================================================================================

h        <- 12     # forecast horizon (3 years)
k_var    <- ncol(dat_var)

# VAR coefficients
B_all    <- Bcoef(bv_est)                      # k x (k*p + 1)
const_v  <- B_all[, ncol(B_all)]               # intercept vector (k x 1)
B_list   <- lapply(seq_len(p_sel), function(j) {
  B_all[, ((j - 1) * k_var + 1):(j * k_var)]  # B_j coefficient matrix
})

# Last p_sel observations as starting conditions (rows = oldest to newest)
last_obs <- tail(as.matrix(dat_var), p_sel)

# Last observed interest rate
last_int <- as.numeric(tail(int, 1))

# ---- Helper: simulate VAR h steps forward ----
# int_constraint: numeric vector of length h for the interest rate path;
#                 NULL means unconstrained (pure VAR forecast)
sim_var_forward <- function(last_obs, h, int_constraint = NULL) {
  p  <- nrow(last_obs)
  fc <- matrix(NA, h, k_var)
  colnames(fc) <- var_names

  get_obs <- function(t_fc, lag) {
    idx <- t_fc - lag
    if (idx <= 0) last_obs[p + idx, ] else fc[idx, ]
  }

  for (t in seq_len(h)) {
    y_new <- const_v
    for (j in seq_len(p)) {
      y_new <- y_new + drop(B_list[[j]] %*% get_obs(t, j))
    }
    # Override interest rate with scenario path if provided
    if (!is.null(int_constraint) && t <= length(int_constraint)) {
      y_new["int"] <- int_constraint[t]
    }
    fc[t, ] <- y_new
  }
  return(fc)
}

# ---- Scenario interest rate paths ----

# Baseline: use predict() to get standard VAR forecast (for CI bands)
base_pred  <- predict(bv_est, n.ahead = h, ci = 0.95)
int_s0     <- base_pred$fcst$int[, "fcst"]       # unconstrained int forecast

# Moderate: +50bp each quarter for 4 quarters, then hold
ramp_mod   <- cumsum(rep(0.50, 4))                # 0.5, 1.0, 1.5, 2.0 pp above last obs
int_s1     <- c(last_int + ramp_mod,
                rep(last_int + ramp_mod[4], h - 4))

# Aggressive: +100bp each quarter for 4 quarters, then hold
ramp_agg   <- cumsum(rep(1.00, 4))                # 1, 2, 3, 4 pp above last obs
int_s2     <- c(last_int + ramp_agg,
                rep(last_int + ramp_agg[4], h - 4))

# ---- Run simulations ----
fc_s0 <- sim_var_forward(last_obs, h, int_constraint = NULL)   # baseline (VAR)
fc_s1 <- sim_var_forward(last_obs, h, int_constraint = int_s1) # moderate
fc_s2 <- sim_var_forward(last_obs, h, int_constraint = int_s2) # aggressive

# 95% CI from baseline predict() object
inf_lo  <- base_pred$fcst$inf[, "lower"]
inf_hi  <- base_pred$fcst$inf[, "upper"]
gdp_lo  <- base_pred$fcst$gdp_gap[, "lower"]
gdp_hi  <- base_pred$fcst$gdp_gap[, "upper"]

# Forecast time index
last_t   <- max(as.numeric(time(dat_var)))
fc_time  <- seq(last_t + 0.25, by = 0.25, length.out = h)

fc_df <- data.frame(
  time          = fc_time,
  int_baseline  = fc_s0[, "int"],
  int_moderate  = int_s1,
  int_aggressive = int_s2,
  inf_baseline  = fc_s0[, "inf"],
  inf_moderate  = fc_s1[, "inf"],
  inf_aggressive = fc_s2[, "inf"],
  inf_lo        = inf_lo,
  inf_hi        = inf_hi,
  gdp_baseline  = fc_s0[, "gdp_gap"],
  gdp_moderate  = fc_s1[, "gdp_gap"],
  gdp_aggressive = fc_s2[, "gdp_gap"],
  gdp_lo        = gdp_lo,
  gdp_hi        = gdp_hi
)


#=========================================================================================
# 8.  PLOT FORECAST SCENARIOS
#=========================================================================================

scen_colours <- c(
  "Baseline"               = "black",
    "Moderate (+50bp x 4Q)"  = "steelblue",
      "Aggressive (+100bp x 4Q)" = "tomato"
      )

      #-- Interest rate path --
      pE <- ggplot(fc_df, aes(x = time)) +
        geom_line(aes(y = int_baseline,   colour = "Baseline"),
                    linewidth = 0.9) +
                      geom_line(aes(y = int_moderate,   colour = "Moderate (+50bp x 4Q)"),
                                  linewidth = 0.9, linetype = "dashed") +
                                    geom_line(aes(y = int_aggressive, colour = "Aggressive (+100bp x 4Q)"),
                                                linewidth = 0.9, linetype = "dotdash") +
                                                  scale_colour_manual(values = scen_colours) +
                                                    labs(title    = "Assumed Interest Rate Path Under Policy Scenarios",
                                                           subtitle = paste("Forecast horizon:", h, "quarters | Starting rate:",
                                                                                   round(last_int, 2), "%"),
                                                                                          x = NULL, y = "Interest Rate (%)", colour = "Scenario") +
                                                                                            theme_bw(base_size = 12) +
                                                                                              theme(legend.position = "bottom",
                                                                                                      plot.title = element_text(face = "bold"))

                                                                                                      #-- Inflation forecast --
                                                                                                      pF <- ggplot(fc_df, aes(x = time)) +
                                                                                                        geom_ribbon(aes(ymin = inf_lo, ymax = inf_hi),
                                                                                                                      fill = "grey70", alpha = 0.35) +
                                                                                                                        geom_line(aes(y = inf_baseline,   colour = "Baseline"),
                                                                                                                                    linewidth = 0.9) +
                                                                                                                                      geom_line(aes(y = inf_moderate,   colour = "Moderate (+50bp x 4Q)"),
                                                                                                                                                  linewidth = 0.9, linetype = "dashed") +
                                                                                                                                                    geom_line(aes(y = inf_aggressive, colour = "Aggressive (+100bp x 4Q)"),
                                                                                                                                                                linewidth = 0.9, linetype = "dotdash") +
                                                                                                                                                                  scale_colour_manual(values = scen_colours) +
                                                                                                                                                                    labs(title    = "Inflation Forecast Under Monetary Policy Tightening Scenarios",
                                                                                                                                                                           subtitle = "Shaded band = 95% CI for baseline | VAR conditional forecast",
                                                                                                                                                                                  x = NULL, y = "Inflation (%, q-o-q)", colour = "Scenario") +
                                                                                                                                                                                    theme_bw(base_size = 12) +
                                                                                                                                                                                      theme(legend.position = "bottom",
                                                                                                                                                                                              plot.title = element_text(face = "bold"))

                                                                                                                                                                                              #-- Output gap forecast --
                                                                                                                                                                                              pG <- ggplot(fc_df, aes(x = time)) +
                                                                                                                                                                                                geom_hline(yintercept = 0, linetype = "dotted", colour = "grey50") +
                                                                                                                                                                                                  geom_ribbon(aes(ymin = gdp_lo, ymax = gdp_hi),
                                                                                                                                                                                                                fill = "grey70", alpha = 0.35) +
                                                                                                                                                                                                                  geom_line(aes(y = gdp_baseline,   colour = "Baseline"),
                                                                                                                                                                                                                              linewidth = 0.9) +
                                                                                                                                                                                                                                geom_line(aes(y = gdp_moderate,   colour = "Moderate (+50bp x 4Q)"),
                                                                                                                                                                                                                                            linewidth = 0.9, linetype = "dashed") +
                                                                                                                                                                                                                                              geom_line(aes(y = gdp_aggressive, colour = "Aggressive (+100bp x 4Q)"),
                                                                                                                                                                                                                                                          linewidth = 0.9, linetype = "dotdash") +
                                                                                                                                                                                                                                                            scale_colour_manual(values = scen_colours) +
                                                                                                                                                                                                                                                              labs(title    = "Output Gap Forecast Under Monetary Policy Tightening Scenarios",
                                                                                                                                                                                                                                                                     subtitle = "Shaded band = 95% CI for baseline",
                                                                                                                                                                                                                                                                            x = NULL, y = "Output Gap (log deviation)", colour = "Scenario") +
                                                                                                                                                                                                                                                                              theme_bw(base_size = 12) +
                                                                                                                                                                                                                                                                                theme(legend.position = "bottom",
                                                                                                                                                                                                                                                                                        plot.title = element_text(face = "bold"))

                                                                                                                                                                                                                                                                                        # Display
                                                                                                                                                                                                                                                                                        grid.arrange(pE, pF, pG, ncol = 1,
                                                                                                                                                                                                                                                                                                     top = "VAR Conditional Forecast: Monetary Policy Tightening Scenarios")

                                                                                                                                                                                                                                                                                                     ggsave("Scenario_IntRate_Path.png",         pE, width = 11, height = 4, dpi = 150)
                                                                                                                                                                                                                                                                                                     ggsave("Scenario_Inflation_Forecast.png",   pF, width = 11, height = 5, dpi = 150)
                                                                                                                                                                                                                                                                                                     ggsave("Scenario_OutputGap_Forecast.png",   pG, width = 11, height = 5, dpi = 150)
                                                                                                                                                                                                                                                                                                     

#=========================================================================================
# 9.  COMBINED CHART: Last 20 Quarters of History + Forecast
#=========================================================================================

hist_w    <- 20
hist_time <- tail(as.numeric(time(inf)), hist_w)
hist_inf  <- tail(as.numeric(inf), hist_w)
hist_gdp  <- tail(as.numeric(gdp_gap), hist_w)

combo_df <- bind_rows(
  data.frame(time    = hist_time,
             segment = "Historical",
             inf_val = hist_inf,
             gdp_val = hist_gdp,
             series  = "Actual"),
  data.frame(time    = fc_time,
             segment = "Forecast",
             inf_val = fc_s0[, "inf"],
             gdp_val = fc_s0[, "gdp_gap"],
             series  = "Baseline"),
  data.frame(time    = fc_time,
             segment = "Forecast",
             inf_val = fc_s1[, "inf"],
             gdp_val = fc_s1[, "gdp_gap"],
             series  = "Moderate (+50bp x 4Q)"),
  data.frame(time    = fc_time,
             segment = "Forecast",
             inf_val = fc_s2[, "inf"],
             gdp_val = fc_s2[, "gdp_gap"],
             series  = "Aggressive (+100bp x 4Q)")
)

# Ribbon only for baseline CI
ribbon_df <- data.frame(time = fc_time, inf_lo = inf_lo, inf_hi = inf_hi,
                        gdp_lo = gdp_lo, gdp_hi = gdp_hi)

pH <- ggplot() +
  geom_vline(xintercept = last_t, linetype = "dashed", colour = "grey40") +
  annotate("text", x = last_t + 0.1, y = Inf, label = "Forecast start",
           vjust = 1.5, hjust = 0, size = 3.2, colour = "grey40") +
  # Historical (black)
  geom_line(data = combo_df[combo_df$segment == "Historical", ],
            aes(x = time, y = inf_val), colour = "black", linewidth = 0.9) +
  # Baseline CI ribbon
  geom_ribbon(data = ribbon_df,
              aes(x = time, ymin = inf_lo, ymax = inf_hi),
              fill = "grey70", alpha = 0.35) +
  # Forecast scenarios
  geom_line(data = combo_df[combo_df$segment == "Forecast", ],
            aes(x = time, y = inf_val, colour = series, linetype = series),
            linewidth = 0.9) +
  scale_colour_manual(values = c(
    "Baseline"                = "black",
    "Moderate (+50bp x 4Q)"  = "steelblue",
    "Aggressive (+100bp x 4Q)" = "tomato")) +
  scale_linetype_manual(values = c(
    "Baseline"                = "solid",
    "Moderate (+50bp x 4Q)"  = "dashed",
    "Aggressive (+100bp x 4Q)" = "dotdash")) +
  labs(title    = "Inflation: Historical Path and Policy Scenario Forecasts",
       subtitle = "South Africa | VAR-Based Conditional Forecast | 95% CI (baseline)",
       x = NULL, y = "Inflation (%, q-o-q)",
       colour = "Series", linetype = "Series") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom",
        plot.title = element_text(face = "bold"))

print(pH)
ggsave("Combined_Inflation_History_Forecast.png", pH, width = 13, height = 5.5, dpi = 150)


#=========================================================================================
# 10.  SUMMARY TABLES
#=========================================================================================

cat("\n")
cat("================================================================\n")
cat("  COUNTERFACTUAL: Effect of Removing All MP Shocks\n")
cat("  (Average deviation over full sample)\n")
cat("================================================================\n")
cat(sprintf("  Inflation avg. MP contribution : %+.3f pp\n",
            mean(HD_mp[, "inf"])))
cat(sprintf("  Output gap avg. MP contribution: %+.4f\n",
            mean(HD_mp[, "gdp_gap"])))
cat(sprintf("  Inflation std dev -- Actual  : %.3f\n",
            sd(actual[, "inf"])))
cat(sprintf("  Inflation std dev -- CF (no MP): %.3f\n",
            sd(cf_path[, "inf"])))

cat("\n")
cat("================================================================\n")
cat("  FORECAST SCENARIO SUMMARY (", h, "quarters ahead)\n")
cat("================================================================\n")

fc_table <- data.frame(
  Quarter          = paste0("Q+", seq_len(h)),
  Int_Baseline     = round(fc_s0[, "int"], 3),
  Int_Moderate     = round(int_s1, 3),
  Int_Aggressive   = round(int_s2, 3),
  Inf_Baseline     = round(fc_s0[, "inf"], 3),
  Inf_Moderate     = round(fc_s1[, "inf"], 3),
  Inf_Aggressive   = round(fc_s2[, "inf"], 3),
  GDP_Baseline     = round(fc_s0[, "gdp_gap"], 4),
  GDP_Moderate     = round(fc_s1[, "gdp_gap"], 4),
  GDP_Aggressive   = round(fc_s2[, "gdp_gap"], 4)
)

print(fc_table, row.names = FALSE)
write.csv(fc_table, "Forecast_Scenario_Summary.csv", row.names = FALSE)

cat("\n")
cat("================================================================\n")
cat("  INFLATION IMPACT OF TIGHTENING vs. BASELINE\n")
cat("  (pp difference from baseline, averaged over forecast)\n")
cat("================================================================\n")
cat(sprintf("  Moderate tightening  : %+.3f pp\n",
            mean(fc_s1[, "inf"] - fc_s0[, "inf"])))
cat(sprintf("  Aggressive tightening: %+.3f pp\n",
            mean(fc_s2[, "inf"] - fc_s0[, "inf"])))

cat("\n")
cat("================================================================\n")
cat("  OUTPUT GAP IMPACT OF TIGHTENING vs. BASELINE\n")
cat("================================================================\n")
cat(sprintf("  Moderate tightening  : %+.4f\n",
            mean(fc_s1[, "gdp_gap"] - fc_s0[, "gdp_gap"])))
cat(sprintf("  Aggressive tightening: %+.4f\n",
            mean(fc_s2[, "gdp_gap"] - fc_s0[, "gdp_gap"])))

cat("\nAll output files saved in working directory.\n")
cat("=== DONE ===\n")
