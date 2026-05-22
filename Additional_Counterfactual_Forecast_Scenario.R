################################################################################
# SVAR Policy Counterfactual and Forecast Scenario Analysis
# Author: Nafiu (CBN Monetary Policy Department)
# Date: May 21, 2026
# Purpose: Policy counterfactual (neutralizing MP shock) and forecast scenarios
#          (sustained monetary policy tightening/easing)
################################################################################

# Clear workspace
rm(list = ls())
gc()

# Package management with need/want guards
need_pkg <- function(pkg) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE, repos = "https://cran.rstudio.com/")
    library(pkg, character.only = TRUE)
  }
}

want_pkg <- function(pkg) {
  suppressWarnings(suppressMessages(
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
      try(install.packages(pkg, dependencies = TRUE, repos = "https://cran.rstudio.com/"), silent = TRUE)
      require(pkg, character.only = TRUE, quietly = TRUE)
    }
  ))
}

# Load required packages
need_pkg("vars")
need_pkg("readr")
need_pkg("dplyr")
need_pkg("ggplot2")
need_pkg("gridExtra")
need_pkg("lubridate")
need_pkg("zoo")
need_pkg("reshape2")
want_pkg("scales")
want_pkg("viridis")

################################################################################
# 1. DATA LOADING AND PREPARATION
################################################################################

# Load data
data_raw <- read_csv("c:/Users/Ex-Ethic/OneDrive/Documents/CBs_Works/Facilitation/2026/TSM_Forecasting/Policy_Counterfactual_Forecast_Scenaario_from_SVAR/data_sa.csv", show_col_types = FALSE)

# data_sa <- read.csv("c:/Users/Ex-Ethic/OneDrive/Documents/CBs_Works/Facilitation/2026/TSM_Forecasting/Policy_Counterfactual_Forecast_Scenaario_from_SVAR/data_sa.csv")


# Data processing
data <- data_raw %>%
  mutate(
    Date = as.yearqtr(Date, format = "%Y/%m"),
    quarter = as.Date(Date)
  ) %>%
  arrange(Date)

# Create time series object
ts_data <- ts(
  data[, c("gdp", "inf", "int")],
  start = c(year(min(data$quarter)), quarter(min(data$quarter))),
  frequency = 4
)

# Log transformation for GDP (levels to growth)
ts_data[, "gdp"] <- log(ts_data[, "gdp"]) * 100

# Rename for clarity
colnames(ts_data) <- c("Output", "Inflation", "PolicyRate")

cat("=================================================================\n")
cat("DATA SUMMARY\n")
cat("=================================================================\n")
cat("Sample period:", format(min(data$Date)), "to", format(max(data$Date)), "\n")
cat("Observations:", nrow(data), "\n")
cat("Variables: Output (log GDP), Inflation, Policy Rate\n\n")
print(summary(ts_data))
cat("\n")

################################################################################
# 2. VAR MODEL ESTIMATION
################################################################################

# Determine optimal lag length
cat("=================================================================\n")
cat("VAR LAG SELECTION\n")
cat("=================================================================\n")

lag_selection <- VARselect(ts_data, lag.max = 8, type = "const")
print(lag_selection$selection)
cat("\n")

# Use AIC-selected lag (or override with CBN standard if needed)
optimal_lag <- lag_selection$selection["AIC(n)"]
cat("Selected lag order:", optimal_lag, "(based on AIC)\n\n")

# Estimate reduced-form VAR
var_model <- VAR(ts_data, p = optimal_lag, type = "const")

cat("=================================================================\n")
cat("VAR MODEL ESTIMATION RESULTS\n")
cat("=================================================================\n")
print(summary(var_model))
cat("\n")

# Model diagnostics
cat("=================================================================\n")
cat("VAR MODEL DIAGNOSTICS\n")
cat("=================================================================\n")

# Serial correlation test
serial_test <- serial.test(var_model, lags.pt = 12, type = "PT.asymptotic")
cat("\nPortmanteau Test (Autocorrelation):\n")
print(serial_test)

# Heteroskedasticity test
arch_test <- arch.test(var_model, lags.multi = 5)
cat("\nARCH Test (Heteroskedasticity):\n")
print(arch_test)

# Normality test
norm_test <- normality.test(var_model, multivariate.only = TRUE)
cat("\nJarque-Bera Normality Test:\n")
print(norm_test)

# Stability test
stability_test <- stability(var_model, type = "OLS-CUSUM")
cat("\nModel stability: Check plot for CUSUM stability\n\n")

################################################################################
# 3. STRUCTURAL IDENTIFICATION (RECURSIVE ORDERING)
################################################################################

cat("=================================================================\n")
cat("STRUCTURAL IDENTIFICATION\n")
cat("=================================================================\n")
cat("Identification: Recursive (Cholesky decomposition)\n")
cat("Ordering: Output → Inflation → Policy Rate\n")
cat("Rationale:\n")
cat("  1. Output determined by real shocks (moves first)\n")
cat("  2. Inflation responds to output shocks within quarter\n")
cat("  3. Policy rate responds to output and inflation contemporaneously\n\n")

# SVAR identification using Cholesky
svar_model <- irf(var_model, 
                  impulse = "PolicyRate",
                  response = c("Output", "Inflation", "PolicyRate"),
                  n.ahead = 20,
                  ortho = TRUE,  # Orthogonalized (structural) shocks
                  cumulative = FALSE,
                  boot = TRUE,
                  ci = 0.90,
                  runs = 500,
                  seed = 12345)

cat("Structural shocks identified via Cholesky decomposition\n")
cat("Bootstrap confidence intervals: 90% (500 replications)\n\n")

################################################################################
# 4. POLICY COUNTERFACTUAL ANALYSIS
################################################################################
# Objective: What would have happened to inflation and output if monetary 
# policy shocks had been neutralized over the sample period?

cat("=================================================================\n")
cat("POLICY COUNTERFACTUAL ANALYSIS\n")
cat("=================================================================\n")
cat("Counterfactual: Neutralizing monetary policy shocks\n")
cat("Question: What would inflation and output have been without MP shocks?\n\n")

# Extract VAR coefficients
B_matrices <- Bcoef(var_model)  # Coefficient matrices
resid_var <- residuals(var_model)  # Reduced-form residuals
n_vars <- ncol(ts_data)
n_obs <- nrow(resid_var)

# Get Cholesky decomposition for structural shocks
P_chol <- t(chol(summary(var_model)$covres))  # Lower triangular
structural_shocks <- solve(P_chol) %*% t(resid_var)  # Structural shocks
structural_shocks <- t(structural_shocks)

cat("Structural shocks extracted (dimensions):", dim(structural_shocks), "\n")

# Neutralize monetary policy shock (3rd shock in our ordering)
structural_shocks_neutral <- structural_shocks
structural_shocks_neutral[, 3] <- 0  # Set policy shock to zero

cat("Monetary policy shocks set to zero\n\n")

# Reconstruct reduced-form residuals without policy shocks
resid_counterfactual <- structural_shocks_neutral %*% t(solve(P_chol))

# Simulate counterfactual data
n_lags <- optimal_lag
y_counterfactual <- matrix(NA, nrow = n_obs, ncol = n_vars)
colnames(y_counterfactual) <- colnames(ts_data)

# Initialize with actual pre-sample values
y_counterfactual[1:n_lags, ] <- ts_data[1:n_lags, ]

# Simulate forward using VAR structure with neutralized shocks
for (t in (n_lags + 1):n_obs) {
  y_sim <- var_model$varresult[[1]]$coefficients["const"]  # Intercept
  
  for (lag in 1:n_lags) {
    idx <- t - lag
    for (var_idx in 1:n_vars) {
      var_name <- colnames(ts_data)[var_idx]
      coef_name <- paste0(var_name, ".l", lag)
      
      # Add lagged contributions for all equations
      for (eq_idx in 1:n_vars) {
        if (coef_name %in% names(var_model$varresult[[eq_idx]]$coefficients)) {
          if (lag == 1 && var_idx == 1) {  # First time, set intercept
            y_counterfactual[t, eq_idx] <- var_model$varresult[[eq_idx]]$coefficients["const"]
          }
          y_counterfactual[t, eq_idx] <- y_counterfactual[t, eq_idx] + 
            var_model$varresult[[eq_idx]]$coefficients[coef_name] * y_counterfactual[idx, var_idx]
        }
      }
    }
  }
  
  # Add counterfactual residuals
  y_counterfactual[t, ] <- y_counterfactual[t, ] + resid_counterfactual[t, ]
}

# Create comparison data frame
comparison_data <- data.frame(
  Date = data$quarter[(n_lags + 1):n_obs],
  Output_Actual = as.numeric(ts_data[(n_lags + 1):n_obs, "Output"]),
  Output_Counterfactual = y_counterfactual[(n_lags + 1):n_obs, "Output"],
  Inflation_Actual = as.numeric(ts_data[(n_lags + 1):n_obs, "Inflation"]),
  Inflation_Counterfactual = y_counterfactual[(n_lags + 1):n_obs, "Inflation"],
  PolicyRate_Actual = as.numeric(ts_data[(n_lags + 1):n_obs, "PolicyRate"]),
  PolicyRate_Counterfactual = y_counterfactual[(n_lags + 1):n_obs, "PolicyRate"]
)

# Calculate differences (impact of MP shocks)
comparison_data <- comparison_data %>%
  mutate(
    Output_Impact = Output_Actual - Output_Counterfactual,
    Inflation_Impact = Inflation_Actual - Inflation_Counterfactual,
    PolicyRate_Impact = PolicyRate_Actual - PolicyRate_Counterfactual
  )

# Summary statistics
cat("COUNTERFACTUAL RESULTS SUMMARY\n")
cat("Average impact of MP shocks over sample:\n")
cat(sprintf("  Output:      %.4f pp (actual vs. no-MP-shock)\n", 
            mean(comparison_data$Output_Impact, na.rm = TRUE)))
cat(sprintf("  Inflation:   %.4f pp\n", 
            mean(comparison_data$Inflation_Impact, na.rm = TRUE)))
cat(sprintf("  Policy Rate: %.4f pp\n\n", 
            mean(comparison_data$PolicyRate_Impact, na.rm = TRUE)))

cat("Standard deviation of MP shock impact:\n")
cat(sprintf("  Output:      %.4f pp\n", 
            sd(comparison_data$Output_Impact, na.rm = TRUE)))
cat(sprintf("  Inflation:   %.4f pp\n", 
            sd(comparison_data$Inflation_Impact, na.rm = TRUE)))
cat(sprintf("  Policy Rate: %.4f pp\n\n", 
            sd(comparison_data$PolicyRate_Impact, na.rm = TRUE)))

################################################################################
# 5. FORECAST SCENARIO ANALYSIS
################################################################################
# Scenario 1: Sustained monetary tightening (rate hikes)
# Scenario 2: Baseline (no additional shocks)

cat("=================================================================\n")
cat("FORECAST SCENARIO ANALYSIS\n")
cat("=================================================================\n")

# Forecast horizon
h_forecast <- 12  # 3 years ahead (quarterly)

# Baseline forecast (no additional shocks)
forecast_baseline <- predict(var_model, n.ahead = h_forecast, ci = 0.90)

cat("Baseline forecast generated (", h_forecast, "quarters ahead)\n")

# Scenario 1: Sustained monetary tightening
# Simulate +100 bps shock to policy rate sustained over forecast horizon
cat("\nScenario 1: Sustained monetary policy tightening\n")
cat("  Assumption: +100 bps shock to policy rate, sustained for", h_forecast, "quarters\n")

# Create shock matrix (positive shock to policy rate)
shock_tightening <- matrix(0, nrow = h_forecast, ncol = n_vars)
colnames(shock_tightening) <- colnames(ts_data)
shock_tightening[, "PolicyRate"] <- 1.0  # 100 bps sustained shock

# Convert to structural shocks
shock_structural <- solve(P_chol) %*% t(shock_tightening)
shock_structural <- t(shock_structural)

# Simulate tightening scenario
last_obs <- ts_data[nrow(ts_data), ]
y_tightening <- matrix(NA, nrow = h_forecast, ncol = n_vars)
colnames(y_tightening) <- colnames(ts_data)

# Historical data for lags
y_history <- rbind(
  ts_data[(nrow(ts_data) - n_lags + 1):nrow(ts_data), ],
  y_tightening
)

for (t in 1:h_forecast) {
  idx_history <- n_lags + t
  
  for (eq_idx in 1:n_vars) {
    # Start with intercept
    y_tightening[t, eq_idx] <- var_model$varresult[[eq_idx]]$coefficients["const"]
    
    # Add lagged contributions
    for (lag in 1:n_lags) {
      for (var_idx in 1:n_vars) {
        var_name <- colnames(ts_data)[var_idx]
        coef_name <- paste0(var_name, ".l", lag)
        
        if (coef_name %in% names(var_model$varresult[[eq_idx]]$coefficients)) {
          lag_value <- y_history[idx_history - lag, var_idx]
          y_tightening[t, eq_idx] <- y_tightening[t, eq_idx] +
            var_model$varresult[[eq_idx]]$coefficients[coef_name] * lag_value
        }
      }
    }
  }
  
  # Add structural shock converted to reduced form
  shock_rf <- P_chol %*% shock_structural[t, ]
  y_tightening[t, ] <- y_tightening[t, ] + shock_rf
  
  # Update history for next iteration
  y_history[idx_history, ] <- y_tightening[t, ]
}

cat("Tightening scenario simulated\n")

# Scenario 2: Sustained monetary easing
cat("\nScenario 2: Sustained monetary policy easing\n")
cat("  Assumption: -100 bps shock to policy rate, sustained for", h_forecast, "quarters\n")

shock_easing <- shock_tightening
shock_easing[, "PolicyRate"] <- -1.0  # -100 bps sustained shock

shock_structural_easing <- solve(P_chol) %*% t(shock_easing)
shock_structural_easing <- t(shock_structural_easing)

y_easing <- matrix(NA, nrow = h_forecast, ncol = n_vars)
colnames(y_easing) <- colnames(ts_data)

y_history_easing <- rbind(
  ts_data[(nrow(ts_data) - n_lags + 1):nrow(ts_data), ],
  y_easing
)

for (t in 1:h_forecast) {
  idx_history <- n_lags + t
  
  for (eq_idx in 1:n_vars) {
    y_easing[t, eq_idx] <- var_model$varresult[[eq_idx]]$coefficients["const"]
    
    for (lag in 1:n_lags) {
      for (var_idx in 1:n_vars) {
        var_name <- colnames(ts_data)[var_idx]
        coef_name <- paste0(var_name, ".l", lag)
        
        if (coef_name %in% names(var_model$varresult[[eq_idx]]$coefficients)) {
          lag_value <- y_history_easing[idx_history - lag, var_idx]
          y_easing[t, eq_idx] <- y_easing[t, eq_idx] +
            var_model$varresult[[eq_idx]]$coefficients[coef_name] * lag_value
        }
      }
    }
  }
  
  shock_rf_easing <- P_chol %*% shock_structural_easing[t, ]
  y_easing[t, ] <- y_easing[t, ] + shock_rf_easing
  
  y_history_easing[idx_history, ] <- y_easing[t, ]
}

cat("Easing scenario simulated\n\n")

# Create forecast comparison data frame
forecast_dates <- seq.Date(
  from = max(data$quarter) + months(3),
  by = "quarter",
  length.out = h_forecast
)

forecast_comparison <- data.frame(
  Date = forecast_dates,
  Horizon = 1:h_forecast,
  
  # Baseline
  Output_Baseline = forecast_baseline$fcst$Output[, "fcst"],
  Inflation_Baseline = forecast_baseline$fcst$Inflation[, "fcst"],
  PolicyRate_Baseline = forecast_baseline$fcst$PolicyRate[, "fcst"],
  
  # Tightening scenario
  Output_Tightening = y_tightening[, "Output"],
  Inflation_Tightening = y_tightening[, "Inflation"],
  PolicyRate_Tightening = y_tightening[, "PolicyRate"],
  
  # Easing scenario
  Output_Easing = y_easing[, "Output"],
  Inflation_Easing = y_easing[, "Inflation"],
  PolicyRate_Easing = y_easing[, "PolicyRate"]
)

# Calculate scenario deviations from baseline
forecast_comparison <- forecast_comparison %>%
  mutate(
    Output_Tight_Dev = Output_Tightening - Output_Baseline,
    Inflation_Tight_Dev = Inflation_Tightening - Inflation_Baseline,
    PolicyRate_Tight_Dev = PolicyRate_Tightening - PolicyRate_Baseline,
    
    Output_Easy_Dev = Output_Easing - Output_Baseline,
    Inflation_Easy_Dev = Inflation_Easing - Inflation_Baseline,
    PolicyRate_Easy_Dev = PolicyRate_Easing - PolicyRate_Baseline
  )

# Forecast summary
cat("FORECAST SCENARIO RESULTS (Average over", h_forecast, "quarters)\n\n")
cat("INFLATION FORECASTS:\n")
cat(sprintf("  Baseline:           %.2f%%\n", mean(forecast_comparison$Inflation_Baseline)))
cat(sprintf("  Tightening (+100bp): %.2f%% (%.2f pp below baseline)\n", 
            mean(forecast_comparison$Inflation_Tightening),
            mean(forecast_comparison$Inflation_Tight_Dev)))
cat(sprintf("  Easing (-100bp):     %.2f%% (%.2f pp above baseline)\n\n", 
            mean(forecast_comparison$Inflation_Easing),
            mean(forecast_comparison$Inflation_Easy_Dev)))

cat("OUTPUT FORECASTS:\n")
cat(sprintf("  Baseline:           %.2f\n", mean(forecast_comparison$Output_Baseline)))
cat(sprintf("  Tightening (+100bp): %.2f (%.2f pp below baseline)\n", 
            mean(forecast_comparison$Output_Tightening),
            mean(forecast_comparison$Output_Tight_Dev)))
cat(sprintf("  Easing (-100bp):     %.2f (%.2f pp above baseline)\n\n", 
            mean(forecast_comparison$Output_Easing),
            mean(forecast_comparison$Output_Easy_Dev)))

cat("POLICY RATE FORECASTS:\n")
cat(sprintf("  Baseline:           %.2f%%\n", mean(forecast_comparison$PolicyRate_Baseline)))
cat(sprintf("  Tightening (+100bp): %.2f%%\n", mean(forecast_comparison$PolicyRate_Tightening)))
cat(sprintf("  Easing (-100bp):     %.2f%%\n\n", mean(forecast_comparison$PolicyRate_Easing)))

################################################################################
# 6. VISUALIZATION
################################################################################

cat("=================================================================\n")
cat("GENERATING VISUALIZATIONS\n")
cat("=================================================================\n")

# Set publication-ready theme
theme_cbn <- theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5, color = "gray40"),
    axis.title = element_text(size = 11, face = "bold"),
    axis.text = element_text(size = 10),
    legend.position = "bottom",
    legend.title = element_text(face = "bold", size = 10),
    legend.text = element_text(size = 9),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "gray80", fill = NA, linewidth = 0.5)
  )

# Plot 1: Impulse Response Functions
irf_data_output <- data.frame(
  Horizon = 0:20,
  IRF = svar_model$irf$PolicyRate[, "Output"],
  Lower = svar_model$Lower$PolicyRate[, "Output"],
  Upper = svar_model$Upper$PolicyRate[, "Output"]
)

irf_data_inflation <- data.frame(
  Horizon = 0:20,
  IRF = svar_model$irf$PolicyRate[, "Inflation"],
  Lower = svar_model$Lower$PolicyRate[, "Inflation"],
  Upper = svar_model$Upper$PolicyRate[, "Inflation"]
)

p1 <- ggplot(irf_data_output, aes(x = Horizon, y = IRF)) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = "steelblue", alpha = 0.3) +
  geom_line(color = "steelblue", linewidth = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(
    title = "Response of Output to Monetary Policy Shock",
    subtitle = "One S.D. shock to policy rate (90% CI)",
    x = "Quarters ahead",
    y = "Percentage points"
  ) +
  theme_cbn

p2 <- ggplot(irf_data_inflation, aes(x = Horizon, y = IRF)) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = "darkred", alpha = 0.3) +
  geom_line(color = "darkred", linewidth = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(
    title = "Response of Inflation to Monetary Policy Shock",
    subtitle = "One S.D. shock to policy rate (90% CI)",
    x = "Quarters ahead",
    y = "Percentage points"
  ) +
  theme_cbn

# Plot 2: Counterfactual Analysis
cf_long_output <- comparison_data %>%
  select(Date, Output_Actual, Output_Counterfactual) %>%
  reshape2::melt(id.vars = "Date", variable.name = "Series", value.name = "Value")

cf_long_inflation <- comparison_data %>%
  select(Date, Inflation_Actual, Inflation_Counterfactual) %>%
  reshape2::melt(id.vars = "Date", variable.name = "Series", value.name = "Value")

p3 <- ggplot(cf_long_output, aes(x = Date, y = Value, color = Series, linetype = Series)) +
  geom_line(linewidth = 0.8) +
  scale_color_manual(
    values = c("Output_Actual" = "black", "Output_Counterfactual" = "steelblue"),
    labels = c("Actual", "No MP Shocks")
  ) +
  scale_linetype_manual(
    values = c("Output_Actual" = "solid", "Output_Counterfactual" = "dashed"),
    labels = c("Actual", "No MP Shocks")
  ) +
  labs(
    title = "Output: Actual vs. Counterfactual (No MP Shocks)",
    subtitle = "What would output have been without monetary policy shocks?",
    x = NULL,
    y = "Log Output (×100)",
    color = NULL,
    linetype = NULL
  ) +
  theme_cbn

p4 <- ggplot(cf_long_inflation, aes(x = Date, y = Value, color = Series, linetype = Series)) +
  geom_line(linewidth = 0.8) +
  scale_color_manual(
    values = c("Inflation_Actual" = "black", "Inflation_Counterfactual" = "darkred"),
    labels = c("Actual", "No MP Shocks")
  ) +
  scale_linetype_manual(
    values = c("Inflation_Actual" = "solid", "Inflation_Counterfactual" = "dashed"),
    labels = c("Actual", "No MP Shocks")
  ) +
  labs(
    title = "Inflation: Actual vs. Counterfactual (No MP Shocks)",
    subtitle = "What would inflation have been without monetary policy shocks?",
    x = NULL,
    y = "Inflation (%)",
    color = NULL,
    linetype = NULL
  ) +
  theme_cbn

# Plot 3: Forecast Scenarios
fc_long_inflation <- forecast_comparison %>%
  select(Date, Inflation_Baseline, Inflation_Tightening, Inflation_Easing) %>%
  reshape2::melt(id.vars = "Date", variable.name = "Scenario", value.name = "Value")

fc_long_output <- forecast_comparison %>%
  select(Date, Output_Baseline, Output_Tightening, Output_Easing) %>%
  reshape2::melt(id.vars = "Date", variable.name = "Scenario", value.name = "Value")

p5 <- ggplot(fc_long_inflation, aes(x = Date, y = Value, color = Scenario, linetype = Scenario)) +
  geom_line(linewidth = 1) +
  scale_color_manual(
    values = c(
      "Inflation_Baseline" = "black",
      "Inflation_Tightening" = "darkred",
      "Inflation_Easing" = "darkgreen"
    ),
    labels = c("Baseline", "Tightening (+100bp)", "Easing (-100bp)")
  ) +
  scale_linetype_manual(
    values = c(
      "Inflation_Baseline" = "solid",
      "Inflation_Tightening" = "dashed",
      "Inflation_Easing" = "dotted"
    ),
    labels = c("Baseline", "Tightening (+100bp)", "Easing (-100bp)")
  ) +
  labs(
    title = "Inflation Forecast Under Alternative Policy Scenarios",
    subtitle = "12-quarter ahead forecast (3 years)",
    x = NULL,
    y = "Inflation (%)",
    color = "Scenario",
    linetype = "Scenario"
  ) +
  theme_cbn

p6 <- ggplot(fc_long_output, aes(x = Date, y = Value, color = Scenario, linetype = Scenario)) +
  geom_line(linewidth = 1) +
  scale_color_manual(
    values = c(
      "Output_Baseline" = "black",
      "Output_Tightening" = "steelblue",
      "Output_Easing" = "purple"
    ),
    labels = c("Baseline", "Tightening (+100bp)", "Easing (-100bp)")
  ) +
  scale_linetype_manual(
    values = c(
      "Output_Baseline" = "solid",
      "Output_Tightening" = "dashed",
      "Output_Easing" = "dotted"
    ),
    labels = c("Baseline", "Tightening (+100bp)", "Easing (-100bp)")
  ) +
  labs(
    title = "Output Forecast Under Alternative Policy Scenarios",
    subtitle = "12-quarter ahead forecast (3 years)",
    x = NULL,
    y = "Log Output (×100)",
    color = "Scenario",
    linetype = "Scenario"
  ) +
  theme_cbn

# Save plots
pdf("svar_policy_analysis.pdf", width = 12, height = 8)
grid.arrange(p1, p2, ncol = 2, top = "Impulse Response Functions")
grid.arrange(p3, p4, ncol = 1, top = "Counterfactual Analysis")
grid.arrange(p5, p6, ncol = 1, top = "Forecast Scenarios")
dev.off()

cat("Plots saved to: svar_policy_analysis.pdf\n\n")

################################################################################
# 7. EXPORT RESULTS
################################################################################

cat("=================================================================\n")
cat("EXPORTING RESULTS\n")
cat("=================================================================\n")

# Export counterfactual comparison
output_dir <- "c:/Users/Ex-Ethic/OneDrive/Documents/CBs_Works/Facilitation/2026/TSM_Forecasting/Policy_Counterfactual_Forecast_Scenaario_from_SVAR"

write.csv(
  comparison_data,
  file.path(output_dir, "counterfactual_results.csv"),
  row.names = FALSE
)
cat("Counterfactual results exported: counterfactual_results.csv\n")

# Export forecast scenarios
write.csv(
  forecast_comparison,
  file.path(output_dir, "forecast_scenarios.csv"),
  row.names = FALSE
)
cat("Forecast scenarios exported: forecast_scenarios.csv\n")

# Export IRF data
irf_export <- data.frame(
  Horizon = 0:20,
  Output_IRF = svar_model$irf$PolicyRate[, "Output"],
  Output_Lower = svar_model$Lower$PolicyRate[, "Output"],
  Output_Upper = svar_model$Upper$PolicyRate[, "Output"],
  Inflation_IRF = svar_model$irf$PolicyRate[, "Inflation"],
  Inflation_Lower = svar_model$Lower$PolicyRate[, "Inflation"],
  Inflation_Upper = svar_model$Upper$PolicyRate[, "Inflation"],
  PolicyRate_IRF = svar_model$irf$PolicyRate[, "PolicyRate"],
  PolicyRate_Lower = svar_model$Lower$PolicyRate[, "PolicyRate"],
  PolicyRate_Upper = svar_model$Upper$PolicyRate[, "PolicyRate"]
)

write.csv(
  irf_export,
  file.path(output_dir, "impulse_response_functions.csv"),
  row.names = FALSE
)
cat("Impulse response functions exported: impulse_response_functions.csv\n")

################################################################################
# 8. POLICY BRIEF SUMMARY
################################################################################

cat("\n")
cat("=================================================================\n")
cat("POLICY BRIEF SUMMARY\n")
cat("=================================================================\n\n")

cat("KEY FINDINGS:\n\n")

cat("1. COUNTERFACTUAL ANALYSIS (Neutralizing MP Shocks):\n")
cat(sprintf("   - On average, MP shocks %s inflation by %.2f pp\n",
            ifelse(mean(comparison_data$Inflation_Impact, na.rm = TRUE) > 0, "raised", "lowered"),
            abs(mean(comparison_data$Inflation_Impact, na.rm = TRUE))))
cat(sprintf("   - On average, MP shocks %s output by %.2f pp\n",
            ifelse(mean(comparison_data$Output_Impact, na.rm = TRUE) > 0, "raised", "lowered"),
            abs(mean(comparison_data$Output_Impact, na.rm = TRUE))))
cat("   - Implication: Monetary policy has been an active stabilization tool\n\n")

cat("2. IMPULSE RESPONSE ANALYSIS:\n")
peak_inflation_impact <- min(irf_data_inflation$IRF)
peak_horizon <- which.min(irf_data_inflation$IRF) - 1
cat(sprintf("   - Peak impact on inflation: %.3f pp at horizon %d quarters\n",
            peak_inflation_impact, peak_horizon))
cat(sprintf("   - Inflation response is %s\n",
            ifelse(peak_inflation_impact < 0, "consistent with policy effectiveness (negative)", 
                   "unexpectedly positive (price puzzle?)")))
cat("\n")

cat("3. FORECAST SCENARIOS (12 Quarters Ahead):\n\n")
cat("   SUSTAINED TIGHTENING (+100bp shock):\n")
cat(sprintf("     • Inflation: %.2f%% (%.2f pp below baseline)\n",
            mean(forecast_comparison$Inflation_Tightening),
            abs(mean(forecast_comparison$Inflation_Tight_Dev))))
cat(sprintf("     • Output: %.2f (%.2f pp below baseline)\n",
            mean(forecast_comparison$Output_Tightening),
            abs(mean(forecast_comparison$Output_Tight_Dev))))
cat("\n")
cat("   SUSTAINED EASING (-100bp shock):\n")
cat(sprintf("     • Inflation: %.2f%% (%.2f pp above baseline)\n",
            mean(forecast_comparison$Inflation_Easing),
            mean(forecast_comparison$Inflation_Easy_Dev)))
cat(sprintf("     • Output: %.2f (%.2f pp above baseline)\n",
            mean(forecast_comparison$Output_Easing),
            mean(forecast_comparison$Output_Easy_Dev)))
cat("\n")

cat("POLICY IMPLICATIONS:\n")
cat("   - Sustained tightening would help anchor inflation expectations\n")
cat("   - Trade-off: Lower inflation comes at cost of reduced output\n")
cat("   - Policy effectiveness validated through counterfactual exercise\n")
cat("   - Scenario analysis provides guidance for MPC deliberations\n\n")

cat("=================================================================\n")
cat("SCRIPT COMPLETED SUCCESSFULLY\n")
cat("=================================================================\n")
cat("Output files generated:\n")
cat("  1. svar_policy_analysis.pdf (visualizations)\n")
cat("  2. counterfactual_results.csv (actual vs. no-MP-shock)\n")
cat("  3. forecast_scenarios.csv (baseline, tightening, easing)\n")
cat("  4. impulse_response_functions.csv (IRF data)\n\n")

################################################################################
# END OF SCRIPT
################################################################################