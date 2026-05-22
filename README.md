# Brief Guide on Policy Counterfactual Analysis and Conditional Forecast Scenarios in SVAR

Policy counterfactual analysis and conditional forecasting are among the most important applications of Econometrics and Structural Vector Autoregression models in central banking and macroeconomic policy analysis. They allow policymakers to answer questions such as:

* *What would inflation have been if the policy rate had not increased?*
* *What if oil prices fall by 30% next quarter?*
* *What would happen if exchange rate pressures intensify?*
* *How would monetary policy react under alternative scenarios?*

These frameworks are heavily used in institutions such as the International Monetary Fund, Bank of England, European Central Bank, and Federal Reserve System.

---

# 1. Conceptual Foundations

## A. Structural VAR (SVAR)

An SVAR extends the standard VAR by imposing economic structure that allows identification of structural shocks.

Standard reduced-form VAR:

Y_t = A_1Y_{t-1}+A_2Y_{t-2}+\cdots+A_pY_{t-p}+u_t

Structural representation:

B_0Y_t = B_1Y_{t-1}+\cdots+B_pY_{t-p}+\varepsilon_t

Where:

* (Y_t): vector of endogenous macro variables
* (u_t): reduced-form residuals
* (\varepsilon_t): structural shocks
* (B_0): contemporaneous structural relationships

Typical macro variables include:

* Inflation
* Output gap
* Policy rate
* Exchange rate
* Oil price
* Money supply
* Credit
* External reserves

---

# 2. Policy Counterfactual Analysis in SVAR

## Definition

Counterfactual analysis asks:

> “What would macroeconomic outcomes have been under an alternative policy path or alternative shock realization?”

The idea is to alter one structural shock or policy instrument while keeping the rest of the system unchanged.

---

# 3. Types of Counterfactual Analysis

## A. Monetary Policy Counterfactual

Example:

* Actual policy:

  * MPR increased from 18% to 27.5%
* Counterfactual:

  * Assume no tightening

Questions answered:

* Would inflation have been higher?
* Would exchange rate depreciation worsen?
* What would happen to output growth?

---

## B. External Shock Counterfactual

Examples:

* Oil price collapse
* Global food inflation shock
* US interest rate hike
* Capital outflow episode

---

## C. Exchange Rate Policy Counterfactual

Examples:

* Fixed exchange rate vs managed float
* Faster FX pass-through containment
* Alternative reserve intervention strategies

---

# 4. General Procedure for Counterfactual Analysis in SVAR

## Step 1: Estimate the SVAR

Estimate the model using:

* Bayesian SVAR
* Classical SVAR
* Sign-restricted SVAR
* Recursive (Cholesky) SVAR
* Proxy SVAR

---

## Step 2: Identify Structural Shocks

Common identification approaches:

| Method                 | Description                  |
| ---------------------- | ---------------------------- |
| Cholesky decomposition | Recursive ordering           |
| Sign restrictions      | Economic theory restrictions |
| Long-run restrictions  | Blanchard-Quah style         |
| External instruments   | Proxy SVAR                   |
| Narrative restrictions | Historical episodes          |

---

## Step 3: Recover Structural Shocks

Recover:

\varepsilon_t = B_0u_t

These shocks represent:

* Monetary policy shocks
* Supply shocks
* Demand shocks
* Exchange rate shocks
* Oil shocks

---

## Step 4: Shut Down or Modify Specific Shocks

For example:

### Monetary policy counterfactual

Set:

\varepsilon_t^{MP}=0

Meaning:

* No monetary policy tightening shock occurred.

Alternative:

* Reduce shock magnitude by 50%
* Replace with expansionary shock
* Impose alternative policy path

---

## Step 5: Simulate Counterfactual Paths

Generate simulated series:

* Inflation
* GDP
* Exchange rate
* Credit
* Reserves

using recursive SVAR dynamics.

---

# 5. Conditional Forecast Scenario Analysis

## Definition

Conditional forecasting asks:

> “Given a specified future path for some variables, what is the implied forecast for the remaining variables?”

This is widely used in:

* Inflation targeting
* Stress testing
* Policy scenario analysis
* FPAS/QPM systems

---

# 6. Unconditional vs Conditional Forecast

| Forecast Type          | Meaning                                 |
| ---------------------- | --------------------------------------- |
| Unconditional forecast | Pure model-based forecast               |
| Conditional forecast   | Forecast subject to imposed assumptions |

---

# 7. Examples of Conditional Scenarios

## Scenario 1: Oil Price Shock

Assume:

* Brent oil price falls by 25%

Then forecast:

* Inflation
* Exchange rate
* Fiscal balance
* GDP growth

---

## Scenario 2: Monetary Tightening

Assume:

* Policy rate rises by 300 bps

Forecast implications for:

* Core inflation
* Credit
* Output gap
* FX pressures

---

## Scenario 3: Exchange Rate Shock

Assume:

* 20% depreciation of naira

Forecast:

* Imported inflation
* Food inflation
* Reserve depletion
* Pass-through effects

---

# 8. Mathematical Representation of Conditional Forecast

State-space representation:

Y_{t+h}=\mu + \sum_{i=1}^{p}A_iY_{t+h-i}+\varepsilon_{t+h}

Conditional forecast imposes restrictions such as:

[
i_t = i_t^*
]

or

[
e_t = e_t^*
]

where:

* (i_t^*): imposed policy rate path
* (e_t^*): imposed exchange rate path

The model solves for all remaining endogenous variables consistently.

---

# 9. Techniques for Conditional Forecasting

## A. Waggoner-Zha Conditional Forecasting

Classic approach used in many central banks.

Main reference:

* Waggoner & Zha (1999)

Features:

* Bayesian framework
* Conditional paths
* Shock decomposition
* Scenario consistency

---

## B. Kalman Filter / State-Space Methods

Used especially in:

* DSGE models
* Time-varying parameter SVARs
* FPAS systems

---

## C. Bayesian Conditional Forecasting

Advantages:

* Incorporates uncertainty
* Density forecasts
* Fan charts
* Scenario probabilities

---

# 10. Key Outputs from Counterfactual and Conditional Analysis

## A. Impulse Response Functions (IRFs)

Show dynamic response to shocks.

Example:

* Inflation response to monetary tightening

---

## B. Historical Decomposition

Explains which shocks drove inflation historically.

Example:

* Exchange rate shocks contributed 40% of inflation surge.

---

## C. Forecast Error Variance Decomposition (FEVD)

Measures relative importance of shocks.

---

## D. Scenario Fan Charts

Probabilistic forecast ranges.

---

# 11. Applications to Nigeria

For the Central Bank of Nigeria, important applications include:

## Monetary Policy

* MPR tightening counterfactual
* Inflation targeting transition
* Exchange rate stabilization

## External Sector

* Oil price scenarios
* Global financial tightening
* Capital flow reversals

## Fiscal-Monetary Coordination

* Fuel subsidy removal
* Fiscal dominance analysis
* Debt monetization scenarios

## Inflation Analysis

* Food inflation decomposition
* FX pass-through analysis
* Imported inflation scenarios

---

# 12. Recommended SVAR Variables for Nigeria

| Block       | Variables                               |
| ----------- | --------------------------------------- |
| Monetary    | MPR, money supply, treasury bill rate   |
| Inflation   | Headline CPI, core CPI, food CPI        |
| Real sector | GDP, PMI, industrial production         |
| External    | Oil price, reserves, exchange rate      |
| Financial   | Credit, stock market, bond yields       |
| Global      | Fed Funds rate, global commodity prices |

---

# 13. Software Commonly Used

## R Packages

* [vars package](https://cran.r-project.org/package=vars?utm_source=chatgpt.com)
* [BVAR package](https://cran.r-project.org/package=BVAR?utm_source=chatgpt.com)
* [svars package](https://cran.r-project.org/package=svars?utm_source=chatgpt.com)
* [bvartools package](https://cran.r-project.org/package=bvartools?utm_source=chatgpt.com)

---

## Python Libraries

* [statsmodels](https://www.statsmodels.org?utm_source=chatgpt.com)
* [PyMC](https://www.pymc.io?utm_source=chatgpt.com)
* [NumPy](https://numpy.org?utm_source=chatgpt.com)

---

## MATLAB Toolboxes

* BEAR Toolbox
* Dynare
* IRIS Toolbox

---

# 14. Important Practical Challenges

| Challenge              | Explanation                          |
| ---------------------- | ------------------------------------ |
| Identification problem | Structural shocks may not be unique  |
| Lucas critique         | Policy regime changes alter behavior |
| Parameter instability  | Nigerian macro data often unstable   |
| Small samples          | Quarterly data limitations           |
| Structural breaks      | FX reforms, subsidy reforms, crises  |

---

# 15. Best-Practice Recommendations

For Nigeria-specific SVAR policy analysis:

1. Use Bayesian SVAR rather than purely classical SVAR.
2. Include oil-price and exchange-rate channels explicitly.
3. Allow for structural breaks.
4. Incorporate stochastic volatility if possible.
5. Use density forecasts rather than only point forecasts.
6. Combine SVAR with judgmental adjustments.
7. Complement SVAR with QPM/FPAS framework.

---

# Key References

## Foundational SVAR References

### Christopher Sims

* Sims, C. A. (1980)
* *Macroeconomics and Reality*
* Econometrica

---

### Ben Bernanke

* Bernanke, Gertler & Watson (1997)
* *Systematic Monetary Policy and the Effects of Oil Price Shocks*

---

### Olivier Blanchard and Danny Quah

* Blanchard & Quah (1989)
* *The Dynamic Effects of Aggregate Demand and Supply Disturbances*

---

# Counterfactual and Conditional Forecasting References

### Daniel Waggoner and Tao Zha

* Waggoner, D. & Zha, T. (1999)
* *Conditional Forecasts in Dynamic Multivariate Models*

---

### Tao Zha

* Zha (1999)
* *Block Recursion and Structural Vector Autoregressions*

---

### Frank Smets and Raf Wouters

* DSGE and scenario forecasting applications

---

# Bayesian SVAR References

### Gary Koop and Dimitris Korobilis

* *Bayesian Multivariate Time Series Methods for Empirical Macroeconomics*

---

### Helmut Lütkepohl

* *New Introduction to Multiple Time Series Analysis*

---

# Central Bank and IMF References

## [IMF FPAS Notes](https://www.imf.org/en/Publications/WP?utm_source=chatgpt.com)

Search for:

* FPAS
* Conditional forecasting
* Scenario analysis
* Inflation targeting

---

## [BEAR Toolbox](https://www.ecb.europa.eu/pub/research/working-papers/html/bear-toolbox.en.html?utm_source=chatgpt.com)

Excellent for:

* BVAR
* Conditional forecasts
* Counterfactual simulations
* Shock decomposition

---

# Suggested Advanced Extensions

You may later extend the framework to:

* TVP-SVAR
* SVAR-SV
* Factor-Augmented SVAR (FAVAR)
* Proxy SVAR
* Markov-switching SVAR
* DSGE-VAR hybrid models
* Semi-structural gap models
* FPAS/QPM integration

These are especially useful for Nigeria due to:

* Regime shifts
* Oil dependence
* Exchange-rate volatility
* Structural breaks
* Monetary policy transition dynamics
