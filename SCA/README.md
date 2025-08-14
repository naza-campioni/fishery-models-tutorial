# Statistical Catch-at-Age (SCA)

---

## Overview
A **Statistical Catch-at-Age (SCA) model** is a type of fish stock assessment model that uses both **total catch data** and **age composition data** (proportions-at-age or numbers-at-age) to estimate key population parameters such as:

- Recruitment (number of new fish entering the population each year)
- Fishing mortality (F)
- Abundance-at-age (N)
- Selectivity-at-age (s)
- Catchability (q) for survey indices

SCA models combine **population dynamics equations** with **statistical likelihoods** to find the parameter set that best explains observed data.

---

## Key Concepts

### 1. Population Dynamics
We model numbers-at-age using the **cohort survival equations**:

\[
N_{a+1, y+1} =
\begin{cases}
N_{a, y} e^{-(s_a F_y + M_{a, y})}, & a < A_{\text{max}} \\
N_{a, y} e^{-(s_a F_y + M_{a, y})} + N_{a, y} e^{-(s_a F_y + M_{a, y})}, & a = A_{\text{max}}
\end{cases}
\]

Where:
- \(N_{a,y}\) = numbers-at-age \(a\) in year \(y\)
- \(s_a\) = selectivity-at-age
- \(F_y\) = fishing mortality in year \(y\)
- \(M_{a,y}\) = natural mortality

---

### 2. Catch Equations
Catch-at-age is computed using **Baranov's equation**:

\[
C_{a,y} = \frac{F_{a,y}}{Z_{a,y}} N_{a,y} \left( 1 - e^{-Z_{a,y}} \right),
\]
where \(Z_{a,y} = F_{a,y} + M_{a,y}\) is total mortality.

---

### 3. Data Types Used
- **Total catch** (\(C_{\text{tot},y}\)): sum over ages
- **Proportions-at-age** (\(p_{a,y}\)): \(C_{a,y} / C_{\text{tot},y}\)
- **Survey indices** (\(I_{a,y}\)): abundance estimates from independent surveys

---

## Likelihood Components

An SCA model combines multiple likelihood terms:

1. **Lognormal for total catch**:
\[
C_{\text{tot},y}^{\text{obs}} = C_{\text{tot},y}^{\text{pred}} \times \varepsilon, \quad \ln \varepsilon \sim \mathcal{N}(0, \sigma^2_{\text{tot}})
\]

2. **Lognormal for survey indices**:
\[
I_{a,y}^{\text{obs}} = I_{a,y}^{\text{pred}} \times \varepsilon, \quad \ln \varepsilon \sim \mathcal{N}(0, \sigma^2_I)
\]

3. **Multinomial for proportions-at-age**:
\[
\mathbf{C}_y \sim \text{Multinomial}(n_y, \mathbf{p}_y)
\]
where \(n_y\) is the effective sample size (number of fish aged), and \(\mathbf{p}_y\) are predicted proportions.

---

## Parameters Estimated
Typical estimated parameters include:
- **Selectivity parameters** (\(a_{50}\), \(a_{95}\)) for logistic selectivity curve
- **Annual fishing mortality** (\(F_y\))
- **Recruitment (initial N)** for each year
- **Catchability coefficient** (\(q\))
- **Observation error variances** (\(\sigma_{\text{tot}}, \sigma_I\))
- **(Optional)** Effective sample size for proportions (\(n_y\)) if not fixed

---

## Workflow
1. **Input data**:
   - Catches-at-age
   - Total catch
   - Survey indices
   - Natural mortality
   - Timing of survey relative to fishing (ts)

2. **Model population numbers forward in time** using recruitment and mortality estimates.

3. **Predict observations** (catch, proportions, survey index).

4. **Compute likelihood** comparing predictions to observed data.

5. **Optimize parameters** by minimizing negative log-likelihood.

---

## Tips for Stability
- If estimating \(n_y\), use a weak prior or fix it to the sampling design value.
- Start with reasonable initial guesses for F and N.
- Check that selectivity parameters produce biologically realistic shapes.
- Be cautious of parameter confounding (especially N and q).

---

## References
- Haddon, M. (2001). *Modelling and Quantitative Methods in Fisheries*. Chapman & Hall/CRC.
- Quinn, T. J., & Deriso, R. B. (1999). *Quantitative Fish Dynamics*. Oxford University Press.

---


