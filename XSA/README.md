# Extended Survivors Analysis (XSA)

---

## What is XSA?

**Extended Survivors Analysis (XSA)** is an extension of **Virtual Population Analysis (VPA)** used in fisheries science to reconstruct population-at-age (`N_ay`) from **catch-at-age** (`C_ay`) and **survey indices** (`u_ay`).

Unlike VPA, XSA incorporates:
- **Abundance indices** to tune fishing mortality estimates,
- **Cohort-based survivors** (`P_t(k)`) to reconstruct population-at-age (`N_ay`),
- **Weighted fitting** to downweight old data.

It is particularly useful when:
- There’s uncertainty in terminal fishing mortality,
- Abundance surveys vary in reliability over time,
- You need a more robust population reconstruction.

---

## How it Works
This implementation is based on the 1999 Shepherd's paper *Extended Survivors Analysis: an improved method for the estimation of stock abundance from catch and survey data*. Thus, the method does not include shrinkage of age - qage, rage - nor other more recent tunings.

XSA estimates `N_ay` by:
1. Starting with a VPA estimate of `N_ay`.
2. Using abundance indices to tune the catchability coefficient.
3. Minimizing the difference between the VPA-derived and index-derived `N_ay`.

Step 3 leads to an estimate of the survivors for each cohort `P_t(k)` which produces a diagonal reconstruction of `N_ay` which is then completed **using VPA** for missing entries.

The fishing mortality matrix F_ay is reconstructed via the population dynamics equation. Since this estimate depends on `N_ay`, which may be noisy, we specify a minimum fishing mortality threshold to avoid unrealistic growth of population abundance.
An alternative approach would be to calculate F_ay through Baranov's equation - not implemented here.

---

## Understanding P_t(k)
`P_t(k)` is the terminal population (hence subscript t) at the end of the final year of the cohort, i.e. the final *age* of the cohort, namely the **survivors of that cohort** after all ages have been observed.
The cohort index is `k = y - a`, where y is year and a is age. In the example in the notebook we have `a = 5`, `y = 7`, meaning ages range from age 0 to age 4 and we assume years range from 1990 to 1996. We therefore see in our data 11 cohorts, from 1986 to 1996; however, we only observe cohorts 1990-1996 from the start, so `P_t(k)` estimates the survivors of these cohorts. 
For example, `P_t[6]` is the terminal population of cohort 1996 at the age of 4, meaning at year 2000. This is why using P_t(k) leads to a diagonal estimate of `N_ay`: we only get values where the age and year line up to form a cohort, i.e., `y - a = k`. Therefore, cohorts before 1990 will not have their `P_t(k)` estimated, and their diagonals in `N_ay` will remain unfilled unless reconstructed via VPA.


**Author:** Nazareno Campioni  
**Based on:** Shepherd (1999) – *Extended Survivors Analysis: an improved method for the estimation of stock abundance from catch and survey data*  
**Language:** Python 3  
**Repository:** https://github.com/naza-campioni/fishery-models-tutorial/XSA
