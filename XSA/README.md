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


**Author:** Nazareno Campioni  
**Based on:** Shepherd (1999) – *Extended Survivors Analysis: an improved method for the estimation of stock abundance from catch and survey data*  
**Language:** Python 3  
**Repository:** https://github.com/naza-campioni/fishery-models-tutorial/XSA
