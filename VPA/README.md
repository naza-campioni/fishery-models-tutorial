# Virtual Population Analysis

Virtual Population Analysis (VPA) is a cohort-based method used in fisheries science to estimate population-at-age `N_ay` and fishing mortality `F_ay` from observed catch-at-age `C_ay`. It is the foundation of many stock assessment models.

### Key assumptions
1. The total mortality `Z` is split into natural `M` and fishing `F` mortality: `Z_ay = M_ay + F_ay`.
2. Terminal fishing mortality `F_AY` is known.

### VPA equations
1. Terminal year estimate (from Baranov's equation):
     `N_AY = C_AY/(1 - exp(-Z_AY) * (Z_AY/F_AY)`,
   where subscript `AY` indicates ginal year.

2. Recursive backward estimation assuming catches are taken mid-year (Pope's approximation):
     `N_ay = N_{a+1,y+1}*e^Z + C_ay*e^{Z/2}`

