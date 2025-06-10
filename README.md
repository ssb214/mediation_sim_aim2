# mediation_sim_aim2
This repository contains code for a simulation study evaluating the performance of three high-dimensional mediation analysis methods (HIMA, HDMA, and Meet-in-the-Middle [MITM]) in the context of untargeted metabolomics data.

# Overview
The goal of this study is to assess how each method estimates:

- Component Indirect Effects (CIEs)
- Total Indirect Effects (TIEs)
- Sensitivity and Specificity in identifying true mediators

Simulations were conducted under 36 scenarios for each of three mediator set sizes (p = 200, 400, 600), varying by:
- Mediator correlation structure (independent vs. correlated)
- Sample size (n = 200, 500, 1000)
- Proportion of true mediators (2%, 5%, 10%)
- Mediator effect size (β = 0.1, 0.3)
