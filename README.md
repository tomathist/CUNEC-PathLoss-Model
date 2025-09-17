# CUNEC: Corner-Urban Neighborhood Channel Path-Loss Model (MATLAB)

CUNEC is a modular MATLAB implementation of a **corner-aware** urban path-loss model with **jointly correlated log-normal shadowing** across UEs and APs. It provides:
- **Zeroth-order** baseline loss (distance-only)
- **First-order** corner contribution (with near/far gating)
- **Second-order** further corner/leg contribution
- UE–AP **joint correlation** of shadowing via anisotropic fields

The code is designed for cell-free/dense-AP studies and generates paper-quality figures with Monte-Carlo (MC) bands.

---

## Repository Structure
.
├─ demos/
│ ├─ demo_cunec_zeroth_shadowing.m
│ ├─ demo_cunec_firstorder_shadowing.m
│ └─ demo_cunec_secondorder_shadowing.m
├─ src/
│ ├─ calc_PL_CUNEC_0th.m % [D, PL0] = ...
│ ├─ calc_PL_CUNEC_1st.m % [D, PL1] = ...
│ └─ calc_PL_CUNEC_2nd.m % [D, PL2] = ...
├─ utils/
│ └─ gen_shadowing_joint_aniso.m % UE–AP jointly correlated shadowing (dB)
├─ params/
│ ├─ load_model_parameters.m % FSPL_1m, mu_, sigma_ (user-provided)
│ └─ load_correlations.m % C_* correlation matrices (user-provided)
├─ LICENSE
└─ README.md
