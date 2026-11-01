---
title: "Pressure growth curve"
author: "Jay T. Lennon"
date: "`r format(Sys.time(), '%d %B, %Y')`"
header-includes:
  - \usepackage{array}
output: pdf_document
geometry: margin=2.54cm
---
  
  ## Clear workspace and set directory
```{r setup}
rm(list=ls())
getwd()
setwd("~/GitHub/MetaboPressure")

# Coupled bacterial growth, respiration, gas transfer, and headspace pressure
# -------------------------------------------------------------------------
# Illustrative model for a closed, aerobic batch culture. It uses deSolve for
# the ODE integration and base R graphics for optional plotting.
#
# State variables and units
#   time       h
#   X          g dry biomass / L
#   S          g carbon substrate / L
#   O          mol dissolved O2 / L
#   DIC        mol dissolved inorganic carbon / L (CO2* + HCO3- + CO3--)
#   nO2, nCO2  mol in headspace
#   a          dimensionless lag activation (0 to 1)
#   O2_used    cumulative mol O2 consumed
#   CO2_made   cumulative mol CO2 produced
#
# IMPORTANT: Defaults are demonstration values, not fitted Bacillus parameters.
# Replace them with measurements or fit them to your data before interpreting
# simulated values quantitatively. Time is in hours, pressure in atm, volumes
# in L, and concentrations in mol/L unless the variable names say otherwise.

if (!requireNamespace("deSolve", quietly = TRUE)) {
  stop("Install deSolve first with install.packages('deSolve').")
}

# ------------------------------ User settings -----------------------------
cfg <- list(
  # Vessel: 125 mL total jar with approximately 75 mL headspace.
  # Include all connected syringe/tubing/sensor dead volume in Vgas_eff_L.
  Vgas_eff_L = 0.075,
  Vliquid_L = 0.050,
  shaker_rpm = 200,
  temperature_K = 298.15,
  pressure_initial_atm = 1.0,
  water_vapor_pressure_atm = 0.0313, # approximate at 25 C

  # Initial dry headspace composition (air, approximately).
  yO2_initial = 0.2095,
  yCO2_initial = 0.0004,

  # Fixed-pH approximation for the carbonate system.
  pH_fixed = 7.5,
  pKa1 = 6.35,
  pKa2 = 10.33,

  # Henry solubility coefficients at approximately 25 C:
  # dissolved concentration (mol/L) = H * gas partial pressure (atm).
  H_O2_mol_L_atm = 0.0013,
  H_CO2_mol_L_atm = 0.033,
  # Effective values for the actual bottle geometry, fill volume, and shaking
  # condition (about 200 rpm here). Shaking generally raises gas-liquid
  # transfer, but rpm alone does not determine kLa; measure or fit these.
  kLa_O2_per_h = 1.0,
  kLa_CO2_per_h = 1.0,

  # Growth and substrate use. These are generic starting values only.
  mu_max_per_h = 0.60,
  K_substrate_g_L = 0.05,
  K_O2_mol_L = 1.0e-5,
  yield_XS_gX_per_gS = 0.50,
  maintenance_S_gS_gX_h = 0.005,
  lag_time_constant_h = 2.0,
  decay_per_h = 0.0,

  # Respiratory demand. qO2_growth is mol O2 per g biomass formed;
  # maintenance_O2 is mol O2 per g biomass per h.
  qO2_growth_mol_gX = 0.08,
  maintenance_O2_mol_gX_h = 0.001,
  respiratory_quotient_molCO2_per_molO2 = 1.0,

  # Initial culture state.
  X0_g_L = 0.01,
  S0_g_L = 2.0,
  activation0 = 0.0,

  # Integration and optional output.
  time_end_h = 48,
  time_step_h = 0.05,
  make_plots = TRUE,
  save_plot_pdf = FALSE,
  plot_pdf_file = "microbial_pressure_model.pdf",
  save_csv = FALSE,
  csv_file = "microbial_pressure_model_output.csv"
)

R_gas_L_atm_mol_K <- 0.082057366

# Carbonate speciation multiplier at the specified fixed pH. The dissolved
# molecular CO2 concentration is DIC / alpha. At pH 7.5 most DIC is HCO3-.
carbonate_alpha <- function(pH, pKa1, pKa2) {
  1 + 10^(pH - pKa1) + 10^(2 * pH - pKa1 - pKa2)
}
alpha_DIC <- carbonate_alpha(cfg$pH_fixed, cfg$pKa1, cfg$pKa2)

# Initial dry-gas pressure excludes water vapor; the constant vapor term is
# added back when calculating total pressure.
Pdry0_atm <- cfg$pressure_initial_atm - cfg$water_vapor_pressure_atm
if (cfg$yO2_initial + cfg$yCO2_initial >= 1) {
  stop("Initial O2 and CO2 fractions must sum to less than one.")
}
nO2_0_mol <- cfg$yO2_initial * Pdry0_atm * cfg$Vgas_eff_L /
  (R_gas_L_atm_mol_K * cfg$temperature_K)
nCO2_0_mol <- cfg$yCO2_initial * Pdry0_atm * cfg$Vgas_eff_L /
  (R_gas_L_atm_mol_K * cfg$temperature_K)
nN2_other_0_mol <- (1 - cfg$yO2_initial - cfg$yCO2_initial) * Pdry0_atm *
  cfg$Vgas_eff_L / (R_gas_L_atm_mol_K * cfg$temperature_K)

# Assume the liquid starts in gas-liquid equilibrium with the initial headspace.
O0_mol_L <- cfg$H_O2_mol_L_atm * cfg$yO2_initial * Pdry0_atm
DIC0_mol_L <- alpha_DIC * cfg$H_CO2_mol_L_atm *
  cfg$yCO2_initial * Pdry0_atm

initial_state <- c(
  X = cfg$X0_g_L,
  S = cfg$S0_g_L,
  O = O0_mol_L,
  DIC = DIC0_mol_L,
  nO2 = nO2_0_mol,
  nCO2 = nCO2_0_mol,
  activation = cfg$activation0,
  O2_used = 0,
  CO2_made = 0
)

parameters <- cfg
parameters$nN2_other_mol <- nN2_other_0_mol
parameters$alpha_DIC <- alpha_DIC
parameters$R_gas_L_atm_mol_K <- R_gas_L_atm_mol_K

# ------------------------------- ODE system -------------------------------
microbe_gas_model <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    # Keep rates physical if a numerical solver briefly evaluates a state just
    # below zero near complete substrate or oxygen depletion.
    X_pos <- max(X, 0)
    S_pos <- max(S, 0)
    O_pos <- max(O, 0)
    DIC_pos <- max(DIC, 0)
    nO2_pos <- max(nO2, 0)
    nCO2_pos <- max(nCO2, 0)
    a_pos <- min(max(activation, 0), 1)

    # Headspace partial pressures follow the ideal gas law.
    pO2_atm <- nO2_pos * R_gas_L_atm_mol_K * temperature_K / Vgas_eff_L
    pCO2_atm <- nCO2_pos * R_gas_L_atm_mol_K * temperature_K / Vgas_eff_L

    # Monod limitation by carbon and dissolved oxygen; activation generates a
    # gradual lag-to-growth transition.
    fS <- S_pos / (K_substrate_g_L + S_pos)
    fO <- O_pos / (K_O2_mol_L + O_pos)
    mu <- mu_max_per_h * a_pos * fS * fO

    # Growth-associated plus maintenance oxygen demand. O2 limitation also
    # scales maintenance, preventing uptake after oxygen is exhausted.
    rO2_mol_h <- (qO2_growth_mol_gX * mu +
                    maintenance_O2_mol_gX_h * fO) * X_pos * Vliquid_L
    rCO2_mol_h <- respiratory_quotient_molCO2_per_molO2 * rO2_mol_h

    # CO2 gas transfer responds to molecular dissolved CO2, not total DIC.
    CO2_star_mol_L <- DIC_pos / alpha_DIC
    J_O2_gas_to_liquid_mol_h <- kLa_O2_per_h * Vliquid_L *
      (H_O2_mol_L_atm * pO2_atm - O_pos)
    J_CO2_gas_to_liquid_mol_h <- kLa_CO2_per_h * Vliquid_L *
      (H_CO2_mol_L_atm * pCO2_atm - CO2_star_mol_L)

    # Carbon substrate is consumed for growth and substrate maintenance.
    rS_g_L_h <- if (S_pos > 0) {
      (mu / yield_XS_gX_per_gS + maintenance_S_gS_gX_h * fO) * X_pos
    } else {
      0
    }

    dX <- (mu - decay_per_h) * X_pos
    dS <- -rS_g_L_h
    dO <- (J_O2_gas_to_liquid_mol_h - rO2_mol_h) / Vliquid_L
    dDIC <- (J_CO2_gas_to_liquid_mol_h + rCO2_mol_h) / Vliquid_L
    dnO2 <- -J_O2_gas_to_liquid_mol_h
    dnCO2 <- -J_CO2_gas_to_liquid_mol_h
    dActivation <- (1 - a_pos) / lag_time_constant_h

    list(c(dX, dS, dO, dDIC, dnO2, dnCO2, dActivation,
           rO2_mol_h, rCO2_mol_h),
         mu_per_h = mu,
         fO2 = fO,
         rO2_mol_h = rO2_mol_h,
         rCO2_mol_h = rCO2_mol_h,
         rS_g_L_h = rS_g_L_h,
         pO2_headspace_atm = pO2_atm,
         pCO2_headspace_atm = pCO2_atm,
         CO2_aq_mol_L = CO2_star_mol_L)
  })
}

times_h <- seq(0, cfg$time_end_h, by = cfg$time_step_h)
solution <- deSolve::ode(
  y = initial_state,
  times = times_h,
  func = microbe_gas_model,
  parms = parameters,
  method = "lsoda"
)
out <- as.data.frame(solution)

# Reconstruct measured total pressure from all headspace gases plus a constant
# water-vapor pressure. The nitrogen/other inert gas amount is held constant.
out$pO2_atm <- out$nO2 * R_gas_L_atm_mol_K * cfg$temperature_K /
  cfg$Vgas_eff_L
out$pCO2_atm <- out$nCO2 * R_gas_L_atm_mol_K * cfg$temperature_K /
  cfg$Vgas_eff_L
out$pN2_other_atm <- nN2_other_0_mol * R_gas_L_atm_mol_K *
  cfg$temperature_K / cfg$Vgas_eff_L
out$pressure_atm <- out$pO2_atm + out$pCO2_atm + out$pN2_other_atm +
  cfg$water_vapor_pressure_atm
out$pressure_kPa <- out$pressure_atm * 101.325
out$O_mg_L <- out$O * 32.00 * 1000
out$DIC_mmol_L <- out$DIC * 1000

# ------------------------------ Optional plots ----------------------------
if (isTRUE(cfg$make_plots)) {
  old_par <- par(no.readonly = TRUE)
  par(mfrow = c(3, 2), mar = c(4, 4, 2.5, 1))

  plot(out$time, out$X, type = "l", lwd = 2, col = "forestgreen",
       xlab = "Time (h)", ylab = "Biomass (g/L)",
       main = "Biomass")
  plot(out$time, out$S, type = "l", lwd = 2, col = "darkgoldenrod3",
       xlab = "Time (h)", ylab = "Substrate (g/L)",
       main = "Carbon substrate")
  plot(out$time, out$O_mg_L, type = "l", lwd = 2, col = "steelblue",
       xlab = "Time (h)", ylab = "Dissolved O2 (mg/L)",
       main = "Dissolved oxygen")
  plot(out$time, out$pO2_atm, type = "l", lwd = 2, col = "steelblue",
       ylim = range(c(out$pO2_atm, out$pCO2_atm)),
       xlab = "Time (h)", ylab = "Partial pressure (atm)",
       main = "Headspace gases")
  lines(out$time, out$pCO2_atm, lwd = 2, col = "firebrick")
  legend("right", legend = c("O2", "CO2"),
         col = c("steelblue", "firebrick"), lwd = 2, bty = "n")
  plot(out$time, out$pressure_kPa, type = "l", lwd = 2, col = "black",
       xlab = "Time (h)", ylab = "Total pressure (kPa)",
       main = "Predicted sensor pressure")
  abline(h = out$pressure_kPa[1], lty = 2, col = "grey50")
  plot(out$time, out$rO2_mol_h, type = "l", lwd = 2, col = "purple4",
       xlab = "Time (h)", ylab = "O2 uptake (mol/h)",
       main = "Respiration, including maintenance")

  if (isTRUE(cfg$save_plot_pdf)) {
    grDevices::pdf(cfg$plot_pdf_file, width = 10, height = 12)
    # Re-draw all panels on the file device.
    par(mfrow = c(3, 2), mar = c(4, 4, 2.5, 1))
    plot(out$time, out$X, type = "l", lwd = 2, col = "forestgreen",
         xlab = "Time (h)", ylab = "Biomass (g/L)", main = "Biomass")
    plot(out$time, out$S, type = "l", lwd = 2, col = "darkgoldenrod3",
         xlab = "Time (h)", ylab = "Substrate (g/L)", main = "Carbon substrate")
    plot(out$time, out$O_mg_L, type = "l", lwd = 2, col = "steelblue",
         xlab = "Time (h)", ylab = "Dissolved O2 (mg/L)", main = "Dissolved oxygen")
    plot(out$time, out$pO2_atm, type = "l", lwd = 2, col = "steelblue",
         ylim = range(c(out$pO2_atm, out$pCO2_atm)),
         xlab = "Time (h)", ylab = "Partial pressure (atm)", main = "Headspace gases")
    lines(out$time, out$pCO2_atm, lwd = 2, col = "firebrick")
    legend("right", legend = c("O2", "CO2"),
           col = c("steelblue", "firebrick"), lwd = 2, bty = "n")
    plot(out$time, out$pressure_kPa, type = "l", lwd = 2, col = "black",
         xlab = "Time (h)", ylab = "Total pressure (kPa)", main = "Predicted sensor pressure")
    abline(h = out$pressure_kPa[1], lty = 2, col = "grey50")
    plot(out$time, out$rO2_mol_h, type = "l", lwd = 2, col = "purple4",
         xlab = "Time (h)", ylab = "O2 uptake (mol/h)",
         main = "Respiration, including maintenance")
    grDevices::dev.off()
  }
  par(old_par)
}

if (isTRUE(cfg$save_csv)) {
  utils::write.csv(out, cfg$csv_file, row.names = FALSE)
}

cat("Simulation complete.\n")
cat(sprintf("Initial pressure: %.3f kPa\n", out$pressure_kPa[1]))
cat(sprintf("Final pressure:   %.3f kPa\n", tail(out$pressure_kPa, 1)))
cat(sprintf("Final biomass:    %.4g g/L\n", tail(out$X, 1)))
cat(sprintf("Final substrate:  %.4g g/L\n", tail(out$S, 1)))
cat(sprintf("O2 consumed:      %.4g mmol\n", tail(out$O2_used, 1) * 1000))
cat(sprintf("CO2 produced:     %.4g mmol\n", tail(out$CO2_made, 1) * 1000))

# Model interpretation and limitations
# -------------------------------------
# 1. The pH-fixed carbonate approximation makes DIC speciation algebraic. It
#    does not predict pH drift or buffer exhaustion; replace it with alkalinity
#    and carbonate-equilibrium equations if those matter to your experiment.
# 2. RQ is fixed here. Organic-acid/byproduct formation, overflow metabolism,
#    sporulation, or oxygen-limited metabolism can change CO2/O2 stoichiometry.
# 3. There are no explicit volatile products, evaporation, leaks, sensor drift,
#    or temperature dynamics. These can matter in a sensitive pressure trace.
# 4. Pressure alone cannot uniquely identify biomass, oxygen uptake, CO2
#    dissolution, and RQ. Fit against independent OD/CFU, substrate, pH, or gas
#    composition data where possible.
# 5. The lag state is a simple phenomenological ramp, not a mechanistic model
#    of spore germination. For Bacillus spores, a germination/vegetative-cell
#    transition may need explicit states.
# 6. The stated ~200 rpm shaking condition is recorded as metadata. Use kLa
#    values measured or fitted under that rpm, orbit, vessel, fill volume, and
#    closure; shaker rpm by itself does not specify oxygen-transfer capacity.
