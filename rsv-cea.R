# RSV Vaccincation Cost Effectivensess Analysis
#
# This analysis attempts to understand the cost effectinvess of two RSV vaccines.
# It leverages a compartmentalized SIR model vairent.
# Simulation Period: Week
#

# ==============================================================================
# Step 1: Load Libraries
# ==============================================================================

library(data.table) # Fast data manipulation and storage for simulation results
library(writexl)    # Excel file export
library(ggplot2)    # Tornado diagram visualisation
library(parallel)   # Parallel PSA via fork-based mclapply (macOS / Linux)

# ==============================================================================
# Step 2: Define variables
# ==============================================================================

# Import variable definitions
input_parameters <- read.csv('rsv_input_parameters.csv')

# ==============================================================================
# Step 2a: Initial Population Distribution
# ==============================================================================

# Total U.S. Population Age 60+
p_60 <- input_parameters$value[input_parameters$parameter == "p_60"]

# Vaccinated
v_0_p <- input_parameters$value[input_parameters$parameter == "v_0_p"] # Percentage of p_60 that represents the vaccinated population for vaccine scenario at time 0

# Infectious
i_0_p <- input_parameters$value[input_parameters$parameter == "i_0_p"] # Percentage of p_60 that represents the infectious population at time 0

# Hospitalized
h_0_p <- input_parameters$value[input_parameters$parameter == "h_0_p"] # Percentage of p_60 that represents the hospitalized population at time 0

# Suceptible
s_vn_0_p <- (1 - i_0_p - h_0_p) # Percentage of p_60 that represents the initial susceptible, vaccine-naive, no vaccine scenario at time 0
s_vn_0_v_p <- (1 - i_0_p - h_0_p - v_0_p) # Percentage of p_60 that represents the initial susceptible, vaccine-naive, vaccine scenario at time 0

# ==============================================================================
# Step 2b: Transition Probabilities
# ==============================================================================

# Vaccine Breakthrough Rate (Kappa)
kappa_arexvy <- input_parameters$value[input_parameters$parameter == "kappa_arexvy"] # Vaccine Breakthrough Rate for Arexvy. Range from 20.2% to 22.0%
kappa_abrysvo <- input_parameters$value[input_parameters$parameter == "kappa_abrysvo"] # Vaccine Breakthrough Rate for Abrysvo. Range from 21.1% to 23.3%

# RSV Vaccination Rate in US
psi <- v_0_p / 52 # Weekly vaccination rate derived from annual rate of 16.2%

# Waning immunity
gamma_v <- input_parameters$value[input_parameters$parameter == "gamma_v"] # Waning vaccine confered immunity
gamma_r <- input_parameters$value[input_parameters$parameter == "gamma_r"] # Waning natual immunity

# Conversion rate
sigma <- input_parameters$value[input_parameters$parameter == "sigma"] # Rate of conversion from exposed to infected

# Hospitalization rate
phi <- input_parameters$value[input_parameters$parameter == "phi"] # Rate of hospitalization

# Recovery rate
tau_2 <- input_parameters$value[input_parameters$parameter == "tau_2"] # Recovery Rate from hospitalization
tau_3 <- input_parameters$value[input_parameters$parameter == "tau_3"] # Recovery from sequelae

# Chronic State
rho <- input_parameters$value[input_parameters$parameter == "rho"] # Rate of transition to chronic state

# Mortality
mu <- input_parameters$value[input_parameters$parameter == "mu"] # Age adjusted weekly mortality

# RSV Mortality
delta_m <- input_parameters$value[input_parameters$parameter == "delta_m"] # Mortality rate from mild RSV
delta_s <- input_parameters$value[input_parameters$parameter == "delta_s"] # Mortality rate from severe RSV
delta_sq <- input_parameters$value[input_parameters$parameter == "delta_sq"] # Mortality rate from RSV sequelae

# ==============================================================================
# Step 2c: Health Utility For Each State
# ==============================================================================

# Each is represented as a unit of utility that must then be multiplied by the time in that state.
u_s_vn <- input_parameters$value[input_parameters$parameter == "u_s_vn"] # Susceptible, vaccine naive
u_v <- input_parameters$value[input_parameters$parameter == "u_v"] # Vaccinated (healthy, protected) - same as susceptible
u_s_ve <- input_parameters$value[input_parameters$parameter == "u_s_ve"] # Susceptible, vaccine experienced
u_e <- input_parameters$value[input_parameters$parameter == "u_e"] # Exposed
u_i <- input_parameters$value[input_parameters$parameter == "u_i"] # Infected (baseline utility for ongoing state)
u_h <- input_parameters$value[input_parameters$parameter == "u_h"] # Hospitalized (baseline utility for ongoing state)
u_sq <- input_parameters$value[input_parameters$parameter == "u_sq"] # sequalae
u_r <- input_parameters$value[input_parameters$parameter == "u_r"] # recovered
u_d <- input_parameters$value[input_parameters$parameter == "u_d"] # Death, absorbing

# Transient disutilities (applied only when entering acute states)
u_adverse_events <- input_parameters$value[input_parameters$parameter == "u_adverse_events"] # Annual utility weight decrement during adverse events (EQ-5D scale)
p_adverse_events <- input_parameters$value[input_parameters$parameter == "p_adverse_events"]  # Proportion of vaccinated people who experience an adverse event
u_disutility_infection <- input_parameters$value[input_parameters$parameter == "u_disutility_infection"] # Weekly QALY disutility during acute infection phase
u_disutility_hospitalization <- input_parameters$value[input_parameters$parameter == "u_disutility_hospitalization"] # Weekly QALY disutility during acute hospitalization phase
duration_adverse_events <- input_parameters$value[input_parameters$parameter == "duration_adverse_events"] # Duration of adverse events in weeks

# ==============================================================================
# Step 2d: Costs
# ==============================================================================

# Vaccine Costs
c_abrysvo <- input_parameters$value[input_parameters$parameter == "c_abrysvo"] # Cost of Abrysvo vaccine dose
c_arexvy <- input_parameters$value[input_parameters$parameter == "c_arexvy"] # Cost of Abrexvy dose
c_vac_admin <- input_parameters$value[input_parameters$parameter == "c_vac_admin"] # Cost to administer vaccines

# RSV Costs
c_office <- input_parameters$value[input_parameters$parameter == "c_office"] # Cost per office visit
c_ed <- input_parameters$value[input_parameters$parameter == "c_ed"] # Cost per ed day
c_hosp <- input_parameters$value[input_parameters$parameter == "c_hosp"] # Cost per hospital day
c_icu_v <- input_parameters$value[input_parameters$parameter == "c_icu_v"] # Cost per icu day with mechanincal vent
c_icu <- input_parameters$value[input_parameters$parameter == "c_icu"] # Cost per icu day without mechanical vent

# Indirect Costs
c_pl <- input_parameters$value[input_parameters$parameter == "c_pl"] # Cost per case of productivity loss
c_care <- input_parameters$value[input_parameters$parameter == "c_care"] # Cost per case of caregiver time
c_hfr <- input_parameters$value[input_parameters$parameter == "c_hfr"] # Cost of household financial risk (middle of distribution)

# Death Costs
c_death <- input_parameters$value[input_parameters$parameter == "c_death"] # Cost to medicare for patient death (middle of distribution)

# ==============================================================================
# Step 2e: Outcomes
# ==============================================================================

o_infected    <- input_parameters$value[input_parameters$parameter == "o_infected"]    # Days of infection without medical attention
o_Infected_ma <- input_parameters$value[input_parameters$parameter == "o_Infected_ma"] # Days of infection with medical attention
o_infected_h  <- input_parameters$value[input_parameters$parameter == "o_infected_h"]  # Days of infection to hospitalization
o_los_h       <- input_parameters$value[input_parameters$parameter == "o_los_h"]       # Days of hospitalization without ICU
o_los_icu     <- input_parameters$value[input_parameters$parameter == "o_los_icu"]     # Days in ICU
# o_los_h_pre_icu <- # TODO: Days in hospital prior to ICU admission (value not yet defined)
# o_los_h_icu <-    # TODO: Days in hospital after being admitted to the ICU (value not yet defined)

# ==============================================================================
# Step 2f: Discount
# ==============================================================================

p_discount_yr <- input_parameters$value[input_parameters$parameter == "p_discount_yr"]
p_discount_wk <- p_discount_yr / 52  # Weekly discount rate derived from annual rate

# ==============================================================================
# Step 2g: Simulation Controls
# ==============================================================================

max_periods <- 20000            # Safety limit to prevent infinite loops
susceptible_threshold <- 1      # Threshold below which susceptible population is considered depleted

# ==============================================================================
# Step 2h: Derived Transition Probabilities
# ==============================================================================

tau_1 <- 7/o_infected # Recovery Rate from infected

# Infection rate: derived from the RSV-ARI annual attack rate using the SIR
# final size equation. R0 = -ln(1 - A) / A, where A = i_0_p.
# beta = tau_1 * R0
attack_rate <- input_parameters$value[input_parameters$parameter == "attack_rate"]
beta <- tau_1 * (-log(1 - attack_rate) / attack_rate)

# ==============================================================================
# Step 3: Validation Functions
# ==============================================================================

# Function to perform quality checks during simulation
validate_period <- function(scenario_name, period, row_data, initial_pop, prev_deaths = 0) {

  issues <- character(0)

  # Extract populations (ensure scalar values)
  s_vn    <- row_data$s_vn_end[1]
  s_ve    <- row_data$s_ve_end[1]
  v       <- row_data$v_end[1]
  v_ae    <- row_data$v_ae_end[1]
  s_ve_ae <- row_data$s_ve_ae_end[1]
  e       <- row_data$e_end[1]
  i       <- row_data$i_end[1]
  h       <- row_data$h_end[1]
  r       <- row_data$r_end[1]
  sq      <- row_data$sq_end[1]
  d       <- row_data$d_end[1]

  # Check 1: Population Conservation
  total_pop <- s_vn + s_ve + v + v_ae + s_ve_ae + e + i + h + r + sq + d
  pop_diff <- abs(total_pop - initial_pop)
  if (pop_diff > 1) {  # Allow for small rounding errors
    issues <- c(issues, sprintf("Population not conserved: %.2f (expected %.2f, diff %.2f)",
                                total_pop, initial_pop, pop_diff))
  }

  # Check 2: Non-Negativity (with tolerance for floating-point rounding errors)
  # Allow small negative values (< -0.01) which are just rounding errors
  tolerance <- 0.01  # Tolerance for rounding errors (much less than 1 person)
  if (s_vn    < -tolerance) issues <- c(issues, sprintf("S_VN negative: %.2f", s_vn))
  if (s_ve    < -tolerance) issues <- c(issues, sprintf("S_VE negative: %.2f", s_ve))
  if (v       < -tolerance) issues <- c(issues, sprintf("V negative: %.2f", v))
  if (v_ae    < -tolerance) issues <- c(issues, sprintf("V_AE negative: %.2f", v_ae))
  if (s_ve_ae < -tolerance) issues <- c(issues, sprintf("S_VE_AE negative: %.2f", s_ve_ae))
  if (e < -tolerance) issues <- c(issues, sprintf("E negative: %.2f", e))
  if (i < -tolerance) issues <- c(issues, sprintf("I negative: %.2f", i))
  if (h < -tolerance) issues <- c(issues, sprintf("H negative: %.2f", h))
  if (r < -tolerance) issues <- c(issues, sprintf("R negative: %.2f", r))
  if (sq < -tolerance) issues <- c(issues, sprintf("SQ negative: %.2f", sq))
  if (d < -tolerance) issues <- c(issues, sprintf("D negative: %.2f", d))

  # Check 3: Monotonicity of Deaths (deaths should only increase)
  if (period > 0 && d < prev_deaths) {
    issues <- c(issues, sprintf("Deaths decreased: %.2f -> %.2f", prev_deaths, d))
  }

  # Return validation results
  if (length(issues) > 0) {
    warning(sprintf("[%s Period %d] Validation issues:\n  %s",
                    scenario_name, period, paste(issues, collapse = "\n  ")))
    return(FALSE)
  }
  return(TRUE)
}

# Function to generate comprehensive validation report
generate_validation_report <- function() {

  cat("\n")
  cat("==============================================================================\n")
  cat("                      VALIDATION REPORT\n")
  cat("==============================================================================\n")
  cat("\n")

  validation_results <- list()
  all_passed <- TRUE

  # Helper function to check a condition and record result
  check <- function(test_name, condition, details = "") {
    passed <- condition
    status <- if (passed) "PASS" else "FAIL"
    cat(sprintf("%-50s [%s] %s\n", test_name, status, details))
    if (!passed) all_passed <<- FALSE
    validation_results[[test_name]] <<- list(passed = passed, details = details)
    return(passed)
  }

  cat("--- Population Conservation ---\n")

  # Check final population for each scenario
  initial_pop <- p_60

  final_soc <- results_soc_dt[nrow(results_soc_dt)]
  total_soc <- final_soc$s_vn_end +
    final_soc$s_ve_end +
    final_soc$v_end +
    final_soc$v_ae_end +
    final_soc$s_ve_ae_end +
    final_soc$e_end +
    final_soc$i_end +
    final_soc$h_end +
    final_soc$r_end +
    final_soc$sq_end +
    final_soc$d_end
  check("SOC: Population conserved",
        abs(total_soc - initial_pop) < 100,
        sprintf("Total: %.0f, Expected: %.0f, Diff: %.0f", total_soc, initial_pop, total_soc - initial_pop))

  final_arexvy <- results_arexvy_dt[nrow(results_arexvy_dt)]
  total_arexvy <- final_arexvy$s_vn_end +
    final_arexvy$s_ve_end +
    final_arexvy$v_end +
    final_arexvy$v_ae_end +
    final_arexvy$s_ve_ae_end +
    final_arexvy$e_end +
    final_arexvy$i_end +
    final_arexvy$h_end +
    final_arexvy$r_end +
    final_arexvy$sq_end +
    final_arexvy$d_end
  check("Arexvy: Population conserved",
        abs(total_arexvy - initial_pop) < 100,
        sprintf("Total: %.0f, Expected: %.0f, Diff: %.0f", total_arexvy, initial_pop, total_arexvy - initial_pop))

  final_abrysvo <- results_abrysvo_dt[nrow(results_abrysvo_dt)]
  total_abrysvo <- final_abrysvo$s_vn_end +
    final_abrysvo$s_ve_end +
    final_abrysvo$v_end +
    final_abrysvo$v_ae_end +
    final_abrysvo$s_ve_ae_end +
    final_abrysvo$e_end +
    final_abrysvo$i_end +
    final_abrysvo$h_end +
    final_abrysvo$r_end +
    final_abrysvo$sq_end +
    final_abrysvo$d_end
  check("Abrysvo: Population conserved",
        abs(total_abrysvo - initial_pop) < 100,
        sprintf("Total: %.0f, Expected: %.0f, Diff: %.0f", total_abrysvo, initial_pop, total_abrysvo - initial_pop))

  cat("\n--- Convergence ---\n")

  # Check if simulations converged naturally
  soc_periods <- nrow(results_soc_dt) - 1
  arexvy_periods <- nrow(results_arexvy_dt) - 1
  abrysvo_periods <- nrow(results_abrysvo_dt) - 1

  check("SOC: Converged before max periods",
        soc_periods < max_periods,
        sprintf("%d periods", soc_periods))
  check("Arexvy: Converged before max periods",
        arexvy_periods < max_periods,
        sprintf("%d periods", arexvy_periods))
  check("Abrysvo: Converged before max periods",
        abrysvo_periods < max_periods,
        sprintf("%d periods", abrysvo_periods))

  cat("\n--- Monotonicity ---\n")

  # Check that costs are non-decreasing (cumulative)
  check("SOC: Costs are non-negative",
        all(results_soc_dt$cost_total >= 0),
        "")
  check("Arexvy: Costs are non-negative",
        all(results_arexvy_dt$cost_total >= 0),
        "")
  check("Abrysvo: Costs are non-negative",
        all(results_abrysvo_dt$cost_total >= 0),
        "")

  # Check that deaths only increase
  check("SOC: Deaths monotonically increasing",
        all(diff(results_soc_dt$d_end) >= 0),
        "")
  check("Arexvy: Deaths monotonically increasing",
        all(diff(results_arexvy_dt$d_end) >= 0),
        "")
  check("Abrysvo: Deaths monotonically increasing",
        all(diff(results_abrysvo_dt$d_end) >= 0),
        "")

  cat("\n--- Scenario Comparisons ---\n")

  # Vaccine scenarios should have different costs than SOC (could be higher or lower)
  arexvy_cost_diff <- summary_arexvy$total_cost - summary_soc$total_cost
  arexvy_cost_status <- ifelse(arexvy_cost_diff < 0, " (Cost-Saving)", " (Higher Cost)")
  check("Arexvy: Different total cost than SOC",
        summary_arexvy$total_cost != summary_soc$total_cost,
        sprintf("Arexvy: $%.2fM, SOC: $%.2fM, Diff: $%.2fM%s",
                summary_arexvy$total_cost / 1e6, summary_soc$total_cost / 1e6,
                arexvy_cost_diff / 1e6, arexvy_cost_status))

  abrysvo_cost_diff <- summary_abrysvo$total_cost - summary_soc$total_cost
  abrysvo_cost_status <- ifelse(abrysvo_cost_diff < 0, " (Cost-Saving)", " (Higher Cost)")
  check("Abrysvo: Different total cost than SOC",
        summary_abrysvo$total_cost != summary_soc$total_cost,
        sprintf("Abrysvo: $%.2fM, SOC: $%.2fM, Diff: $%.2fM%s",
                summary_abrysvo$total_cost / 1e6, summary_soc$total_cost / 1e6,
                abrysvo_cost_diff / 1e6, abrysvo_cost_status))

  # Vaccine scenarios should have different QALYs than SOC
  check("Arexvy: Different QALYs than SOC",
        summary_arexvy$total_qaly != summary_soc$total_qaly,
        sprintf("Diff: %.2fM QALYs", (summary_arexvy$total_qaly - summary_soc$total_qaly) / 1e6))
  check("Abrysvo: Different QALYs than SOC",
        summary_abrysvo$total_qaly != summary_soc$total_qaly,
        sprintf("Diff: %.2fM QALYs", (summary_abrysvo$total_qaly - summary_soc$total_qaly) / 1e6))

  cat("\n--- ICER Validity ---\n")

  # ICERs should be finite
  check("Arexvy vs SOC: ICER is finite",
        is.finite(icer_arexvy_vs_soc),
        sprintf("ICER: $%.2f/QALY", icer_arexvy_vs_soc))
  check("Abrysvo vs SOC: ICER is finite",
        is.finite(icer_abrysvo_vs_soc),
        sprintf("ICER: $%.2f/QALY", icer_abrysvo_vs_soc))

  cat("\n--- Internal Consistency ---\n")

  # Period 0 should have start = end for all compartments
  check("SOC: Period 0 start equals end",
        all(results_soc_dt[1, .(s_vn_start, s_ve_start, v_start, e_start, i_start,
                                h_start, r_start, sq_start, d_start)] ==
              results_soc_dt[1, .(s_vn_end, s_ve_end, v_end, e_end, i_end,
                                  h_end, r_end, sq_end, d_end)]),
        "")

  # SOC should have no vaccination
  check("SOC: No vaccination occurred",
        all(results_soc_dt$v_start == 0) && all(results_soc_dt$v_end == 0),
        "")

  cat("\n--- Extreme Value Checks ---\n")

  # Check for extreme populations (none should exceed initial population)
  max_pop_any_compartment <- max(
    max(results_soc_dt$s_vn_end, results_soc_dt$s_ve_end, results_soc_dt$v_end,
        results_soc_dt$e_end, results_soc_dt$i_end, results_soc_dt$h_end,
        results_soc_dt$r_end, results_soc_dt$sq_end),
    max(results_arexvy_dt$s_vn_end, results_arexvy_dt$s_ve_end, results_arexvy_dt$v_end,
        results_arexvy_dt$e_end, results_arexvy_dt$i_end, results_arexvy_dt$h_end,
        results_arexvy_dt$r_end, results_arexvy_dt$sq_end),
    max(results_abrysvo_dt$s_vn_end, results_abrysvo_dt$s_ve_end, results_abrysvo_dt$v_end,
        results_abrysvo_dt$e_end, results_abrysvo_dt$i_end, results_abrysvo_dt$h_end,
        results_abrysvo_dt$r_end, results_abrysvo_dt$sq_end)
  )
  check("No compartment exceeds initial population",
        max_pop_any_compartment <= initial_pop,
        sprintf("Max: %.0f, Initial: %.0f", max_pop_any_compartment, initial_pop))

  # Check for reasonable total costs
  check("SOC: Total cost is reasonable",
        summary_soc$total_cost > 0 && summary_soc$total_cost < 1e12,
        sprintf("$%.2fM", summary_soc$total_cost / 1e6))
  check("Arexvy: Total cost is reasonable",
        summary_arexvy$total_cost > 0 && summary_arexvy$total_cost < 1e12,
        sprintf("$%.2fM", summary_arexvy$total_cost / 1e6))
  check("Abrysvo: Total cost is reasonable",
        summary_abrysvo$total_cost > 0 && summary_abrysvo$total_cost < 1e12,
        sprintf("$%.2fM", summary_abrysvo$total_cost / 1e6))

  cat("\n")
  cat("==============================================================================\n")
  if (all_passed) {
    cat("                     ALL VALIDATION CHECKS PASSED\n")
  } else {
    cat("                   SOME VALIDATION CHECKS FAILED\n")
    cat("                   Please review issues above\n")
  }
  cat("==============================================================================\n")
  cat("\n")

  return(validation_results)
}

# ==============================================================================
# Step 4: Unified Simulation Function
# ==============================================================================

# run_simulation() executes a single scenario of the RSV compartmental model.
#
# Parameters:
#   p_60         - Total population age 60+
#   v_0_p        - Initial vaccinated proportion (only used when use_vaccine = TRUE)
#   i_0_p        - Initial infectious proportion
#   h_0_p        - Initial hospitalized proportion
#   beta         - Transmission rate
#   sigma        - Progression rate from exposed to infectious
#   phi          - Hospitalization rate from infectious
#   tau_1        - Recovery rate from mild infection
#   tau_2        - Recovery rate from hospitalization
#   tau_3        - Recovery rate from sequelae
#   rho          - Rate of transition to sequelae from hospitalization
#   mu           - Background mortality rate
#   delta_m      - RSV mortality rate (mild)
#   delta_s      - RSV mortality rate (severe/hospitalized)
#   delta_sq     - RSV mortality rate (sequelae)
#   gamma_r      - Waning natural immunity rate
#   psi          - Vaccination rate (only used when use_vaccine = TRUE)
#   kappa        - Vaccine breakthrough rate (only used when use_vaccine = TRUE)
#   gamma_v      - Waning vaccine-conferred immunity rate (only used when use_vaccine = TRUE)
#   u_s_vn       - Utility: susceptible, vaccine-naive
#   u_v          - Utility: vaccinated
#   u_s_ve       - Utility: susceptible, vaccine-experienced
#   u_e          - Utility: exposed
#   u_i          - Utility: infectious
#   u_h          - Utility: hospitalized
#   u_sq         - Utility: sequelae
#   u_r          - Utility: recovered
#   u_d          - Utility: deceased
#   u_adverse_events              - Weekly QALY disutility when an adverse event occurs
#   p_adverse_events              - Proportion of vaccinated people who experience an adverse event
#   u_disutility_infection        - Disutility during acute infection
#   u_disutility_hospitalization  - Disutility during acute hospitalization
#   duration_adverse_events       - Duration of adverse events in weeks
#   c_vaccine    - Cost per vaccine dose (set to 0 for SOC)
#   c_vac_admin  - Vaccine administration cost per dose
#   c_hosp       - Cost per hospital day
#   c_pl         - Cost per case of productivity loss
#   c_care       - Cost per case of caregiver time
#   c_hfr        - Cost of household financial risk per case
#   c_death      - Cost per death
#   o_infected   - Days of infection (used for disutility scaling)
#   o_los_h      - Average hospital length of stay in days
#   p_discount_wk - Weekly discount rate
#   max_periods          - Maximum number of simulation periods (safety limit)
#   susceptible_threshold - Population threshold for stopping condition
#   use_vaccine  - Logical; TRUE runs a vaccine scenario, FALSE runs standard of care
#   vaccine_name - Label for this scenario (used in validation messages and output)
#
# Returns: data.table with one row per period containing start/end compartment
#          populations, period costs, and period utilities.

run_simulation <- function(
# Population
p_60,
v_0_p,
i_0_p,
h_0_p,
# Epidemic transitions
beta,
sigma,
phi,
tau_1,
tau_2,
tau_3,
rho,
mu,
delta_m,
delta_s,
delta_sq,
gamma_r,
# Vaccine transitions (ignored when use_vaccine = FALSE)
psi,
kappa,
gamma_v,
# Utilities
u_s_vn,
u_v,
u_s_ve,
u_e,
u_i,
u_h,
u_sq,
u_r,
u_d,
u_adverse_events,
  p_adverse_events,
  u_disutility_infection,
  u_disutility_hospitalization,
  duration_adverse_events,
  # Costs
c_vaccine,
c_vac_admin,
c_hosp,
c_pl,
c_care,
c_hfr,
c_death,
# Outcomes
o_infected,
o_los_h,
# Discounting
p_discount_wk,
# Simulation controls
max_periods = 20000,
susceptible_threshold = 1,
# Scenario identification
use_vaccine = FALSE,
vaccine_name = "Standard of Care"
) {

  # ----------------------------------------------------------------------------
  # Initialize results storage
  # ----------------------------------------------------------------------------

  results_dt <- data.table(
    period = integer(),
    s_vn_start = numeric(), s_ve_start = numeric(),
    v_start = numeric(), v_ae_start = numeric(), s_ve_ae_start = numeric(),
    e_start = numeric(), i_start = numeric(), h_start = numeric(),
    r_start = numeric(), sq_start = numeric(), d_start = numeric(),
    s_vn_end = numeric(), s_ve_end = numeric(),
    v_end = numeric(), v_ae_end = numeric(), s_ve_ae_end = numeric(),
    e_end = numeric(), i_end = numeric(), h_end = numeric(),
    r_end = numeric(), sq_end = numeric(), d_end = numeric(),
    cost_total = numeric(),
    utility_total = numeric()
  )

  # ----------------------------------------------------------------------------
  # Period 0: Initial conditions
  # ----------------------------------------------------------------------------

  if (use_vaccine) {
    s_vn_init <- p_60 * (1 - i_0_p - h_0_p - v_0_p)
    v_init <- p_60 * v_0_p
  } else {
    s_vn_init <- p_60 * (1 - i_0_p - h_0_p)
    v_init <- 0
  }
  i_init <- p_60 * i_0_p
  h_init <- p_60 * h_0_p

  results_dt <- rbindlist(list(
    results_dt,
    data.table(
      period = 0L,
      s_vn_start = s_vn_init, s_ve_start = 0,
      v_start = v_init, v_ae_start = 0, s_ve_ae_start = 0,
      e_start = 0, i_start = i_init, h_start = h_init,
      r_start = 0, sq_start = 0, d_start = 0,
      s_vn_end = s_vn_init, s_ve_end = 0,
      v_end = v_init, v_ae_end = 0, s_ve_ae_end = 0,
      e_end = 0, i_end = i_init, h_end = h_init,
      r_end = 0, sq_end = 0, d_end = 0,
      cost_total = 0,
      utility_total = 0
    )
  ), use.names = TRUE)

  # ----------------------------------------------------------------------------
  # Simulation loop
  # ----------------------------------------------------------------------------

  period <- 1
  no_susceptible_count <- 0

  while (no_susceptible_count < 2 && period <= max_periods) {

    # Get previous period's ending values as current period's starting values
    prev_row <- results_dt[period]

    s_vn_t    <- prev_row$s_vn_end
    s_ve_t    <- prev_row$s_ve_end
    v_t       <- prev_row$v_end
    v_ae_t    <- prev_row$v_ae_end
    s_ve_ae_t <- prev_row$s_ve_ae_end
    e_t       <- prev_row$e_end
    i_t <- prev_row$i_end
    h_t <- prev_row$h_end
    r_t <- prev_row$r_end
    sq_t <- prev_row$sq_end
    d_t <- prev_row$d_end

    # Total living population for force of infection
    n_t <- s_vn_t +
      s_ve_t +
      v_t +
      v_ae_t +
      s_ve_ae_t +
      e_t +
      i_t +
      h_t +
      r_t +
      sq_t
    if (n_t <= 0) break

    # Force of infection
    foi <- beta * (i_t / n_t)

    # --------------------------------------------------------------------------
    # S_VN transitions
    # Vaccine scenario: vaccination (psi), breakthrough (kappa), infection (foi), death (mu)
    # SOC scenario:     infection (foi), death (mu) only
    # --------------------------------------------------------------------------

    if (use_vaccine) {
      total_exit_s_vn <- psi + kappa + foi + mu
      if (total_exit_s_vn > 0) {
        n_exit_s_vn    <- min(s_vn_t, s_vn_t * total_exit_s_vn)
        n_vaccinated   <- n_exit_s_vn * (psi / total_exit_s_vn)
        n_protected    <- n_vaccinated * (1 - kappa) # Vaccine worked
        n_unprotected  <- n_vaccinated * kappa       # Vaccine failed (breakthrough)

        # Protected: split by whether AE occurs
        s_vn_to_v      <- n_protected   * (1 - p_adverse_events)
        s_vn_to_v_ae   <- n_protected   * p_adverse_events

        # Unprotected: split by whether AE occurs
        s_vn_to_s_ve    <- n_unprotected * (1 - p_adverse_events)
        s_vn_to_s_ve_ae <- n_unprotected * p_adverse_events

        s_vn_to_e <- n_exit_s_vn * (foi / total_exit_s_vn)
        s_vn_to_d <- n_exit_s_vn * (mu  / total_exit_s_vn)
      } else {
        s_vn_to_v <- s_vn_to_v_ae <- s_vn_to_s_ve <- s_vn_to_s_ve_ae <-
          s_vn_to_e <- s_vn_to_d <- 0
      }
    } else {
      total_exit_s_vn <- foi + mu
      if (total_exit_s_vn > 0) {
        n_exit_s_vn <- min(s_vn_t, s_vn_t * total_exit_s_vn)
        s_vn_to_e <- n_exit_s_vn * (foi / total_exit_s_vn)
        s_vn_to_d <- n_exit_s_vn * (mu  / total_exit_s_vn)
      } else {
        s_vn_to_e <- s_vn_to_d <- 0
      }
      s_vn_to_v <- s_vn_to_v_ae <- s_vn_to_s_ve <- s_vn_to_s_ve_ae <- 0
    }

    # --------------------------------------------------------------------------
    # V_AE transitions: 1-period staging compartment.
    # All V_AE move to V at the start of the next period (rate = 1), with
    # a small fraction dying from background mortality.
    # Zero in SOC scenario.
    # --------------------------------------------------------------------------

    if (use_vaccine) {
      total_exit_v_ae <- 1 + mu
      if (total_exit_v_ae > 0) {
        n_exit_v_ae <- min(v_ae_t, v_ae_t * total_exit_v_ae)
        v_ae_to_v <- n_exit_v_ae * (1  / total_exit_v_ae)
        v_ae_to_d <- n_exit_v_ae * (mu / total_exit_v_ae)
      } else {
        v_ae_to_v <- v_ae_to_d <- 0
      }
    } else {
      v_ae_to_v <- v_ae_to_d <- 0
    }

    # --------------------------------------------------------------------------
    # S_VE_AE transitions: 1-period staging compartment.
    # All S_VE_AE move to S_VE at the start of the next period (rate = 1),
    # with a small fraction dying from background mortality.
    # Zero in SOC scenario.
    # --------------------------------------------------------------------------

    if (use_vaccine) {
      total_exit_s_ve_ae <- 1 + mu
      if (total_exit_s_ve_ae > 0) {
        n_exit_s_ve_ae  <- min(s_ve_ae_t, s_ve_ae_t * total_exit_s_ve_ae)
        s_ve_ae_to_s_ve <- n_exit_s_ve_ae * (1  / total_exit_s_ve_ae)
        s_ve_ae_to_d    <- n_exit_s_ve_ae * (mu / total_exit_s_ve_ae)
      } else {
        s_ve_ae_to_s_ve <- s_ve_ae_to_d <- 0
      }
    } else {
      s_ve_ae_to_s_ve <- s_ve_ae_to_d <- 0
    }

    # --------------------------------------------------------------------------
    # V transitions: waning immunity (gamma_v) to S_VE, death (mu)
    # Zero in SOC scenario.
    # --------------------------------------------------------------------------

    if (use_vaccine) {
      total_exit_v <- gamma_v + mu
      if (total_exit_v > 0) {
        n_exit_v <- min(v_t, v_t * total_exit_v)
        v_to_s_ve <- n_exit_v * (gamma_v / total_exit_v)
        v_to_d <- n_exit_v * (mu / total_exit_v)
      } else {
        v_to_s_ve <- v_to_d <- 0
      }
    } else {
      v_to_s_ve <- v_to_d <- 0
    }

    # --------------------------------------------------------------------------
    # S_VE transitions: infection (foi), death (mu)
    # Zero in SOC scenario.
    # --------------------------------------------------------------------------

    if (use_vaccine) {
      total_exit_s_ve <- foi + mu
      if (total_exit_s_ve > 0) {
        n_exit_s_ve <- min(s_ve_t, s_ve_t * total_exit_s_ve)
        s_ve_to_e <- n_exit_s_ve * (foi / total_exit_s_ve)
        s_ve_to_d <- n_exit_s_ve * (mu / total_exit_s_ve)
      } else {
        s_ve_to_e <- s_ve_to_d <- 0
      }
    } else {
      s_ve_to_e <- s_ve_to_d <- 0
    }

    # --------------------------------------------------------------------------
    # E transitions: progression to infectious (sigma), death (mu)
    # --------------------------------------------------------------------------

    total_exit_e <- sigma + mu
    if (total_exit_e > 0) {
      n_exit_e <- min(e_t, e_t * total_exit_e)
      e_to_i <- n_exit_e * (sigma / total_exit_e)
      e_to_d <- n_exit_e * (mu / total_exit_e)
    } else {
      e_to_i <- e_to_d <- 0
    }

    # --------------------------------------------------------------------------
    # I transitions: hospitalization (phi), recovery (tau_1), death (mu + delta_m)
    # --------------------------------------------------------------------------

    total_exit_i <- phi + tau_1 + mu + delta_m
    if (total_exit_i > 0) {
      n_exit_i <- min(i_t, i_t * total_exit_i)
      i_to_h <- n_exit_i * (phi / total_exit_i)
      i_to_r <- n_exit_i * (tau_1 / total_exit_i)
      i_to_d <- n_exit_i * ((mu + delta_m) / total_exit_i)
    } else {
      i_to_h <- i_to_r <- i_to_d <- 0
    }

    # --------------------------------------------------------------------------
    # H transitions: sequelae (rho), recovery (tau_2), death (mu + delta_s)
    # --------------------------------------------------------------------------

    total_exit_h <- rho + tau_2 + mu + delta_s
    if (total_exit_h > 0) {
      n_exit_h <- min(h_t, h_t * total_exit_h)
      h_to_sq <- n_exit_h * (rho / total_exit_h)
      h_to_r <- n_exit_h * (tau_2 / total_exit_h)
      h_to_d <- n_exit_h * ((mu + delta_s) / total_exit_h)
    } else {
      h_to_sq <- h_to_r <- h_to_d <- 0
    }

    # --------------------------------------------------------------------------
    # SQ transitions: recovery (tau_3), death (mu + delta_sq)
    # --------------------------------------------------------------------------

    total_exit_sq <- tau_3 + mu + delta_sq
    if (total_exit_sq > 0) {
      n_exit_sq <- min(sq_t, sq_t * total_exit_sq)
      sq_to_r <- n_exit_sq * (tau_3 / total_exit_sq)
      sq_to_d <- n_exit_sq * ((mu + delta_sq) / total_exit_sq)
    } else {
      sq_to_r <- sq_to_d <- 0
    }

    # --------------------------------------------------------------------------
    # R transitions: waning natural immunity (gamma_r) to S_VN, death (mu)
    # --------------------------------------------------------------------------

    total_exit_r <- gamma_r + mu
    if (total_exit_r > 0) {
      n_exit_r <- min(r_t, r_t * total_exit_r)
      r_to_s_vn <- n_exit_r * (gamma_r / total_exit_r)
      r_to_d <- n_exit_r * (mu / total_exit_r)
    } else {
      r_to_s_vn <- r_to_d <- 0
    }

    # --------------------------------------------------------------------------
    # End-of-period populations
    # --------------------------------------------------------------------------

    s_vn_end <- s_vn_t -
      s_vn_to_v - s_vn_to_v_ae -
      s_vn_to_s_ve - s_vn_to_s_ve_ae -
      s_vn_to_e - s_vn_to_d + r_to_s_vn

    # V_AE and S_VE_AE: receive new entrants, all prior residents exit
    v_ae_end    <- v_ae_t    - v_ae_to_v    - v_ae_to_d    + s_vn_to_v_ae
    s_ve_ae_end <- s_ve_ae_t - s_ve_ae_to_s_ve - s_ve_ae_to_d + s_vn_to_s_ve_ae

    # V: receives new direct vaccinations + graduates from V_AE
    v_end <- v_t + s_vn_to_v + v_ae_to_v - v_to_s_ve - v_to_d

    # S_VE: receives vaccine failures + V_AE graduates + S_VE_AE graduates + V waning
    s_ve_end <- s_ve_t +
      s_vn_to_s_ve + s_ve_ae_to_s_ve + v_to_s_ve -
      s_ve_to_e - s_ve_to_d

    e_end  <- e_t  + s_vn_to_e + s_ve_to_e - e_to_i - e_to_d
    i_end  <- i_t  + e_to_i - i_to_h - i_to_r - i_to_d
    h_end  <- h_t  + i_to_h - h_to_sq - h_to_r - h_to_d
    r_end  <- r_t  + i_to_r + h_to_r + sq_to_r - r_to_s_vn - r_to_d
    sq_end <- sq_t + h_to_sq - sq_to_r - sq_to_d
    d_end  <- d_t  +
      s_vn_to_d + s_ve_to_d + v_to_d + v_ae_to_d + s_ve_ae_to_d +
      e_to_d + i_to_d + h_to_d + sq_to_d + r_to_d

    # --------------------------------------------------------------------------
    # Costs
    # --------------------------------------------------------------------------

    # Vaccination costs: incurred for all who receive the vaccine this period
    # (V, V_AE, S_VE, S_VE_AE — every vaccinated person pays the dose + admin cost)
    cost_vaccine_period <- if (use_vaccine)
      (s_vn_to_v + s_vn_to_v_ae + s_vn_to_s_ve + s_vn_to_s_ve_ae) *
        (c_vaccine + c_vac_admin) else 0

    # Hospitalization costs: incurred when entering H compartment
    cost_hosp <- i_to_h * o_los_h * c_hosp

    # Indirect costs: incurred when entering I compartment (per case)
    cost_indirect <- e_to_i * (c_pl + c_care + c_hfr)

    # Death costs: incurred when entering D compartment (one-time per death)
    cost_death_period <- (s_vn_to_d +
      s_ve_to_d +
      v_to_d + v_ae_to_d + s_ve_ae_to_d +
      e_to_d +
      i_to_d +
      h_to_d +
      sq_to_d +
      r_to_d) * c_death

    cost_total <- cost_vaccine_period +
      cost_hosp +
      cost_indirect +
      cost_death_period

    # --------------------------------------------------------------------------
    # Utilities
    # --------------------------------------------------------------------------

    # Ongoing state utilities (people currently in each state).
    # All utility weights are annual and divided by 52 to get weekly QALYs.
    # V_AE and S_VE_AE use (u_v + u_adverse_events) and
    # (u_s_ve + u_adverse_events) respectively, capturing the AE disutility
    # for exactly the one period these people spend in the staging compartment.
    # Infection and hospitalization disutilities are weekly values applied
    # outside the /52, scaled by the actual fraction of the period spent sick.
    utility_states <- (s_vn_t    * u_s_vn +
      s_ve_t    * u_s_ve +
      v_t       * u_v +
      v_ae_t    * (u_v   + u_adverse_events) +
      s_ve_ae_t * (u_s_ve + u_adverse_events) +
      e_t    * u_e +
      i_t    * u_i +
      h_t    * u_h +
      r_t    * u_r +
      sq_t   * u_sq +
      d_t    * u_d) / 52 +
      i_t * u_disutility_infection       * (o_infected / 7) +
      h_t * u_disutility_hospitalization * (o_los_h    / 7)

    utility_total <- utility_states

    # --------------------------------------------------------------------------
    # Apply discounting
    # --------------------------------------------------------------------------

    discount_factor <- 1 / (1 + p_discount_wk)^period
    cost_total <- cost_total * discount_factor
    utility_total <- utility_total * discount_factor

    # --------------------------------------------------------------------------
    # Store results for current period
    # --------------------------------------------------------------------------

    results_dt <- rbindlist(list(
      results_dt,
      data.table(
        period = period,
        s_vn_start = s_vn_t, s_ve_start = s_ve_t,
        v_start = v_t, v_ae_start = v_ae_t, s_ve_ae_start = s_ve_ae_t,
        e_start = e_t, i_start = i_t, h_start = h_t,
        r_start = r_t, sq_start = sq_t, d_start = d_t,
        s_vn_end = s_vn_end, s_ve_end = s_ve_end,
        v_end = v_end, v_ae_end = v_ae_end, s_ve_ae_end = s_ve_ae_end,
        e_end = e_end, i_end = i_end, h_end = h_end,
        r_end = r_end, sq_end = sq_end, d_end = d_end,
        cost_total = cost_total,
        utility_total = utility_total
      )
    ), use.names = TRUE)

    # --------------------------------------------------------------------------
    # Real-time validation
    # --------------------------------------------------------------------------

    validate_period(vaccine_name, period, results_dt[.N], p_60,
                    if (period > 0) results_dt[.N - 1]$d_end[1] else 0)

    # --------------------------------------------------------------------------
    # Stopping condition
    # --------------------------------------------------------------------------

    # Vaccine scenario: require both susceptible pools to be depleted
    # SOC scenario:     only the vaccine-naive pool exists
    if (use_vaccine) {
      if (s_vn_end < susceptible_threshold && s_ve_end < susceptible_threshold) {
        no_susceptible_count <- no_susceptible_count + 1
      } else {
        no_susceptible_count <- 0
      }
    } else {
      if (s_vn_end < susceptible_threshold) {
        no_susceptible_count <- no_susceptible_count + 1
      } else {
        no_susceptible_count <- 0
      }
    }

    period <- period + 1
  }

  cat(vaccine_name, "simulation completed after", period - 1, "periods\n")
  cat("Final population distribution:\n")
  cat("  Deceased:", d_end, "\n")
  cat("  Living:", n_t - (d_end - d_t), "\n")

  return(results_dt)
}

# ==============================================================================
# Step 5: Run Simulations
# ==============================================================================

# Standard of Care (no vaccine)
results_soc_dt <- run_simulation(
  p_60 = p_60, v_0_p = v_0_p, i_0_p = i_0_p, h_0_p = h_0_p,
  beta = beta, sigma = sigma, phi = phi,
  tau_1 = tau_1, tau_2 = tau_2, tau_3 = tau_3,
  rho = rho, mu = mu, delta_m = delta_m, delta_s = delta_s, delta_sq = delta_sq,
  gamma_r = gamma_r,
  psi = psi, kappa = 0, gamma_v = gamma_v,
  u_s_vn = u_s_vn, u_v = u_v, u_s_ve = u_s_ve,
  u_e = u_e, u_i = u_i, u_h = u_h, u_sq = u_sq, u_r = u_r, u_d = u_d,
  u_adverse_events = u_adverse_events,
  p_adverse_events = p_adverse_events,
  u_disutility_infection = u_disutility_infection,
  u_disutility_hospitalization = u_disutility_hospitalization,
  duration_adverse_events = duration_adverse_events,
  c_vaccine = 0, c_vac_admin = c_vac_admin, c_hosp = c_hosp,
  c_pl = c_pl, c_care = c_care, c_hfr = c_hfr, c_death = c_death,
  o_infected = o_infected, o_los_h = o_los_h,
  p_discount_wk = p_discount_wk,
  max_periods = max_periods, susceptible_threshold = susceptible_threshold,
  use_vaccine = FALSE,
  vaccine_name = "Standard of Care"
)

# Arexvy vaccine
results_arexvy_dt <- run_simulation(
  p_60 = p_60, v_0_p = v_0_p, i_0_p = i_0_p, h_0_p = h_0_p,
  beta = beta, sigma = sigma, phi = phi,
  tau_1 = tau_1, tau_2 = tau_2, tau_3 = tau_3,
  rho = rho, mu = mu, delta_m = delta_m, delta_s = delta_s, delta_sq = delta_sq,
  gamma_r = gamma_r,
  psi = psi, kappa = kappa_arexvy, gamma_v = gamma_v,
  u_s_vn = u_s_vn, u_v = u_v, u_s_ve = u_s_ve,
  u_e = u_e, u_i = u_i, u_h = u_h, u_sq = u_sq, u_r = u_r, u_d = u_d,
  u_adverse_events = u_adverse_events,
  p_adverse_events = p_adverse_events,
  u_disutility_infection = u_disutility_infection,
  u_disutility_hospitalization = u_disutility_hospitalization,
  duration_adverse_events = duration_adverse_events,
  c_vaccine = c_arexvy, c_vac_admin = c_vac_admin, c_hosp = c_hosp,
  c_pl = c_pl, c_care = c_care, c_hfr = c_hfr, c_death = c_death,
  o_infected = o_infected, o_los_h = o_los_h,
  p_discount_wk = p_discount_wk,
  max_periods = max_periods, susceptible_threshold = susceptible_threshold,
  use_vaccine = TRUE,
  vaccine_name = "Arexvy"
)

# Abrysvo vaccine
results_abrysvo_dt <- run_simulation(
  p_60 = p_60, v_0_p = v_0_p, i_0_p = i_0_p, h_0_p = h_0_p,
  beta = beta, sigma = sigma, phi = phi,
  tau_1 = tau_1, tau_2 = tau_2, tau_3 = tau_3,
  rho = rho, mu = mu, delta_m = delta_m, delta_s = delta_s, delta_sq = delta_sq,
  gamma_r = gamma_r,
  psi = psi, kappa = kappa_abrysvo, gamma_v = gamma_v,
  u_s_vn = u_s_vn, u_v = u_v, u_s_ve = u_s_ve,
  u_e = u_e, u_i = u_i, u_h = u_h, u_sq = u_sq, u_r = u_r, u_d = u_d,
  u_adverse_events = u_adverse_events,
  p_adverse_events = p_adverse_events,
  u_disutility_infection = u_disutility_infection,
  u_disutility_hospitalization = u_disutility_hospitalization,
  duration_adverse_events = duration_adverse_events,
  c_vaccine = c_abrysvo, c_vac_admin = c_vac_admin, c_hosp = c_hosp,
  c_pl = c_pl, c_care = c_care, c_hfr = c_hfr, c_death = c_death,
  o_infected = o_infected, o_los_h = o_los_h,
  p_discount_wk = p_discount_wk,
  max_periods = max_periods, susceptible_threshold = susceptible_threshold,
  use_vaccine = TRUE,
  vaccine_name = "Abrysvo"
)

# ==============================================================================
# Step 6: Aggregate and Present Results
# ==============================================================================

# Calculate summary statistics for each scenario
summary_soc <- data.table(
  scenario = "Standard of Care",
  periods = nrow(results_soc_dt) - 1,        # Exclude period 0
  total_cost = sum(results_soc_dt$cost_total),
  total_qaly = sum(results_soc_dt$utility_total)
)

summary_arexvy <- data.table(
  scenario = "Arexvy",
  periods = nrow(results_arexvy_dt) - 1,     # Exclude period 0
  total_cost = sum(results_arexvy_dt$cost_total),
  total_qaly = sum(results_arexvy_dt$utility_total)
)

summary_abrysvo <- data.table(
  scenario = "Abrysvo",
  periods = nrow(results_abrysvo_dt) - 1,    # Exclude period 0
  total_cost = sum(results_abrysvo_dt$cost_total),
  total_qaly = sum(results_abrysvo_dt$utility_total)
)

# Calculate ICER for both vaccines vs SOC (for validation and display)
inc_cost_arexvy_vs_soc <- summary_arexvy$total_cost - summary_soc$total_cost
inc_qaly_arexvy_vs_soc <- summary_arexvy$total_qaly - summary_soc$total_qaly
icer_arexvy_vs_soc <- inc_cost_arexvy_vs_soc / inc_qaly_arexvy_vs_soc

inc_cost_abrysvo_vs_soc <- summary_abrysvo$total_cost - summary_soc$total_cost
inc_qaly_abrysvo_vs_soc <- summary_abrysvo$total_qaly - summary_soc$total_qaly
icer_abrysvo_vs_soc <- inc_cost_abrysvo_vs_soc / inc_qaly_abrysvo_vs_soc

# Determine row ordering by TOTAL QALY ascending (least to most effective).
# This is the standard CEA ordering required for correct dominance analysis.
# Row 2: Vaccine with fewer total QALYs (less effective)
# Row 3: Vaccine with more total QALYs (more effective)
if (summary_arexvy$total_qaly < summary_abrysvo$total_qaly) {
  row2_summary <- summary_arexvy
  row2_inc_cost <- inc_cost_arexvy_vs_soc
  row2_inc_qaly <- inc_qaly_arexvy_vs_soc
  row2_icer <- icer_arexvy_vs_soc
  row2_dominated <- row2_inc_cost > 0 && row2_inc_qaly < 0  # Costs more, fewer QALYs

  row3_summary <- summary_abrysvo
} else {
  row2_summary <- summary_abrysvo
  row2_inc_cost <- inc_cost_abrysvo_vs_soc
  row2_inc_qaly <- inc_qaly_abrysvo_vs_soc
  row2_icer <- icer_abrysvo_vs_soc
  row2_dominated <- row2_inc_cost > 0 && row2_inc_qaly < 0  # Costs more, fewer QALYs

  row3_summary <- summary_arexvy
}

# Calculate incrementals for row 3 (compared to row 2)
row3_inc_cost <- row3_summary$total_cost - row2_summary$total_cost
row3_inc_qaly <- row3_summary$total_qaly - row2_summary$total_qaly
row3_icer <- row3_inc_cost / row3_inc_qaly
row3_dominated <- row3_inc_cost > 0 && row3_inc_qaly < 0  # Costs more, fewer QALYs

# Helper function to format currency with $ sign, commas, and 2 decimal places
format_currency <- function(x) {
  sprintf("$%s", formatC(x, format = "f", digits = 2, big.mark = ","))
}

# Helper function to format numbers with commas and 2 decimal places
format_number <- function(x) {
  formatC(x, format = "f", digits = 2, big.mark = ",")
}

# Helper function to format ICER with dominated notation if applicable
format_icer <- function(icer, dominated) {
  if (dominated) {
    paste0(format_currency(icer), " (Dominated)")
  } else {
    format_currency(icer)
  }
}

# Build final results table with proper formatting
# Convert to millions for readability
results_summary <- data.table(
  Scenario = c(
    summary_soc$scenario,
    row2_summary$scenario,
    row3_summary$scenario
  ),
  Total_Periods = c(
    summary_soc$periods,
    row2_summary$periods,
    row3_summary$periods
  ),
  `Total Cost (Millions)` = c(
    format_currency(summary_soc$total_cost / 1e6),
    format_currency(row2_summary$total_cost / 1e6),
    format_currency(row3_summary$total_cost / 1e6)
  ),
  `Total QALY (Millions)` = c(
    format_number(summary_soc$total_qaly / 1e6),
    format_number(row2_summary$total_qaly / 1e6),
    format_number(row3_summary$total_qaly / 1e6)
  ),
  `Incremental Cost (Millions)` = c(
    "-",
    format_currency(row2_inc_cost / 1e6),
    format_currency(row3_inc_cost / 1e6)
  ),
  `Incremental QALY (Millions)` = c(
    "-",
    format_number(row2_inc_qaly / 1e6),
    format_number(row3_inc_qaly / 1e6)
  ),
  `ICER ($/QALY)` = c(
    "-",
    format_icer(row2_icer, row2_dominated),
    format_icer(row3_icer, row3_dominated)
  )
)

# Export results to Excel
write_xlsx(results_summary, "rsv_cea_results.xlsx")
cat("Results exported to: rsv_cea_results.xlsx\n")

# Print formatted results table to console
cat("\n")
cat("==============================================================================\n")
cat("                   Cost-Effectiveness Analysis Results\n")
cat("==============================================================================\n")
cat("\n")
print(results_summary)
cat("\n")
cat("Note: ICER = Incremental Cost-Effectiveness Ratio (cost per QALY gained)\n")
cat("      Rows ordered by total QALYs ascending (least to most effective)\n")
cat("      Row 2 compared to Standard of Care\n")
cat("      Row 3 compared to Row 2\n")
cat("\n")

# ==============================================================================
# Step 7: Deterministic One-Way Sensitivity Analysis (OWSA)
# ==============================================================================

# Only parameters that (a) have both bounds defined in the CSV and (b) map
# directly to an argument of run_simulation() are eligible.  o_Infected_ma has
# bounds in the CSV but is not used by run_simulation(), so it is excluded.
owsa_eligible_params <- c(
  "kappa_arexvy", "kappa_abrysvo",
  "u_s_vn", "u_v", "u_s_ve", "u_e", "u_i", "u_h",
  "u_sq", "u_r", "u_adverse_events",
  "u_disutility_infection", "u_disutility_hospitalization",
  "c_hfr", "c_death", "o_infected"
)

# Willingness-to-pay threshold used to compute Net Monetary Benefit (NMB).
# NMB = WTP * delta_QALY - delta_Cost is always continuous, so it is used
# for the tornado diagram instead of ICER.  ICER can be discontinuous when
# the incremental QALY denominator passes through zero across the parameter
# range (which it does for several utility values including u_s_vn), which
# causes bars that do not bracket the baseline ICER.  NMB avoids this pathology entirely.
wtp_threshold <- 100000  # USD per QALY (standard US benchmark)

# Build table of eligible parameters that actually have numeric bounds
owsa_param_table <- input_parameters[
  input_parameters$parameter %in% owsa_eligible_params,
]
owsa_param_table$lower_bound <- suppressWarnings(
  as.numeric(owsa_param_table$lower_bound)
)
owsa_param_table$upper_bound <- suppressWarnings(
  as.numeric(owsa_param_table$upper_bound)
)
owsa_param_table <- owsa_param_table[
  !is.na(owsa_param_table$lower_bound) &
    !is.na(owsa_param_table$upper_bound),
]

# Helper: run all three scenarios with one parameter overridden.
# Console output is redirected to a temp file via sink() so that assignments
# are made directly in this function's frame (avoiding the scoping ambiguity
# of capture.output() when called inside suppressWarnings()).
run_owsa_scenario <- function(param_name, param_value) {

  # tau_1 and beta are both derived from o_infected (tau_1 = 7 / o_infected).
  # When o_infected is the parameter being varied, propagate the change through
  # to tau_1 and beta so the epidemic dynamics are consistent.
  local_o_infected <- if (param_name == "o_infected") param_value else o_infected
  local_tau_1      <- 7 / local_o_infected
  local_beta       <- local_tau_1 * (-log(1 - attack_rate) / attack_rate)

  # Return overridden value for the target parameter, baseline otherwise
  pv <- function(pname, base_val) {
    if (param_name == pname) param_value else base_val
  }

  # Redirect stdout to suppress verbose simulation logs; on.exit ensures the
  # sink is always restored even if an error occurs mid-simulation.
  sink(tempfile())
  on.exit(sink(), add = TRUE)

  res_soc <- suppressWarnings(run_simulation(
    p_60 = p_60, v_0_p = v_0_p, i_0_p = i_0_p, h_0_p = h_0_p,
    beta = local_beta, sigma = sigma, phi = phi,
    tau_1 = local_tau_1, tau_2 = tau_2, tau_3 = tau_3,
    rho = rho, mu = mu,
    delta_m = delta_m, delta_s = delta_s, delta_sq = delta_sq,
    gamma_r = gamma_r, psi = psi, kappa = 0, gamma_v = gamma_v,
    u_s_vn = pv("u_s_vn", u_s_vn), u_v = pv("u_v", u_v),
    u_s_ve = pv("u_s_ve", u_s_ve), u_e = pv("u_e", u_e),
    u_i = pv("u_i", u_i), u_h = pv("u_h", u_h),
    u_sq = pv("u_sq", u_sq), u_r = pv("u_r", u_r), u_d = u_d,
    u_adverse_events = pv("u_adverse_events", u_adverse_events),
    p_adverse_events = p_adverse_events,
    u_disutility_infection = pv(
      "u_disutility_infection", u_disutility_infection
    ),
    u_disutility_hospitalization = pv(
      "u_disutility_hospitalization", u_disutility_hospitalization
    ),
    duration_adverse_events = duration_adverse_events,
    c_vaccine = 0, c_vac_admin = c_vac_admin, c_hosp = c_hosp,
    c_pl = c_pl, c_care = c_care,
    c_hfr = pv("c_hfr", c_hfr), c_death = pv("c_death", c_death),
    o_infected = pv("o_infected", o_infected), o_los_h = o_los_h,
    p_discount_wk = p_discount_wk,
    max_periods = max_periods,
    susceptible_threshold = susceptible_threshold,
    use_vaccine = FALSE, vaccine_name = "SOC_OWSA"
  ))

  res_arexvy <- suppressWarnings(run_simulation(
    p_60 = p_60, v_0_p = v_0_p, i_0_p = i_0_p, h_0_p = h_0_p,
    beta = local_beta, sigma = sigma, phi = phi,
    tau_1 = local_tau_1, tau_2 = tau_2, tau_3 = tau_3,
    rho = rho, mu = mu,
    delta_m = delta_m, delta_s = delta_s, delta_sq = delta_sq,
    gamma_r = gamma_r, psi = psi,
    kappa = pv("kappa_arexvy", kappa_arexvy), gamma_v = gamma_v,
    u_s_vn = pv("u_s_vn", u_s_vn), u_v = pv("u_v", u_v),
    u_s_ve = pv("u_s_ve", u_s_ve), u_e = pv("u_e", u_e),
    u_i = pv("u_i", u_i), u_h = pv("u_h", u_h),
    u_sq = pv("u_sq", u_sq), u_r = pv("u_r", u_r), u_d = u_d,
    u_adverse_events = pv("u_adverse_events", u_adverse_events),
    p_adverse_events = p_adverse_events,
    u_disutility_infection = pv(
      "u_disutility_infection", u_disutility_infection
    ),
    u_disutility_hospitalization = pv(
      "u_disutility_hospitalization", u_disutility_hospitalization
    ),
    duration_adverse_events = duration_adverse_events,
    c_vaccine = c_arexvy, c_vac_admin = c_vac_admin, c_hosp = c_hosp,
    c_pl = c_pl, c_care = c_care,
    c_hfr = pv("c_hfr", c_hfr), c_death = pv("c_death", c_death),
    o_infected = pv("o_infected", o_infected), o_los_h = o_los_h,
    p_discount_wk = p_discount_wk,
    max_periods = max_periods,
    susceptible_threshold = susceptible_threshold,
    use_vaccine = TRUE, vaccine_name = "Arexvy_OWSA"
  ))

  res_abrysvo <- suppressWarnings(run_simulation(
    p_60 = p_60, v_0_p = v_0_p, i_0_p = i_0_p, h_0_p = h_0_p,
    beta = local_beta, sigma = sigma, phi = phi,
    tau_1 = local_tau_1, tau_2 = tau_2, tau_3 = tau_3,
    rho = rho, mu = mu,
    delta_m = delta_m, delta_s = delta_s, delta_sq = delta_sq,
    gamma_r = gamma_r, psi = psi,
    kappa = pv("kappa_abrysvo", kappa_abrysvo), gamma_v = gamma_v,
    u_s_vn = pv("u_s_vn", u_s_vn), u_v = pv("u_v", u_v),
    u_s_ve = pv("u_s_ve", u_s_ve), u_e = pv("u_e", u_e),
    u_i = pv("u_i", u_i), u_h = pv("u_h", u_h),
    u_sq = pv("u_sq", u_sq), u_r = pv("u_r", u_r), u_d = u_d,
    u_adverse_events = pv("u_adverse_events", u_adverse_events),
    p_adverse_events = p_adverse_events,
    u_disutility_infection = pv(
      "u_disutility_infection", u_disutility_infection
    ),
    u_disutility_hospitalization = pv(
      "u_disutility_hospitalization", u_disutility_hospitalization
    ),
    duration_adverse_events = duration_adverse_events,
    c_vaccine = c_abrysvo, c_vac_admin = c_vac_admin, c_hosp = c_hosp,
    c_pl = c_pl, c_care = c_care,
    c_hfr = pv("c_hfr", c_hfr), c_death = pv("c_death", c_death),
    o_infected = pv("o_infected", o_infected), o_los_h = o_los_h,
    p_discount_wk = p_discount_wk,
    max_periods = max_periods,
    susceptible_threshold = susceptible_threshold,
    use_vaccine = TRUE, vaccine_name = "Abrysvo_OWSA"
  ))

  cost_soc     <- sum(res_soc$cost_total)
  qaly_soc     <- sum(res_soc$utility_total)
  cost_arexvy  <- sum(res_arexvy$cost_total)
  qaly_arexvy  <- sum(res_arexvy$utility_total)
  cost_abrysvo <- sum(res_abrysvo$cost_total)
  qaly_abrysvo <- sum(res_abrysvo$utility_total)

  inc_cost_arexvy  <- cost_arexvy  - cost_soc
  inc_cost_abrysvo <- cost_abrysvo - cost_soc
  inc_qaly_arexvy  <- qaly_arexvy  - qaly_soc
  inc_qaly_abrysvo <- qaly_abrysvo - qaly_soc

  list(
    icer_arexvy  = inc_cost_arexvy  / inc_qaly_arexvy,
    icer_abrysvo = inc_cost_abrysvo / inc_qaly_abrysvo,
    # NMB = WTP * delta_QALY - delta_Cost (always continuous)
    nmb_arexvy   = wtp_threshold * inc_qaly_arexvy  - inc_cost_arexvy,
    nmb_abrysvo  = wtp_threshold * inc_qaly_abrysvo - inc_cost_abrysvo
  )
}

# Baseline NMB values (computed from the Step 6 incremental cost/QALY)
nmb_baseline_arexvy  <- wtp_threshold * inc_qaly_arexvy_vs_soc  -
  inc_cost_arexvy_vs_soc
nmb_baseline_abrysvo <- wtp_threshold * inc_qaly_abrysvo_vs_soc -
  inc_cost_abrysvo_vs_soc

# Iterate over each eligible parameter, testing at lower and upper bound
cat("Running One-Way Sensitivity Analysis...\n")
cat(sprintf("  Testing %d parameters\n\n", nrow(owsa_param_table)))

owsa_results <- data.table(
  parameter             = character(),
  description           = character(),
  base_value            = numeric(),
  lower_bound           = numeric(),
  upper_bound           = numeric(),
  icer_arexvy_at_lower  = numeric(),
  icer_arexvy_at_upper  = numeric(),
  icer_abrysvo_at_lower = numeric(),
  icer_abrysvo_at_upper = numeric(),
  nmb_arexvy_at_lower   = numeric(),
  nmb_arexvy_at_upper   = numeric(),
  nmb_abrysvo_at_lower  = numeric(),
  nmb_abrysvo_at_upper  = numeric()
)

for (i in seq_len(nrow(owsa_param_table))) {
  pname   <- owsa_param_table$parameter[i]
  p_lower <- owsa_param_table$lower_bound[i]
  p_upper <- owsa_param_table$upper_bound[i]
  pbase   <- as.numeric(owsa_param_table$value[i])
  pdesc   <- owsa_param_table$description[i]

  cat(sprintf("  [%d/%d] %s\n", i, nrow(owsa_param_table), pname))

  res_lower <- run_owsa_scenario(pname, p_lower)
  res_upper <- run_owsa_scenario(pname, p_upper)

  owsa_results <- rbindlist(list(
    owsa_results,
    data.table(
      parameter             = pname,
      description           = pdesc,
      base_value            = pbase,
      lower_bound           = p_lower,
      upper_bound           = p_upper,
      icer_arexvy_at_lower  = res_lower$icer_arexvy,
      icer_arexvy_at_upper  = res_upper$icer_arexvy,
      icer_abrysvo_at_lower = res_lower$icer_abrysvo,
      icer_abrysvo_at_upper = res_upper$icer_abrysvo,
      nmb_arexvy_at_lower   = res_lower$nmb_arexvy,
      nmb_arexvy_at_upper   = res_upper$nmb_arexvy,
      nmb_abrysvo_at_lower  = res_lower$nmb_abrysvo,
      nmb_abrysvo_at_upper  = res_upper$nmb_abrysvo
    )
  ))
}

cat("\n")

# Print OWSA summary tables
cat("==============================================================================\n")
cat("        One-Way Sensitivity Analysis Results (sorted by ICER range)\n")
cat("==============================================================================\n\n")

owsa_arexvy_display <- owsa_results[, .(
  Parameter       = parameter,
  `Base Value`    = round(base_value, 4),
  `Lower Bound`   = round(lower_bound, 4),
  `ICER at Lower` = round(icer_arexvy_at_lower, 0),
  `Upper Bound`   = round(upper_bound, 4),
  `ICER at Upper` = round(icer_arexvy_at_upper, 0),
  `ICER Range`    = round(
    abs(icer_arexvy_at_upper - icer_arexvy_at_lower), 0
  )
)]

owsa_abrysvo_display <- owsa_results[, .(
  Parameter       = parameter,
  `Base Value`    = round(base_value, 4),
  `Lower Bound`   = round(lower_bound, 4),
  `ICER at Lower` = round(icer_abrysvo_at_lower, 0),
  `Upper Bound`   = round(upper_bound, 4),
  `ICER at Upper` = round(icer_abrysvo_at_upper, 0),
  `ICER Range`    = round(
    abs(icer_abrysvo_at_upper - icer_abrysvo_at_lower), 0
  )
)]

cat("Arexvy vs Standard of Care:\n")
print(owsa_arexvy_display[order(-`ICER Range`)])
cat("\n")
cat("Abrysvo vs Standard of Care:\n")
print(owsa_abrysvo_display[order(-`ICER Range`)])
cat("\n")

# Build an ICER tornado diagram.
# Parameters where the baseline ICER falls outside [lower-bound ICER,
# upper-bound ICER] are excluded and listed separately.  This happens when
# the incremental QALY denominator changes sign across the parameter range
# (a mathematical discontinuity in the ICER formula, not a model error).
make_tornado_plot <- function(owsa_dt, col_lower, col_upper,
                              icer_baseline, title_str) {

  plot_dt <- data.table(
    label      = owsa_dt$description,
    icer_lower = owsa_dt[[col_lower]],
    icer_upper = owsa_dt[[col_upper]]
  )

  # Replace non-finite ICERs (0/0) with the baseline so the bar collapses
  # to zero width and is removed by the range filter below
  plot_dt[!is.finite(icer_lower), icer_lower := icer_baseline]
  plot_dt[!is.finite(icer_upper), icer_upper := icer_baseline]

  plot_dt[, bar_min   := pmin(icer_lower, icer_upper)]
  plot_dt[, bar_max   := pmax(icer_lower, icer_upper)]
  plot_dt[, icer_range := bar_max - bar_min]

  # Separate parameters that bracket the baseline from those that do not
  valid_dt   <- plot_dt[bar_min <= icer_baseline & icer_baseline <= bar_max &
                          icer_range > 0]
  excluded_dt <- plot_dt[!(bar_min <= icer_baseline & icer_baseline <= bar_max) |
                           icer_range == 0]

  if (nrow(excluded_dt) > 0) {
    cat(sprintf("\n%s -- excluded parameters:\n", title_str))
    cat("  (baseline ICER outside [lower-bound ICER, upper-bound ICER];\n")
    cat("  incremental QALY changes sign across the parameter range)\n")
    for (j in seq_len(nrow(excluded_dt))) {
      cat(sprintf("  %-45s range [%10.0f, %10.0f]  baseline %10.0f\n",
                  excluded_dt$label[j],
                  excluded_dt$bar_min[j], excluded_dt$bar_max[j],
                  icer_baseline))
    }
    cat("\n")
  }

  if (nrow(valid_dt) == 0) {
    cat("No valid parameters for tornado diagram:", title_str, "\n")
    return(NULL)
  }

  # Sort ascending so the widest bar plots at the top in ggplot
  valid_dt <- valid_dt[order(icer_range)]
  valid_dt[, label := factor(label, levels = label)]

  ggplot(valid_dt, aes(y = label)) +
    geom_segment(
      aes(x = bar_min, xend = bar_max, yend = label),
      linewidth = 6, color = "#2B7BB9", alpha = 0.75
    ) +
    geom_vline(
      xintercept = icer_baseline,
      color = "red", linetype = "dashed", linewidth = 0.9
    ) +
    scale_x_continuous(labels = scales::comma) +
    labs(
      title    = title_str,
      subtitle = sprintf(
        "Baseline ICER = $%s / QALY",
        formatC(icer_baseline, format = "f", digits = 0, big.mark = ",")
      ),
      x       = "ICER ($ / QALY)",
      y       = NULL,
      caption = paste0(
        "Red dashed line = baseline ICER. ",
        "Parameters sorted by influence (widest bar at top). ",
        "Parameters excluded where delta-QALY changes sign."
      )
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title    = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      axis.text.y   = element_text(size = 9),
      plot.caption  = element_text(size = 8, hjust = 0)
    )
}

plot_tornado_arexvy <- make_tornado_plot(
  owsa_results,
  "icer_arexvy_at_lower",
  "icer_arexvy_at_upper",
  icer_arexvy_vs_soc,
  "Tornado Diagram: Arexvy vs Standard of Care"
)

plot_tornado_abrysvo <- make_tornado_plot(
  owsa_results,
  "icer_abrysvo_at_lower",
  "icer_abrysvo_at_upper",
  icer_abrysvo_vs_soc,
  "Tornado Diagram: Abrysvo vs Standard of Care"
)

pdf("rsv_owsa_tornado.pdf", width = 11, height = 7)
if (!is.null(plot_tornado_arexvy))  print(plot_tornado_arexvy)
if (!is.null(plot_tornado_abrysvo)) print(plot_tornado_abrysvo)
dev.off()
cat("Tornado diagrams saved to: rsv_owsa_tornado.pdf\n")

if (!is.null(plot_tornado_arexvy)) {
  ggsave("rsv_owsa_tornado_arexvy.png", plot = plot_tornado_arexvy,
         width = 11, height = 7, dpi = 300)
}
if (!is.null(plot_tornado_abrysvo)) {
  ggsave("rsv_owsa_tornado_abrysvo.png", plot = plot_tornado_abrysvo,
         width = 11, height = 7, dpi = 300)
}
cat("Tornado diagrams saved to: rsv_owsa_tornado_arexvy.png,",
    "rsv_owsa_tornado_abrysvo.png\n")

write_xlsx(
  list(
    "Arexvy OWSA"  = as.data.frame(
      owsa_arexvy_display[order(-`ICER Range`)]
    ),
    "Abrysvo OWSA" = as.data.frame(
      owsa_abrysvo_display[order(-`ICER Range`)]
    )
  ),
  "rsv_owsa_results.xlsx"
)
cat("OWSA results exported to: rsv_owsa_results.xlsx\n\n")

# ==============================================================================
# Step 8: Probabilistic Sensitivity Analysis (PSA)
# ==============================================================================

# ------------------------------------------------------------------------------
# Step 8a: Derive sampling parameters for "normal" distribution variables
# ------------------------------------------------------------------------------
#
# Sampling method selection for parameters labelled "normal" in the CSV:
#
#   rlnorm — used for parameters that are strictly positive with asymmetric
#     95% CIs, indicating the log-scale is the natural scale for uncertainty
#     (small rates, probabilities, costs, right-skewed durations). mu_log is
#     the midpoint of log(lower) and log(upper); sigma_log is the half log-
#     range divided by 1.96.
#
#   rnorm — used where the 95% CI is symmetric around the mean, the parameter
#     is bounded well away from zero, and the source literature reports a
#     symmetric interval on the natural scale:
#       kappa_arexvy / kappa_abrysvo: perfectly symmetric CIs (±0.009 /
#         ±0.011) from vaccine trial efficacy on the natural scale.
#       o_infected: perfectly symmetric CI (±3 days).
#       p_adverse_events: near-symmetric CI (-0.04 / +0.05), above zero.
#     Standard deviation is the half CI width divided by 1.96.
#
#   Fixed — tau_2 has no published CI and is held at its point estimate.
#     Note: o_infected_h bounds [2, 6.75] are IQR (25th-75th percentile),
#     not a 95% CI. The sigma_log denominator is 2 * qnorm(0.75) = 1.3490.
# ------------------------------------------------------------------------------

# Extract all parameters labelled "normal" in the distribution column
psa_normal_raw <- input_parameters[
  !is.na(input_parameters$distribution) &
    trimws(input_parameters$distribution) == "normal",
]

psa_normal_raw$value       <- as.numeric(psa_normal_raw$value)
psa_normal_raw$lower_bound <- suppressWarnings(
  as.numeric(psa_normal_raw$lower_bound)
)
psa_normal_raw$upper_bound <- suppressWarnings(
  as.numeric(psa_normal_raw$upper_bound)
)

# Split: parameters with bounds vs. those held fixed (tau_2: no bounds)
psa_normal_fixed  <- psa_normal_raw[
  is.na(psa_normal_raw$lower_bound) | is.na(psa_normal_raw$upper_bound),
]
psa_normal_bounds <- psa_normal_raw[
  !is.na(psa_normal_raw$lower_bound) & !is.na(psa_normal_raw$upper_bound),
]

# Split bounded parameters by sampling distribution
psa_rnorm_names <- c(
  "kappa_arexvy", "kappa_abrysvo", "o_infected", "p_adverse_events"
)
psa_normal_rnorm <- psa_normal_bounds[
  psa_normal_bounds$parameter %in% psa_rnorm_names,
]
psa_normal_rlnorm <- psa_normal_bounds[
  !psa_normal_bounds$parameter %in% psa_rnorm_names,
]

# rnorm: SD derived from 95% CI half-width
psa_normal_rnorm$sd <- with(psa_normal_rnorm,
  (upper_bound - lower_bound) / (2 * 1.96)
)

# rlnorm: log-scale mean is midpoint of log(lower) and log(upper).
# sigma_log denominator is 2*qnorm(0.75) for o_infected_h (IQR bounds)
# and 2*1.96 for all others (95% CI bounds).
psa_normal_rlnorm$mu_log <- with(psa_normal_rlnorm,
  (log(lower_bound) + log(upper_bound)) / 2
)
psa_normal_rlnorm$sigma_log <- ifelse(
  psa_normal_rlnorm$parameter == "o_infected_h",
  (log(psa_normal_rlnorm$upper_bound) -
     log(psa_normal_rlnorm$lower_bound)) / (2 * qnorm(0.75)),
  (log(psa_normal_rlnorm$upper_bound) -
     log(psa_normal_rlnorm$lower_bound)) / (2 * 1.96)
)

cat("=== PSA Step 8a: Normal distribution parameter derivations ===\n\n")

cat("--- rlnorm (lognormal; asymmetric CI or strictly positive) ---\n")
print(psa_normal_rlnorm[
  , c("parameter", "value", "lower_bound", "upper_bound",
      "mu_log", "sigma_log")
])

cat("\n--- rnorm (normal; symmetric CI, well-bounded from zero) ---\n")
print(psa_normal_rnorm[
  , c("parameter", "value", "lower_bound", "upper_bound", "sd")
])

cat("\n--- Fixed in PSA (no CI bounds provided) ---\n")
print(psa_normal_fixed[, c("parameter", "value")])

# ------------------------------------------------------------------------------
# Step 8b: Sample one draw for "normal" distribution variables
# ------------------------------------------------------------------------------
# For rlnorm parameters: draw on the log scale with rnorm(), then exponentiate
# to recover a natural-scale value. This two-step process is mathematically
# equivalent to rlnorm() and will be consistent with how gamma and beta draws
# are assembled in Steps 8c and 8d before the full loop in Step 8e.
# For rnorm parameters: draw directly on the natural scale.
# Fixed parameters return their point estimate unchanged.
# ------------------------------------------------------------------------------

draw_normal_params <- function() {

  # rlnorm: draw on log scale then exponentiate
  rlnorm_draws <- setNames(
    mapply(
      function(mu, sigma) exp(rnorm(1, mu, sigma)),
      psa_normal_rlnorm$mu_log,
      psa_normal_rlnorm$sigma_log
    ),
    psa_normal_rlnorm$parameter
  )

  # rnorm: draw directly on natural scale
  rnorm_draws <- setNames(
    mapply(
      function(mu, sd) rnorm(1, mu, sd),
      psa_normal_rnorm$value,
      psa_normal_rnorm$sd
    ),
    psa_normal_rnorm$parameter
  )

  # Fixed: return point estimate unchanged
  fixed_draws <- setNames(
    psa_normal_fixed$value,
    psa_normal_fixed$parameter
  )

  c(rlnorm_draws, rnorm_draws, fixed_draws)
}

# ------------------------------------------------------------------------------
# Step 8c: Derive sampling parameters for "gamma" distribution variables
# ------------------------------------------------------------------------------
# Gamma distributions are used for strictly positive costs and length-of-stay
# parameters.
#
# Three cases are handled:
#
#   Direct shape/scale — o_los_h: the lower_bound and upper_bound columns
#     contain the gamma shape and scale parameters directly (not CI bounds).
#     Verified: shape * scale = 1.2258 * 5.0582 = 6.2, matching the mean.
#
#   Method of moments from 95% CI — c_hfr, c_death: SE is the half CI
#     width divided by 1.96. Alpha equals (mean/SE) squared; rate lambda
#     equals mean divided by SE squared.
#
#   Method of moments, CV = 1 — all other gamma parameters have no bounds.
#     SE is assumed equal to the mean (coefficient of variation = 1), which
#     implies alpha = 1, yielding an exponential distribution. This is the
#     standard diffuse assumption when no variance data are available.
# ------------------------------------------------------------------------------

psa_gamma_raw <- input_parameters[
  !is.na(input_parameters$distribution) &
    trimws(input_parameters$distribution) == "gamma",
]

psa_gamma_raw$value       <- as.numeric(psa_gamma_raw$value)
psa_gamma_raw$lower_bound <- suppressWarnings(
  as.numeric(psa_gamma_raw$lower_bound)
)
psa_gamma_raw$upper_bound <- suppressWarnings(
  as.numeric(psa_gamma_raw$upper_bound)
)

# Separate o_los_h (direct shape/scale) from method-of-moments parameters
psa_gamma_direct <- psa_gamma_raw[
  psa_gamma_raw$parameter == "o_los_h",
]
psa_gamma_mom <- psa_gamma_raw[
  psa_gamma_raw$parameter != "o_los_h",
]

# SE from 95% CI where bounds exist; SE = mean (CV = 1) otherwise
psa_gamma_mom$se <- ifelse(
  !is.na(psa_gamma_mom$lower_bound) & !is.na(psa_gamma_mom$upper_bound),
  (psa_gamma_mom$upper_bound - psa_gamma_mom$lower_bound) / (2 * 1.96),
  psa_gamma_mom$value
)
psa_gamma_mom$alpha  <- with(psa_gamma_mom, (value / se)^2)
psa_gamma_mom$lambda <- with(psa_gamma_mom, value / se^2)

cat("=== PSA Step 8c: Gamma distribution parameter derivations ===\n\n")

cat("--- Direct shape/scale (o_los_h only) ---\n")
print(psa_gamma_direct[
  , c("parameter", "value", "lower_bound", "upper_bound")
])

cat("\n--- Method of moments with CI bounds (SE from bounds) ---\n")
psa_gamma_mom_bounds <- psa_gamma_mom[
  !is.na(psa_gamma_mom$lower_bound) & !is.na(psa_gamma_mom$upper_bound),
]
print(psa_gamma_mom_bounds[
  , c("parameter", "value", "lower_bound", "upper_bound",
      "se", "alpha", "lambda")
])

cat("\n--- Method of moments, no bounds (SE = mean; alpha = 1,")
cat(" exponential) ---\n")
psa_gamma_mom_no_bounds <- psa_gamma_mom[
  is.na(psa_gamma_mom$lower_bound) | is.na(psa_gamma_mom$upper_bound),
]
print(psa_gamma_mom_no_bounds[
  , c("parameter", "value", "se", "alpha", "lambda")
])

# ------------------------------------------------------------------------------
# Step 8d: Sample one draw for "gamma" distribution variables
# ------------------------------------------------------------------------------
# o_los_h: draw with shape and scale taken directly from the parameter table.
# All others: draw with method-of-moments shape and rate.
# ------------------------------------------------------------------------------

draw_gamma_params <- function() {

  # o_los_h: direct shape and scale
  direct_draws <- setNames(
    mapply(
      function(shape, scale) rgamma(1, shape = shape, scale = scale),
      psa_gamma_direct$lower_bound,
      psa_gamma_direct$upper_bound
    ),
    psa_gamma_direct$parameter
  )

  # All others: method-of-moments shape and rate
  mom_draws <- setNames(
    mapply(
      function(alpha, lambda) rgamma(1, shape = alpha, rate = lambda),
      psa_gamma_mom$alpha,
      psa_gamma_mom$lambda
    ),
    psa_gamma_mom$parameter
  )

  c(direct_draws, mom_draws)
}

# ------------------------------------------------------------------------------
# Step 8e: Derive sampling parameters for "beta" distribution variables
# ------------------------------------------------------------------------------
# Beta distributions are used for utilities and disutilities.
#
# Two cases are handled:
#
#   Natural [0, 1] — positive utility parameters: beta is fit directly on
#     [0, 1] via method of moments. SE is the half CI width divided by 1.96.
#     Concentration phi equals mu times (1 - mu) divided by SE squared,
#     minus 1. Alpha equals mu times phi; beta_shape equals (1-mu) times phi.
#
#   Scaled [L, U] — negative disutility parameters (u_adverse_events,
#     u_disutility_infection, u_disutility_hospitalization): the parameter
#     space [L, U] is transformed to [0, 1] for fitting. mu_s is the
#     position of the mean within [L, U]. Because SE equals (U-L)/(2*1.96),
#     the scaled SE simplifies to 1/(2*1.96) for all scaled parameters.
#     Method-of-moments is applied on the scaled space. Draws are then
#     back-transformed by multiplying by (U - L) and adding L.
# ------------------------------------------------------------------------------

psa_beta_raw <- input_parameters[
  !is.na(input_parameters$distribution) &
    trimws(input_parameters$distribution) == "beta",
]

psa_beta_raw$value       <- as.numeric(psa_beta_raw$value)
psa_beta_raw$lower_bound <- suppressWarnings(
  as.numeric(psa_beta_raw$lower_bound)
)
psa_beta_raw$upper_bound <- suppressWarnings(
  as.numeric(psa_beta_raw$upper_bound)
)

# Natural [0, 1]: positive utility parameters
psa_beta_natural <- psa_beta_raw[psa_beta_raw$value >= 0, ]
# Scaled [L, U]: negative disutility parameters
psa_beta_scaled  <- psa_beta_raw[psa_beta_raw$value <  0, ]

# Natural: method of moments from mean and SE
psa_beta_natural$se <- with(psa_beta_natural,
  (upper_bound - lower_bound) / (2 * 1.96)
)
psa_beta_natural$phi <- with(psa_beta_natural,
  value * (1 - value) / se^2 - 1
)
psa_beta_natural$alpha      <- with(psa_beta_natural, value * phi)
psa_beta_natural$beta_shape <- with(psa_beta_natural, (1 - value) * phi)

# Scaled: transform mu to [0, 1]; SE_s = 1 / (2 * 1.96) for all
psa_beta_scaled$mu_s  <- with(psa_beta_scaled,
  (value - lower_bound) / (upper_bound - lower_bound)
)
psa_beta_scaled$se_s  <- 1 / (2 * 1.96)
psa_beta_scaled$phi_s <- with(psa_beta_scaled,
  mu_s * (1 - mu_s) / se_s^2 - 1
)
psa_beta_scaled$alpha      <- with(psa_beta_scaled, mu_s * phi_s)
psa_beta_scaled$beta_shape <- with(psa_beta_scaled, (1 - mu_s) * phi_s)

cat("=== PSA Step 8e: Beta distribution parameter derivations ===\n\n")

cat("--- Natural [0, 1] parameters ---\n")
print(psa_beta_natural[
  , c("parameter", "value", "lower_bound", "upper_bound",
      "alpha", "beta_shape")
])

cat("\n--- Scaled [L, U] parameters (negative disutilities) ---\n")
print(psa_beta_scaled[
  , c("parameter", "value", "lower_bound", "upper_bound",
      "mu_s", "alpha", "beta_shape")
])

# ------------------------------------------------------------------------------
# Step 8f: Sample one draw for "beta" distribution variables
# ------------------------------------------------------------------------------
# Natural parameters: draw directly from beta on [0, 1].
# Scaled parameters: draw from beta on [0, 1], back-transform to [L, U].
# ------------------------------------------------------------------------------

draw_beta_params <- function() {

  # Natural [0, 1]: draw directly
  natural_draws <- setNames(
    mapply(
      function(alpha, beta_shape) rbeta(1, alpha, beta_shape),
      psa_beta_natural$alpha,
      psa_beta_natural$beta_shape
    ),
    psa_beta_natural$parameter
  )

  # Scaled [L, U]: draw on [0, 1] then back-transform
  scaled_draws <- setNames(
    mapply(
      function(alpha, beta_shape, lb, ub) {
        rbeta(1, alpha, beta_shape) * (ub - lb) + lb
      },
      psa_beta_scaled$alpha,
      psa_beta_scaled$beta_shape,
      psa_beta_scaled$lower_bound,
      psa_beta_scaled$upper_bound
    ),
    psa_beta_scaled$parameter
  )

  c(natural_draws, scaled_draws)
}

# ==============================================================================
# Step 8g: PSA Loop (parallelised — 4 workers via mclapply)
# ==============================================================================
# run_one_psa() encapsulates one complete simulation iteration and is mapped
# across iterations by mclapply(). Fork-based parallelism (macOS / Linux)
# copies the full parent R environment to each worker, so all parameter tables
# and model functions are available without explicit export.
#
# RNG: RNGkind("L'Ecuyer-CMRG") assigns each forked worker its own
#   independent substream of the same seed, giving reproducible results
#   regardless of the number of cores used.
#
# Console logging: cat() output from child processes is written to the
#   parent's stdout in real time. Lines from different workers may interleave
#   but each line is atomic so readability is preserved.
# ==============================================================================

n_psa   <- 1000  # Number of PSA simulations
n_cores <- 4     # Parallel workers (adjust to available cores)

RNGkind("L'Ecuyer-CMRG")  # Reproducible parallel RNG
set.seed(42)

# run_one_psa() — one PSA iteration, returns a single data.table row
run_one_psa <- function(sim) {

  # Suppress run_simulation() console output within this worker
  run_silent <- function(...) {
    sink(tempfile())
    on.exit(sink(), add = TRUE)
    suppressWarnings(run_simulation(...))
  }

  # ------------------------------------------------------------------
  # Draw one complete parameter set
  # ------------------------------------------------------------------
  p_norm  <- draw_normal_params()
  p_gamma <- draw_gamma_params()
  p_beta  <- draw_beta_params()

  o_infected_s  <- p_norm[["o_infected"]]
  attack_rate_s <- p_norm[["attack_rate"]]
  tau_1_s      <- 7 / o_infected_s
  beta_trans_s <- tau_1_s * (-log(1 - attack_rate_s) / attack_rate_s)

  # ------------------------------------------------------------------
  # Build shared arguments (identical across all three scenarios)
  # ------------------------------------------------------------------
  shared_args <- list(
    p_60 = p_60, v_0_p = v_0_p,
    i_0_p = p_norm[["i_0_p"]], h_0_p = p_norm[["h_0_p"]],
    beta = beta_trans_s, sigma = sigma,
    phi = p_norm[["phi"]],
    tau_1 = tau_1_s,
    tau_2 = p_norm[["tau_2"]], tau_3 = p_norm[["tau_3"]],
    rho = p_norm[["rho"]], mu = mu,
    delta_m = delta_m,
    delta_s = p_norm[["delta_s"]], delta_sq = delta_sq,
    gamma_r = gamma_r, psi = psi, gamma_v = gamma_v,
    u_s_vn = p_beta[["u_s_vn"]], u_v = p_beta[["u_v"]],
    u_s_ve = p_beta[["u_s_ve"]], u_e = p_beta[["u_e"]],
    u_i = p_beta[["u_i"]], u_h = p_beta[["u_h"]],
    u_sq = p_beta[["u_sq"]], u_r = p_beta[["u_r"]], u_d = u_d,
    u_adverse_events = p_beta[["u_adverse_events"]],
    p_adverse_events = p_norm[["p_adverse_events"]],
    u_disutility_infection = p_beta[["u_disutility_infection"]],
    u_disutility_hospitalization =
      p_beta[["u_disutility_hospitalization"]],
    duration_adverse_events = duration_adverse_events,
    c_vac_admin = p_norm[["c_vac_admin"]],
    c_hosp = p_gamma[["c_hosp"]],
    c_pl = p_gamma[["c_pl"]], c_care = p_gamma[["c_care"]],
    c_hfr = p_gamma[["c_hfr"]], c_death = p_gamma[["c_death"]],
    o_infected = o_infected_s, o_los_h = p_gamma[["o_los_h"]],
    p_discount_wk = p_discount_wk,
    max_periods = max_periods,
    susceptible_threshold = susceptible_threshold
  )

  # ------------------------------------------------------------------
  # Run all three scenarios
  # ------------------------------------------------------------------
  res_soc <- do.call(run_silent, c(shared_args, list(
    c_vaccine = 0, kappa = 0,
    use_vaccine = FALSE, vaccine_name = "SOC_PSA"
  )))
  res_arexvy <- do.call(run_silent, c(shared_args, list(
    c_vaccine = c_arexvy,
    kappa = p_norm[["kappa_arexvy"]],
    use_vaccine = TRUE, vaccine_name = "Arexvy_PSA"
  )))
  res_abrysvo <- do.call(run_silent, c(shared_args, list(
    c_vaccine = c_abrysvo,
    kappa = p_norm[["kappa_abrysvo"]],
    use_vaccine = TRUE, vaccine_name = "Abrysvo_PSA"
  )))

  # ------------------------------------------------------------------
  # Aggregate and compute ICERs
  # ------------------------------------------------------------------
  cost_soc     <- sum(res_soc$cost_total)
  qaly_soc     <- sum(res_soc$utility_total)
  cost_arexvy  <- sum(res_arexvy$cost_total)
  qaly_arexvy  <- sum(res_arexvy$utility_total)
  cost_abrysvo <- sum(res_abrysvo$cost_total)
  qaly_abrysvo <- sum(res_abrysvo$utility_total)

  icer_arexvy_s  <- (cost_arexvy  - cost_soc) / (qaly_arexvy  - qaly_soc)
  icer_abrysvo_s <- (cost_abrysvo - cost_soc) / (qaly_abrysvo - qaly_soc)

  cat(sprintf(
    "  Sim %4d / %d | ICER Arexvy = $%s | ICER Abrysvo = $%s\n",
    sim, n_psa,
    format(round(icer_arexvy_s),  big.mark = ",", scientific = FALSE),
    format(round(icer_abrysvo_s), big.mark = ",", scientific = FALSE)
  ))

  # ------------------------------------------------------------------
  # Return one results row
  # ------------------------------------------------------------------
  data.table(
    sim                          = sim,
    # Normal distribution parameters
    i_0_p                        = p_norm[["i_0_p"]],
    h_0_p                        = p_norm[["h_0_p"]],
    phi                          = p_norm[["phi"]],
    rho                          = p_norm[["rho"]],
    delta_s                      = p_norm[["delta_s"]],
    attack_rate                  = attack_rate_s,
    c_vac_admin                  = p_norm[["c_vac_admin"]],
    o_infected                   = o_infected_s,
    o_Infected_ma                = p_norm[["o_Infected_ma"]],
    o_infected_h                 = p_norm[["o_infected_h"]],
    kappa_arexvy                 = p_norm[["kappa_arexvy"]],
    kappa_abrysvo                = p_norm[["kappa_abrysvo"]],
    p_adverse_events             = p_norm[["p_adverse_events"]],
    tau_2                        = p_norm[["tau_2"]],
    tau_3                        = p_norm[["tau_3"]],
    # Derived parameters
    tau_1                        = tau_1_s,
    beta_transmission            = beta_trans_s,
    # Gamma distribution parameters
    o_los_h                      = p_gamma[["o_los_h"]],
    o_los_icu                    = p_gamma[["o_los_icu"]],
    c_office                     = p_gamma[["c_office"]],
    c_ed                         = p_gamma[["c_ed"]],
    c_hosp                       = p_gamma[["c_hosp"]],
    c_icu_v                      = p_gamma[["c_icu_v"]],
    c_icu                        = p_gamma[["c_icu"]],
    c_pl                         = p_gamma[["c_pl"]],
    c_care                       = p_gamma[["c_care"]],
    c_hfr                        = p_gamma[["c_hfr"]],
    c_death                      = p_gamma[["c_death"]],
    # Beta distribution parameters
    u_s_vn                       = p_beta[["u_s_vn"]],
    u_v                          = p_beta[["u_v"]],
    u_s_ve                       = p_beta[["u_s_ve"]],
    u_e                          = p_beta[["u_e"]],
    u_i                          = p_beta[["u_i"]],
    u_h                          = p_beta[["u_h"]],
    u_sq                         = p_beta[["u_sq"]],
    u_r                          = p_beta[["u_r"]],
    u_adverse_events             = p_beta[["u_adverse_events"]],
    u_disutility_infection       = p_beta[["u_disutility_infection"]],
    u_disutility_hospitalization =
      p_beta[["u_disutility_hospitalization"]],
    # Costs and QALYs
    total_cost_soc               = cost_soc,
    total_qaly_soc               = qaly_soc,
    total_cost_arexvy            = cost_arexvy,
    total_qaly_arexvy            = qaly_arexvy,
    total_cost_abrysvo           = cost_abrysvo,
    total_qaly_abrysvo           = qaly_abrysvo,
    # ICERs vs standard of care
    icer_arexvy_vs_soc           = icer_arexvy_s,
    icer_abrysvo_vs_soc          = icer_abrysvo_s
  )
}

# run_psa() dispatches run_one_psa() across n_cores workers and combines rows
run_psa <- function(n_sim, n_cores = 4) {
  cat(sprintf(
    "Running PSA: %d simulations on %d cores\n", n_sim, n_cores
  ))
  rbindlist(mclapply(
    seq_len(n_sim),
    run_one_psa,
    mc.cores     = n_cores,
    mc.set.seed  = TRUE
  ))
}

psa_results_dt <- run_psa(n_psa, n_cores)
cat("PSA complete.\n\n")

write_xlsx(
  list("PSA Results" = as.data.frame(psa_results_dt)),
  "rsv_psa_results.xlsx"
)
cat("PSA results exported to: rsv_psa_results.xlsx\n\n")

# ==============================================================================
# Step 8i: Cost-Effectiveness Plane (PSA)
# ==============================================================================
# Best-practice CE plane elements:
#   Scatter points   — one per simulation, semi-transparent to reveal density
#   95% ellipses     — joint confidence region for incremental cost and effect
#   Mean estimate    — open circle marking the mean incremental pair
#   WTP threshold    — long-dashed line through the origin, slope = $100k/QALY;
#                      points below this line are cost-effective at that WTP
#   Reference axes   — dashed lines at x = 0 and y = 0 to delineate quadrants
#
# NMB-based cost-effectiveness: a simulation is cost-effective when
#   WTP * inc_qaly - inc_cost > 0, which is consistent across all quadrants.
#
# Both vaccines are plotted together for direct visual comparison.
# Colourblind-safe palette (Okabe-Ito blue and vermillion).
# ==============================================================================

# Compute incremental costs and QALYs vs SOC for each simulation
psa_ce_data <- rbind(
  data.table(
    sim      = psa_results_dt$sim,
    vaccine  = "Arexvy",
    inc_qaly = psa_results_dt$total_qaly_arexvy -
      psa_results_dt$total_qaly_soc,
    inc_cost = psa_results_dt$total_cost_arexvy -
      psa_results_dt$total_cost_soc
  ),
  data.table(
    sim      = psa_results_dt$sim,
    vaccine  = "Abrysvo",
    inc_qaly = psa_results_dt$total_qaly_abrysvo -
      psa_results_dt$total_qaly_soc,
    inc_cost = psa_results_dt$total_cost_abrysvo -
      psa_results_dt$total_cost_soc
  )
)

# Mean incremental estimates per vaccine (plotted as open circles)
psa_ce_means <- psa_ce_data[
  , .(inc_qaly = mean(inc_qaly), inc_cost = mean(inc_cost)),
  by = vaccine
]

# Proportion of simulations cost-effective at the WTP threshold via NMB
psa_ce_data[
  , cost_effective := (wtp_threshold * inc_qaly - inc_cost) > 0
]
psa_pct_ce <- psa_ce_data[
  , .(pct_cost_effective = round(mean(cost_effective) * 100, 1)),
  by = vaccine
]

cat("=== PSA Step 8i: Cost-Effectiveness Summary ===\n\n")
cat("Mean incremental estimates vs SOC:\n")
print(psa_ce_means)
cat(sprintf(
  "\nProportion cost-effective at WTP = $%s/QALY:\n",
  format(wtp_threshold, big.mark = ",")
))
print(psa_pct_ce)

# Colourblind-safe palette (Okabe-Ito)
vaccine_colours <- c("Arexvy" = "#0072B2", "Abrysvo" = "#D55E00")

# Build CE plane
ce_plane <- ggplot(
  psa_ce_data,
  aes(x = inc_qaly, y = inc_cost, colour = vaccine)
) +
  # Quadrant reference lines
  geom_hline(
    yintercept = 0, colour = "grey60",
    linewidth = 0.4, linetype = "dashed"
  ) +
  geom_vline(
    xintercept = 0, colour = "grey60",
    linewidth = 0.4, linetype = "dashed"
  ) +
  # WTP threshold through origin (slope = WTP per QALY)
  geom_abline(
    slope = wtp_threshold, intercept = 0,
    colour = "grey20", linetype = "longdash", linewidth = 0.7
  ) +
  # Simulation scatter (semi-transparent to show density)
  geom_point(alpha = 0.15, size = 0.8, shape = 16) +
  # 95% joint confidence ellipse
  stat_ellipse(level = 0.95, linewidth = 1.0) +
  # Mean incremental estimate (open circle, white fill)
  geom_point(
    data = psa_ce_means,
    aes(x = inc_qaly, y = inc_cost),
    size = 4, shape = 21, fill = "white", stroke = 1.5
  ) +
  # Scales
  scale_colour_manual(values = vaccine_colours) +
  scale_y_continuous(labels = scales::label_dollar(big.mark = ",")) +
  scale_x_continuous(labels = scales::label_comma()) +
  # Labels
  labs(
    x       = "Incremental QALYs vs Standard of Care",
    y       = "Incremental Cost (USD) vs Standard of Care",
    colour  = NULL,
    caption = paste0(
      "Long-dashed line: WTP threshold ($",
      format(wtp_threshold, big.mark = ","),
      "/QALY). ",
      "Open circles: mean estimates. ",
      "Ellipses: 95% confidence regions."
    )
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position       = "bottom",
    legend.key.width      = unit(1.5, "lines"),
    panel.grid.minor      = element_blank(),
    plot.caption          = element_text(size = 8, colour = "grey40"),
    plot.caption.position = "plot"
  )

# Export
ggsave(
  "rsv_psa_ce_plane.pdf", plot = ce_plane,
  width = 8, height = 7, device = "pdf"
)
ggsave(
  "rsv_psa_ce_plane.png", plot = ce_plane,
  width = 8, height = 7, dpi = 300
)
cat("CE plane saved to: rsv_psa_ce_plane.pdf, rsv_psa_ce_plane.png\n\n")

# ==============================================================================
# Step 9: Post-Simulation Validation Report
# ==============================================================================

# Run comprehensive validation report
validation_results <- generate_validation_report()