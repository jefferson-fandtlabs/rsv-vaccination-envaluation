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

# ==============================================================================
# Step 2: Define variables
# ==============================================================================

# ==============================================================================
# Step 2a: Initial Population Distribution
# ==============================================================================

# Total U.S. Population Age 60+
p_60 <- 82973748

# Vaccinated
v_0_p <- .162 # Percentage of p_60 that represents the vaccinated population for vaccine scenario at time 0

# Infectious
i_0_p <- .0162 # Percentage of p_60 that represents the infectious population at time 0

# Hospitalized
h_0_p <- 0.00243 # Percentage of p_60 that represents the hospitalized population at time 0

# Suceptible
s_vn_0_p <- (1 - i_0_p - h_0_p) # Percentage of p_60 that represents the initial susceptible, vaccine-naive, no vaccine scenario at time 0
s_vn_0_v_p <- (1 - i_0_p - h_0_p - v_0_p) # Percentage of p_60 that represents the initial susceptible, vaccine-naive, vaccine scenario at time 0

# ==============================================================================
# Step 2b: Transition Probabilities
# ==============================================================================

# Vaccine Breakthrough Rate (Kappa)
kappa_arexvy <- .211 # Vaccine Breakthrough Rate for Arexvy. Range from 20.2% to 22.0%
kappa_abrysvo <- .222 # Vaccine Breakthrough Rate for Abrysvo. Range from 21.1% to 23.3%

# RSV Vaccination Rate in US
psi <- v_0_p / 52 # Weekly vaccination rate derived from annual rate of 16.2%

# Waning immunity
gamma_v <- 0.015 # Waning vaccine confered immunity
gamma_r <- 0.015 # Waning natual immunity

# Conversion rate
sigma <- 1 # Rate of conversion from exposed to infected

# Hospitalization rate
phi <- 0.0015 # Rate of hospitalization

# Recovery rate
tau_2 <- 0.904 # Recovery Rate from hospitalization
tau_3 <- 1     # Recovery from sequelae

# Chronic State
rho <- 0.18 # Rate of transition to chronic state

# Mortality
mu <- 0.00139 # Age adjusted weekly mortality

# RSV Mortality
delta_m <- 0 # Mortality rate from mild RSV
delta_s <- 0.0713 # Mortality rate from severe RSV
delta_sq <- 0 # Mortality rate from RSV sequelae

# ==============================================================================
# Step 2c: Health Utility For Each State
# ==============================================================================

# Each is represented as a unit of utility that must then be multiplied by the time in that state.
u_s_vn <- 0.815 # Susceptible, vaccine naive
u_v <- 0.815 # Vaccinated (healthy, protected) - same as susceptible
u_s_ve <- 0.815 # Susceptible, vaccine experienced
u_e <- 0.815 # Exposed
u_i <- 0.815 # Infected (baseline utility for ongoing state)
u_h <- 0.815 # Hospitalized (baseline utility for ongoing state)
u_sq <- 0.080 # sequalae
u_r <- 0.815 # recovered
u_d <- 0 # Death, absorbing

# Transient disutilities (applied only when entering acute states)
u_adverse_events <- -0.08 # Annual utility weight decrement during adverse events (EQ-5D scale)
p_adverse_events <- 0.20  # Proportion of vaccinated people who experience an adverse event
u_disutility_infection <- -0.010 # Weekly QALY disutility during acute infection phase
u_disutility_hospitalization <- -0.017 # Weekly QALY disutility during acute hospitalization phase
duration_adverse_events <- 1 # Duration of adverse events in weeks

# ==============================================================================
# Step 2d: Costs
# ==============================================================================

# Vaccine Costs
c_abrysvo <- 319 # Cost of Abrysvo vaccine dose
c_arexvy <- 280 # Cost of Abrexvy dose
c_vac_admin <- 50 # Cost to administer vaccines

# RSV Costs
c_office <- 126 # Cost per office visit
c_ed <- 649 # Cost per ed day
c_hosp <- 1116 # Cost per hospital day
c_icu_v <- 4870 # Cost per icu day with mechanincal vent
c_icu <- 3348 # Cost per icu day without mechanical vent

# Indirect Costs
c_pl <- 275 # Cost per case of productivity loss
c_care <- 583 # Cost per case of caregiver time
c_hfr <- 92 # Cost of household financial risk (middle of distribution)

# Death Costs
c_death <- 3799 # Cost to medicare for patient death (middle of distribution)

# ==============================================================================
# Step 2e: Outcomes
# ==============================================================================

o_infected    <- 5  # Days of infection without medical attention
o_Infected_ma <- 10 # Days of infection with medical attention
o_infected_h <- 4 # Days of infection to hospitalization
o_los_h <- 6.2 # Days of hospitalization without ICU
o_los_icu <- 4.5 # Days in ICU
# o_los_h_pre_icu <- # TODO: Days in hospital prior to ICU admission (value not yet defined)
# o_los_h_icu <-    # TODO: Days in hospital after being admitted to the ICU (value not yet defined)

# ==============================================================================
# Step 2f: Discount
# ==============================================================================

p_discount_yr <- 0.03
p_discount_wk <- p_discount_yr / 52  # Weekly discount rate derived from annual rate

# ==============================================================================
# Step 2g: Simulation Controls
# ==============================================================================

max_periods <- 20000            # Safety limit to prevent infinite loops
susceptible_threshold <- 1      # Threshold below which susceptible population is considered depleted

# ==============================================================================
# Step 2h: Derived Transition Probabilities
# ==============================================================================

# Recovery rate from mild illness: expressed as a weekly rate using the
# infectious period in days (o_infected defined in Step 2e).
tau_1 <- 7 / o_infected

# Infection rate: derived from the RSV-ARI annual attack rate using the SIR
# final size equation. R0 = -ln(1 - A) / A, where A = i_0_p.
# beta = tau_1 * R0
attack_rate <- 0.0162
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
# Step 7: Post-Simulation Validation Report
# ==============================================================================

# Run comprehensive validation report
validation_results <- generate_validation_report()