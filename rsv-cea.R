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

# Suceptible
s_vn_0_p <- .98137 # Percentage of p_60 that represents the initial susceptible, vaccine-naive, no vaccine scenario at time 0
s_vn_0_v_p <- .81737 # Percentage of p_60 that represents the initial susceptible, vaccine-naive, vaccine scenario at time 0

# Vaccinated
v_0_p <- .164 # Percentage of p_60 that represents the vaccinated population for vaccine scenario at time 0

# Infectious
i_0_p <- .0162 # Percentage of p_60 that represents the infectious population at time 0

# Hospitalized
h_0_p <- 0.00243 # Percentage of p_60 that represents the hospitalized population at time 0

# ==============================================================================
# Step 2b: Transition Probabilities
# ==============================================================================

# Vaccine Breakthrough Rate (Kappa)
kappa_arexvy <- .0039 # Vaccine Breakthrough Rate for Arexvy.
kappa_abrysvo <- .00395 # Vaccine Breakthrough Rate for Abrysvo

# RSV Vaccination Rate in US
psi <- 0.42 # Vaccination rate for RSV among adults in the use 60+

# Waning immunity
gamma_v <- 0.015 # Waning vaccine confered immunity
gamma_r <- 0.015 # Waning natual immunity

# Infection rate
beta <- 0.00031 # Rate of infection of RSV

# Conversion rate
sigma <- 1 # Rate of conversion from exposed to infected

# Hospitalization rate
phi <- 0.0015 # Rate of hospitalization

# Recovery rate
tau_1 <- 1 # Recovery rate from mild illness
tau_2 <- 0.904 # Recovery Rate from hospitalization
tau_3 <- 1 # Recovery from sequelae

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
u_i <- -0.50 # infected
u_h <- -0.50 # hospitalized
u_sq <- 0.080 # sequalae
u_r <- 0.815 # recovered
u_d <- 0 # Death, absorbing

# Vaccination adverse events (transient disutility)
u_adverse_events <- -0.08 # Disutility during adverse events (side effects)
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

o_infected <- 5 # Days of infection without medical attention
o_Infected_ma <- 10 # Days of infection with medical attention
o_infected_h <- 4 # Days of infection to hospitalization
o_los_h <- 6.2 # Days of hospitalization without ICU
o_los_icu <- 4.5 # Days in ICU
o_los_h_pre_icu <- # Days in hospital prior to ICU admission
o_los_h_icu <- # Days in hospital after being admitted to the ICU.

# ==============================================================================
# Step 3: Model - SOC
# ==============================================================================

# ==============================================================================
# Step 3a: Define Validation Functions
# ==============================================================================

# Function to perform quality checks during simulation
validate_period <- function(scenario_name, period, row_data, initial_pop, prev_deaths = 0) {

  issues <- character(0)

  # Extract populations (ensure scalar values)
  s_vn <- row_data$s_vn_end[1]
  s_ve <- row_data$s_ve_end[1]
  v <- row_data$v_end[1]
  e <- row_data$e_end[1]
  i <- row_data$i_end[1]
  h <- row_data$h_end[1]
  r <- row_data$r_end[1]
  sq <- row_data$sq_end[1]
  d <- row_data$d_end[1]

  # Check 1: Population Conservation
  total_pop <- s_vn + s_ve + v + e + i + h + r + sq + d
  pop_diff <- abs(total_pop - initial_pop)
  if (pop_diff > 1) {  # Allow for small rounding errors
    issues <- c(issues, sprintf("Population not conserved: %.2f (expected %.2f, diff %.2f)",
                                total_pop, initial_pop, pop_diff))
  }

  # Check 2: Non-Negativity
  if (s_vn < 0) issues <- c(issues, sprintf("S_VN negative: %.2f", s_vn))
  if (s_ve < 0) issues <- c(issues, sprintf("S_VE negative: %.2f", s_ve))
  if (v < 0) issues <- c(issues, sprintf("V negative: %.2f", v))
  if (e < 0) issues <- c(issues, sprintf("E negative: %.2f", e))
  if (i < 0) issues <- c(issues, sprintf("I negative: %.2f", i))
  if (h < 0) issues <- c(issues, sprintf("H negative: %.2f", h))
  if (r < 0) issues <- c(issues, sprintf("R negative: %.2f", r))
  if (sq < 0) issues <- c(issues, sprintf("SQ negative: %.2f", sq))
  if (d < 0) issues <- c(issues, sprintf("D negative: %.2f", d))

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
  total_soc <- final_soc$s_vn_end + final_soc$s_ve_end + final_soc$v_end +
               final_soc$e_end + final_soc$i_end + final_soc$h_end +
               final_soc$r_end + final_soc$sq_end + final_soc$d_end
  check("SOC: Population conserved",
        abs(total_soc - initial_pop) < 100,
        sprintf("Total: %.0f, Expected: %.0f, Diff: %.0f", total_soc, initial_pop, total_soc - initial_pop))

  final_arexvy <- results_arexvy_dt[nrow(results_arexvy_dt)]
  total_arexvy <- final_arexvy$s_vn_end + final_arexvy$s_ve_end + final_arexvy$v_end +
                  final_arexvy$e_end + final_arexvy$i_end + final_arexvy$h_end +
                  final_arexvy$r_end + final_arexvy$sq_end + final_arexvy$d_end
  check("Arexvy: Population conserved",
        abs(total_arexvy - initial_pop) < 100,
        sprintf("Total: %.0f, Expected: %.0f, Diff: %.0f", total_arexvy, initial_pop, total_arexvy - initial_pop))

  final_abrysvo <- results_abrysvo_dt[nrow(results_abrysvo_dt)]
  total_abrysvo <- final_abrysvo$s_vn_end + final_abrysvo$s_ve_end + final_abrysvo$v_end +
                   final_abrysvo$e_end + final_abrysvo$i_end + final_abrysvo$h_end +
                   final_abrysvo$r_end + final_abrysvo$sq_end + final_abrysvo$d_end
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
                summary_arexvy$total_cost/1e6, summary_soc$total_cost/1e6,
                arexvy_cost_diff/1e6, arexvy_cost_status))

  abrysvo_cost_diff <- summary_abrysvo$total_cost - summary_soc$total_cost
  abrysvo_cost_status <- ifelse(abrysvo_cost_diff < 0, " (Cost-Saving)", " (Higher Cost)")
  check("Abrysvo: Different total cost than SOC",
        summary_abrysvo$total_cost != summary_soc$total_cost,
        sprintf("Abrysvo: $%.2fM, SOC: $%.2fM, Diff: $%.2fM%s",
                summary_abrysvo$total_cost/1e6, summary_soc$total_cost/1e6,
                abrysvo_cost_diff/1e6, abrysvo_cost_status))

  # Vaccine scenarios should have different QALYs than SOC
  check("Arexvy: Different QALYs than SOC",
        summary_arexvy$total_qaly != summary_soc$total_qaly,
        sprintf("Diff: %.2fM QALYs", (summary_arexvy$total_qaly - summary_soc$total_qaly)/1e6))
  check("Abrysvo: Different QALYs than SOC",
        summary_abrysvo$total_qaly != summary_soc$total_qaly,
        sprintf("Diff: %.2fM QALYs", (summary_abrysvo$total_qaly - summary_soc$total_qaly)/1e6))

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
        sprintf("$%.2fM", summary_soc$total_cost/1e6))
  check("Arexvy: Total cost is reasonable",
        summary_arexvy$total_cost > 0 && summary_arexvy$total_cost < 1e12,
        sprintf("$%.2fM", summary_arexvy$total_cost/1e6))
  check("Abrysvo: Total cost is reasonable",
        summary_abrysvo$total_cost > 0 && summary_abrysvo$total_cost < 1e12,
        sprintf("$%.2fM", summary_abrysvo$total_cost/1e6))

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
# Step 3b: Initialize Results Storage
# ==============================================================================

# Initialize empty data.table to store period-by-period results
# This table will grow dynamically as simulation progresses until no susceptible population remains
results_soc_dt <- data.table(
  period = integer(),                    # Period number (sequential counter)
  s_vn_start = numeric(),                # Starting susceptible, vaccine-naive population
  s_ve_start = numeric(),                # Starting susceptible, vaccine-experienced population
  v_start = numeric(),                   # Starting vaccinated population
  e_start = numeric(),                   # Starting exposed population
  i_start = numeric(),                   # Starting infectious population
  h_start = numeric(),                   # Starting hospitalized population
  r_start = numeric(),                   # Starting recovered population
  sq_start = numeric(),                  # Starting sequelae population
  d_start = numeric(),                   # Starting deceased population
  s_vn_end = numeric(),                  # Ending susceptible, vaccine-naive population
  s_ve_end = numeric(),                  # Ending susceptible, vaccine-experienced population
  v_end = numeric(),                     # Ending vaccinated population
  e_end = numeric(),                     # Ending exposed population
  i_end = numeric(),                     # Ending infectious population
  h_end = numeric(),                     # Ending hospitalized population
  r_end = numeric(),                     # Ending recovered population
  sq_end = numeric(),                    # Ending sequelae population
  d_end = numeric(),                     # Ending deceased population
  cost_total = numeric(),                # Total cost incurred during period
  utility_total = numeric()              # Total utility (QALYs) accrued during period
)

# ==============================================================================
# Step 3b: Initialize Period 0 (Starting Condition)
# ==============================================================================

# Initial Compartment Conditions
s_vn_start <- p_60 * s_vn_0_p # Number of people in suceptable for standard of care
i_start <- p_60 * i_0_p # number of people invected for both scenarios
h_start <- p_60 * h_0_p # number of people hospitalized for both scenarios

# Populate period 0 in results table with initial conditions
# All unspecified compartments default to 0
# Start and end values are identical for period 0 (no transitions in initial state)
results_soc_dt <- rbindlist(list(
  results_soc_dt,
  data.table(
    period = 0L,                         # Initial period
    s_vn_start = s_vn_start,             # Initial susceptible, vaccine-naive population
    s_ve_start = 0,                      # No vaccine-experienced population at start
    v_start = 0,                         # No vaccinated population (standard of care)
    e_start = 0,                         # No exposed population at start
    i_start = i_start,                   # Initial infectious population
    h_start = h_start,                   # Initial hospitalized population
    r_start = 0,                         # No recovered population at start
    sq_start = 0,                        # No sequelae population at start
    d_start = 0,                         # No deceased population at start
    s_vn_end = s_vn_start,               # Ending values same as starting for period 0
    s_ve_end = 0,                        # Ending values same as starting for period 0
    v_end = 0,                           # Ending values same as starting for period 0
    e_end = 0,                           # Ending values same as starting for period 0
    i_end = i_start,                     # Ending values same as starting for period 0
    h_end = h_start,                     # Ending values same as starting for period 0
    r_end = 0,                           # Ending values same as starting for period 0
    sq_end = 0,                          # Ending values same as starting for period 0
    d_end = 0,                           # Ending values same as starting for period 0
    cost_total = 0,                      # No costs incurred in initial state
    utility_total = 0                    # No utility accrued in initial state
  )
), use.names = TRUE)

# ==============================================================================
# Step 3c: Simulation Loop - Standard of Care
# ==============================================================================

# Initialize loop control variables
period <- 1                          # Start from period 1 (period 0 is initial condition)
no_susceptible_count <- 0            # Counter for consecutive periods with no susceptible population
max_periods <- 20000                 # Safety limit to prevent infinite loops
susceptible_threshold <- 1           # Threshold below which susceptible population is considered depleted

# Run simulation until stopping condition is met
while (no_susceptible_count < 2 && period <= max_periods) {

  # Get previous period's ending values as current period's starting values
  prev_row <- results_soc_dt[period]

  # Starting populations for current period
  s_vn_t <- prev_row$s_vn_end        # Susceptible, vaccine-naive
  e_t <- prev_row$e_end              # Exposed
  i_t <- prev_row$i_end              # Infectious
  h_t <- prev_row$h_end              # Hospitalized
  r_t <- prev_row$r_end              # Recovered
  sq_t <- prev_row$sq_end            # Sequelae
  d_t <- prev_row$d_end              # Deceased

  # Calculate total living population for force of infection
  n_t <- s_vn_t + e_t + i_t + h_t + r_t + sq_t

  # Check if everyone is dead
  if (n_t <= 0) break

  # Calculate force of infection (transmission depends on proportion infectious)
  foi <- beta * i_t / n_t

  # ----------------------------------------------------------------------------
  # Calculate transitions from each compartment
  # For competing risks, calculate total exit rate and proportions to each destination
  # ----------------------------------------------------------------------------

  # S_VN transitions: infection (foi) or death (mu)
  total_exit_s_vn <- foi + mu
  if (total_exit_s_vn > 0) {
    n_exit_s_vn <- min(s_vn_t, s_vn_t * total_exit_s_vn)
    s_vn_to_e <- n_exit_s_vn * (foi / total_exit_s_vn)
    s_vn_to_d <- n_exit_s_vn * (mu / total_exit_s_vn)
  } else {
    s_vn_to_e <- 0
    s_vn_to_d <- 0
  }

  # E transitions: progression to infectious (sigma) or death (mu)
  total_exit_e <- sigma + mu
  if (total_exit_e > 0) {
    n_exit_e <- min(e_t, e_t * total_exit_e)
    e_to_i <- n_exit_e * (sigma / total_exit_e)
    e_to_d <- n_exit_e * (mu / total_exit_e)
  } else {
    e_to_i <- 0
    e_to_d <- 0
  }

  # I transitions: hospitalization (phi), recovery (tau_1), or death (mu + delta_m)
  total_exit_i <- phi + tau_1 + mu + delta_m
  if (total_exit_i > 0) {
    n_exit_i <- min(i_t, i_t * total_exit_i)
    i_to_h <- n_exit_i * (phi / total_exit_i)
    i_to_r <- n_exit_i * (tau_1 / total_exit_i)
    i_to_d <- n_exit_i * ((mu + delta_m) / total_exit_i)
  } else {
    i_to_h <- 0
    i_to_r <- 0
    i_to_d <- 0
  }

  # H transitions: sequelae (rho), recovery (tau_2), or death (mu + delta_s)
  total_exit_h <- rho + tau_2 + mu + delta_s
  if (total_exit_h > 0) {
    n_exit_h <- min(h_t, h_t * total_exit_h)
    h_to_sq <- n_exit_h * (rho / total_exit_h)
    h_to_r <- n_exit_h * (tau_2 / total_exit_h)
    h_to_d <- n_exit_h * ((mu + delta_s) / total_exit_h)
  } else {
    h_to_sq <- 0
    h_to_r <- 0
    h_to_d <- 0
  }

  # SQ transitions: recovery (tau_3) or death (mu + delta_sq)
  total_exit_sq <- tau_3 + mu + delta_sq
  if (total_exit_sq > 0) {
    n_exit_sq <- min(sq_t, sq_t * total_exit_sq)
    sq_to_r <- n_exit_sq * (tau_3 / total_exit_sq)
    sq_to_d <- n_exit_sq * ((mu + delta_sq) / total_exit_sq)
  } else {
    sq_to_r <- 0
    sq_to_d <- 0
  }

  # R transitions: waning immunity back to susceptible (gamma_r) or death (mu)
  total_exit_r <- gamma_r + mu
  if (total_exit_r > 0) {
    n_exit_r <- min(r_t, r_t * total_exit_r)
    r_to_s_vn <- n_exit_r * (gamma_r / total_exit_r)
    r_to_d <- n_exit_r * (mu / total_exit_r)
  } else {
    r_to_s_vn <- 0
    r_to_d <- 0
  }

  # ----------------------------------------------------------------------------
  # Calculate ending populations for current period
  # ----------------------------------------------------------------------------

  s_vn_end <- s_vn_t - s_vn_to_e - s_vn_to_d + r_to_s_vn
  e_end <- e_t + s_vn_to_e - e_to_i - e_to_d
  i_end <- i_t + e_to_i - i_to_h - i_to_r - i_to_d
  h_end <- h_t + i_to_h - h_to_sq - h_to_r - h_to_d
  r_end <- r_t + i_to_r + h_to_r + sq_to_r - r_to_s_vn - r_to_d
  sq_end <- sq_t + h_to_sq - sq_to_r - sq_to_d
  d_end <- d_t + s_vn_to_d + e_to_d + i_to_d + h_to_d + sq_to_d + r_to_d

  # ----------------------------------------------------------------------------
  # Calculate costs for current period
  # ----------------------------------------------------------------------------

  # Hospitalization costs: incurred when entering H compartment
  cost_hosp <- i_to_h * o_los_h * c_hosp

  # Indirect costs: incurred when entering I compartment (per case)
  cost_indirect <- e_to_i * (c_pl + c_care + c_hfr)

  # Death costs: incurred when entering D compartment (one-time per death)
  cost_death <- (s_vn_to_d + e_to_d + i_to_d + h_to_d + sq_to_d + r_to_d) * c_death

  # Total costs for period
  cost_total <- cost_hosp + cost_indirect + cost_death

  # ----------------------------------------------------------------------------
  # Calculate utilities for current period
  # ----------------------------------------------------------------------------

  # Weekly utility (QALYs) = (population in state * utility value) / 52 weeks per year
  # Note: V and S_VE compartments are 0 for standard of care
  utility_total <- (s_vn_t * u_s_vn +
                    e_t * u_e +
                    i_t * u_i +
                    h_t * u_h +
                    r_t * u_r +
                    sq_t * u_sq +
                    d_t * u_d) / 52

  # ----------------------------------------------------------------------------
  # Store results for current period
  # ----------------------------------------------------------------------------

  results_soc_dt <- rbindlist(list(
    results_soc_dt,
    data.table(
      period = period,
      s_vn_start = s_vn_t,
      s_ve_start = 0,                  # Not used in SOC model
      v_start = 0,                     # No vaccination in SOC
      e_start = e_t,
      i_start = i_t,
      h_start = h_t,
      r_start = r_t,
      sq_start = sq_t,
      d_start = d_t,
      s_vn_end = s_vn_end,
      s_ve_end = 0,                    # Not used in SOC model
      v_end = 0,                       # No vaccination in SOC
      e_end = e_end,
      i_end = i_end,
      h_end = h_end,
      r_end = r_end,
      sq_end = sq_end,
      d_end = d_end,
      cost_total = cost_total,
      utility_total = utility_total
    )
  ), use.names = TRUE)

  # ----------------------------------------------------------------------------
  # Real-time validation check
  # ----------------------------------------------------------------------------

  validate_period("SOC", period, results_soc_dt[period + 1], p_60,
                  if (period > 0) results_soc_dt[period]$d_end else 0)

  # ----------------------------------------------------------------------------
  # Check stopping conditions
  # ----------------------------------------------------------------------------

  # Count consecutive periods with susceptible population below threshold
  if (s_vn_end < susceptible_threshold) {
    no_susceptible_count <- no_susceptible_count + 1
  } else {
    no_susceptible_count <- 0
  }

  # Increment period counter
  period <- period + 1
}

# Print simulation summary
cat("Simulation completed after", period - 1, "periods\n")
cat("Final population distribution:\n")
cat("  Deceased:", d_end, "\n")
cat("  Living:", n_t - (d_end - d_t), "\n")

# ==============================================================================
# Step 4: Model - Arexvy Vaccine
# ==============================================================================

# ==============================================================================
# Step 4a: Initialize Results Storage
# ==============================================================================

# Initialize empty data.table to store period-by-period results for Arexvy scenario
# This table will grow dynamically as simulation progresses until no susceptible population remains
results_arexvy_dt <- data.table(
  period = integer(),                    # Period number (sequential counter)
  s_vn_start = numeric(),                # Starting susceptible, vaccine-naive population
  s_ve_start = numeric(),                # Starting susceptible, vaccine-experienced population
  v_start = numeric(),                   # Starting vaccinated population
  e_start = numeric(),                   # Starting exposed population
  i_start = numeric(),                   # Starting infectious population
  h_start = numeric(),                   # Starting hospitalized population
  r_start = numeric(),                   # Starting recovered population
  sq_start = numeric(),                  # Starting sequelae population
  d_start = numeric(),                   # Starting deceased population
  s_vn_end = numeric(),                  # Ending susceptible, vaccine-naive population
  s_ve_end = numeric(),                  # Ending susceptible, vaccine-experienced population
  v_end = numeric(),                     # Ending vaccinated population
  e_end = numeric(),                     # Ending exposed population
  i_end = numeric(),                     # Ending infectious population
  h_end = numeric(),                     # Ending hospitalized population
  r_end = numeric(),                     # Ending recovered population
  sq_end = numeric(),                    # Ending sequelae population
  d_end = numeric(),                     # Ending deceased population
  cost_total = numeric(),                # Total cost incurred during period
  utility_total = numeric()              # Total utility (QALYs) accrued during period
)

# ==============================================================================
# Step 4b: Initialize Period 0 (Starting Condition)
# ==============================================================================

# Initial Compartment Conditions for Arexvy scenario
s_vn_start_arexvy <- p_60 * s_vn_0_v_p # Number of people susceptible, vaccine-naive
v_start_arexvy <- p_60 * v_0_p         # Number of people vaccinated at time 0
i_start_arexvy <- p_60 * i_0_p         # Number of people infected
h_start_arexvy <- p_60 * h_0_p         # Number of people hospitalized

# Populate period 0 in results table with initial conditions
# All unspecified compartments default to 0
# Start and end values are identical for period 0 (no transitions in initial state)
results_arexvy_dt <- rbindlist(list(
  results_arexvy_dt,
  data.table(
    period = 0L,                         # Initial period
    s_vn_start = s_vn_start_arexvy,      # Initial susceptible, vaccine-naive population
    s_ve_start = 0,                      # No vaccine-experienced population at start
    v_start = v_start_arexvy,            # Initial vaccinated population
    e_start = 0,                         # No exposed population at start
    i_start = i_start_arexvy,            # Initial infectious population
    h_start = h_start_arexvy,            # Initial hospitalized population
    r_start = 0,                         # No recovered population at start
    sq_start = 0,                        # No sequelae population at start
    d_start = 0,                         # No deceased population at start
    s_vn_end = s_vn_start_arexvy,        # Ending values same as starting for period 0
    s_ve_end = 0,                        # Ending values same as starting for period 0
    v_end = v_start_arexvy,              # Ending values same as starting for period 0
    e_end = 0,                           # Ending values same as starting for period 0
    i_end = i_start_arexvy,              # Ending values same as starting for period 0
    h_end = h_start_arexvy,              # Ending values same as starting for period 0
    r_end = 0,                           # Ending values same as starting for period 0
    sq_end = 0,                          # Ending values same as starting for period 0
    d_end = 0,                           # Ending values same as starting for period 0
    cost_total = 0,                      # No costs incurred in initial state
    utility_total = 0                    # No utility accrued in initial state
  )
), use.names = TRUE)

# ==============================================================================
# Step 4c: Simulation Loop - Arexvy Vaccine
# ==============================================================================

# Initialize loop control variables
period <- 1                          # Start from period 1 (period 0 is initial condition)
no_susceptible_count <- 0            # Counter for consecutive periods with no susceptible population
max_periods <- 20000                 # Safety limit to prevent infinite loops
susceptible_threshold <- 1           # Threshold below which susceptible population is considered depleted

# Run simulation until stopping condition is met
while (no_susceptible_count < 2 && period <= max_periods) {

  # Get previous period's ending values as current period's starting values
  prev_row <- results_arexvy_dt[period]

  # Starting populations for current period
  s_vn_t <- prev_row$s_vn_end        # Susceptible, vaccine-naive
  s_ve_t <- prev_row$s_ve_end        # Susceptible, vaccine-experienced
  v_t <- prev_row$v_end              # Vaccinated
  e_t <- prev_row$e_end              # Exposed
  i_t <- prev_row$i_end              # Infectious
  h_t <- prev_row$h_end              # Hospitalized
  r_t <- prev_row$r_end              # Recovered
  sq_t <- prev_row$sq_end            # Sequelae
  d_t <- prev_row$d_end              # Deceased

  # Calculate total living population for force of infection
  n_t <- s_vn_t + s_ve_t + v_t + e_t + i_t + h_t + r_t + sq_t

  # Check if everyone is dead
  if (n_t <= 0) break

  # Calculate force of infection (transmission depends on proportion infectious)
  foi <- beta * i_t / n_t

  # ----------------------------------------------------------------------------
  # Calculate transitions from each compartment
  # For competing risks, calculate total exit rate and proportions to each destination
  # ----------------------------------------------------------------------------

  # S_VN transitions: vaccination (psi), infection (foi), or death (mu)
  total_exit_s_vn <- psi + foi + mu
  if (total_exit_s_vn > 0) {
    n_exit_s_vn <- min(s_vn_t, s_vn_t * total_exit_s_vn)
    s_vn_to_v <- n_exit_s_vn * (psi / total_exit_s_vn)
    s_vn_to_e <- n_exit_s_vn * (foi / total_exit_s_vn)
    s_vn_to_d <- n_exit_s_vn * (mu / total_exit_s_vn)
  } else {
    s_vn_to_v <- 0
    s_vn_to_e <- 0
    s_vn_to_d <- 0
  }

  # V transitions: waning immunity (gamma_v), breakthrough infection (kappa_arexvy), or death (mu)
  total_exit_v <- gamma_v + kappa_arexvy + mu
  if (total_exit_v > 0) {
    n_exit_v <- min(v_t, v_t * total_exit_v)
    v_to_s_ve <- n_exit_v * (gamma_v / total_exit_v)
    v_to_e <- n_exit_v * (kappa_arexvy / total_exit_v)
    v_to_d <- n_exit_v * (mu / total_exit_v)
  } else {
    v_to_s_ve <- 0
    v_to_e <- 0
    v_to_d <- 0
  }

  # S_VE transitions: infection (foi) or death (mu)
  total_exit_s_ve <- foi + mu
  if (total_exit_s_ve > 0) {
    n_exit_s_ve <- min(s_ve_t, s_ve_t * total_exit_s_ve)
    s_ve_to_e <- n_exit_s_ve * (foi / total_exit_s_ve)
    s_ve_to_d <- n_exit_s_ve * (mu / total_exit_s_ve)
  } else {
    s_ve_to_e <- 0
    s_ve_to_d <- 0
  }

  # E transitions: progression to infectious (sigma) or death (mu)
  total_exit_e <- sigma + mu
  if (total_exit_e > 0) {
    n_exit_e <- min(e_t, e_t * total_exit_e)
    e_to_i <- n_exit_e * (sigma / total_exit_e)
    e_to_d <- n_exit_e * (mu / total_exit_e)
  } else {
    e_to_i <- 0
    e_to_d <- 0
  }

  # I transitions: hospitalization (phi), recovery (tau_1), or death (mu + delta_m)
  total_exit_i <- phi + tau_1 + mu + delta_m
  if (total_exit_i > 0) {
    n_exit_i <- min(i_t, i_t * total_exit_i)
    i_to_h <- n_exit_i * (phi / total_exit_i)
    i_to_r <- n_exit_i * (tau_1 / total_exit_i)
    i_to_d <- n_exit_i * ((mu + delta_m) / total_exit_i)
  } else {
    i_to_h <- 0
    i_to_r <- 0
    i_to_d <- 0
  }

  # H transitions: sequelae (rho), recovery (tau_2), or death (mu + delta_s)
  total_exit_h <- rho + tau_2 + mu + delta_s
  if (total_exit_h > 0) {
    n_exit_h <- min(h_t, h_t * total_exit_h)
    h_to_sq <- n_exit_h * (rho / total_exit_h)
    h_to_r <- n_exit_h * (tau_2 / total_exit_h)
    h_to_d <- n_exit_h * ((mu + delta_s) / total_exit_h)
  } else {
    h_to_sq <- 0
    h_to_r <- 0
    h_to_d <- 0
  }

  # SQ transitions: recovery (tau_3) or death (mu + delta_sq)
  total_exit_sq <- tau_3 + mu + delta_sq
  if (total_exit_sq > 0) {
    n_exit_sq <- min(sq_t, sq_t * total_exit_sq)
    sq_to_r <- n_exit_sq * (tau_3 / total_exit_sq)
    sq_to_d <- n_exit_sq * ((mu + delta_sq) / total_exit_sq)
  } else {
    sq_to_r <- 0
    sq_to_d <- 0
  }

  # R transitions: waning immunity back to susceptible (gamma_r) or death (mu)
  total_exit_r <- gamma_r + mu
  if (total_exit_r > 0) {
    n_exit_r <- min(r_t, r_t * total_exit_r)
    r_to_s_vn <- n_exit_r * (gamma_r / total_exit_r)
    r_to_d <- n_exit_r * (mu / total_exit_r)
  } else {
    r_to_s_vn <- 0
    r_to_d <- 0
  }

  # ----------------------------------------------------------------------------
  # Calculate ending populations for current period
  # ----------------------------------------------------------------------------

  s_vn_end <- s_vn_t - s_vn_to_v - s_vn_to_e - s_vn_to_d + r_to_s_vn
  s_ve_end <- s_ve_t + v_to_s_ve - s_ve_to_e - s_ve_to_d
  v_end <- v_t + s_vn_to_v - v_to_s_ve - v_to_e - v_to_d
  e_end <- e_t + s_vn_to_e + s_ve_to_e + v_to_e - e_to_i - e_to_d
  i_end <- i_t + e_to_i - i_to_h - i_to_r - i_to_d
  h_end <- h_t + i_to_h - h_to_sq - h_to_r - h_to_d
  r_end <- r_t + i_to_r + h_to_r + sq_to_r - r_to_s_vn - r_to_d
  sq_end <- sq_t + h_to_sq - sq_to_r - sq_to_d
  d_end <- d_t + s_vn_to_d + s_ve_to_d + v_to_d + e_to_d + i_to_d + h_to_d + sq_to_d + r_to_d

  # ----------------------------------------------------------------------------
  # Calculate costs for current period
  # ----------------------------------------------------------------------------

  # Vaccination costs: incurred when transitioning from S_VN to V
  cost_vaccine <- s_vn_to_v * (c_arexvy + c_vac_admin)

  # Hospitalization costs: incurred when entering H compartment
  cost_hosp <- i_to_h * o_los_h * c_hosp

  # Indirect costs: incurred when entering I compartment (per case)
  cost_indirect <- e_to_i * (c_pl + c_care + c_hfr)

  # Death costs: incurred when entering D compartment (one-time per death)
  cost_death <- (s_vn_to_d + s_ve_to_d + v_to_d + e_to_d + i_to_d + h_to_d + sq_to_d + r_to_d) * c_death

  # Total costs for period
  cost_total <- cost_vaccine + cost_hosp + cost_indirect + cost_death

  # ----------------------------------------------------------------------------
  # Calculate utilities for current period
  # ----------------------------------------------------------------------------

  # Ongoing state utilities (people currently in each state)
  utility_states <- (s_vn_t * u_s_vn +
                     s_ve_t * u_s_ve +
                     v_t * u_v +
                     e_t * u_e +
                     i_t * u_i +
                     h_t * u_h +
                     r_t * u_r +
                     sq_t * u_sq +
                     d_t * u_d) / 52

  # One-time disutility for NEW vaccinations (adverse events)
  utility_adverse_events <- s_vn_to_v * u_adverse_events * (duration_adverse_events / 52)

  # Total utility for period
  utility_total <- utility_states + utility_adverse_events

  # ----------------------------------------------------------------------------
  # Store results for current period
  # ----------------------------------------------------------------------------

  results_arexvy_dt <- rbindlist(list(
    results_arexvy_dt,
    data.table(
      period = period,
      s_vn_start = s_vn_t,
      s_ve_start = s_ve_t,
      v_start = v_t,
      e_start = e_t,
      i_start = i_t,
      h_start = h_t,
      r_start = r_t,
      sq_start = sq_t,
      d_start = d_t,
      s_vn_end = s_vn_end,
      s_ve_end = s_ve_end,
      v_end = v_end,
      e_end = e_end,
      i_end = i_end,
      h_end = h_end,
      r_end = r_end,
      sq_end = sq_end,
      d_end = d_end,
      cost_total = cost_total,
      utility_total = utility_total
    )
  ), use.names = TRUE)

  # ----------------------------------------------------------------------------
  # Real-time validation check
  # ----------------------------------------------------------------------------

  validate_period("Arexvy", period, results_arexvy_dt[period + 1], p_60,
                  if (period > 0) results_arexvy_dt[period]$d_end else 0)

  # ----------------------------------------------------------------------------
  # Check stopping conditions
  # ----------------------------------------------------------------------------

  # Count consecutive periods with susceptible population below threshold
  if (s_vn_end < susceptible_threshold) {
    no_susceptible_count <- no_susceptible_count + 1
  } else {
    no_susceptible_count <- 0
  }

  # Increment period counter
  period <- period + 1
}

# Print simulation summary
cat("Arexvy simulation completed after", period - 1, "periods\n")
cat("Final population distribution:\n")
cat("  Deceased:", d_end, "\n")
cat("  Living:", n_t - (d_end - d_t), "\n")

# ==============================================================================
# Step 5: Model - Abrysvo Vaccine
# ==============================================================================

# ==============================================================================
# Step 5a: Initialize Results Storage
# ==============================================================================

# Initialize empty data.table to store period-by-period results for Abrysvo scenario
# This table will grow dynamically as simulation progresses until no susceptible population remains
results_abrysvo_dt <- data.table(
  period = integer(),                    # Period number (sequential counter)
  s_vn_start = numeric(),                # Starting susceptible, vaccine-naive population
  s_ve_start = numeric(),                # Starting susceptible, vaccine-experienced population
  v_start = numeric(),                   # Starting vaccinated population
  e_start = numeric(),                   # Starting exposed population
  i_start = numeric(),                   # Starting infectious population
  h_start = numeric(),                   # Starting hospitalized population
  r_start = numeric(),                   # Starting recovered population
  sq_start = numeric(),                  # Starting sequelae population
  d_start = numeric(),                   # Starting deceased population
  s_vn_end = numeric(),                  # Ending susceptible, vaccine-naive population
  s_ve_end = numeric(),                  # Ending susceptible, vaccine-experienced population
  v_end = numeric(),                     # Ending vaccinated population
  e_end = numeric(),                     # Ending exposed population
  i_end = numeric(),                     # Ending infectious population
  h_end = numeric(),                     # Ending hospitalized population
  r_end = numeric(),                     # Ending recovered population
  sq_end = numeric(),                    # Ending sequelae population
  d_end = numeric(),                     # Ending deceased population
  cost_total = numeric(),                # Total cost incurred during period
  utility_total = numeric()              # Total utility (QALYs) accrued during period
)

# ==============================================================================
# Step 5b: Initialize Period 0 (Starting Condition)
# ==============================================================================

# Initial Compartment Conditions for Abrysvo scenario
s_vn_start_abrysvo <- p_60 * s_vn_0_v_p # Number of people susceptible, vaccine-naive
v_start_abrysvo <- p_60 * v_0_p         # Number of people vaccinated at time 0
i_start_abrysvo <- p_60 * i_0_p         # Number of people infected
h_start_abrysvo <- p_60 * h_0_p         # Number of people hospitalized

# Populate period 0 in results table with initial conditions
# All unspecified compartments default to 0
# Start and end values are identical for period 0 (no transitions in initial state)
results_abrysvo_dt <- rbindlist(list(
  results_abrysvo_dt,
  data.table(
    period = 0L,                         # Initial period
    s_vn_start = s_vn_start_abrysvo,     # Initial susceptible, vaccine-naive population
    s_ve_start = 0,                      # No vaccine-experienced population at start
    v_start = v_start_abrysvo,           # Initial vaccinated population
    e_start = 0,                         # No exposed population at start
    i_start = i_start_abrysvo,           # Initial infectious population
    h_start = h_start_abrysvo,           # Initial hospitalized population
    r_start = 0,                         # No recovered population at start
    sq_start = 0,                        # No sequelae population at start
    d_start = 0,                         # No deceased population at start
    s_vn_end = s_vn_start_abrysvo,       # Ending values same as starting for period 0
    s_ve_end = 0,                        # Ending values same as starting for period 0
    v_end = v_start_abrysvo,             # Ending values same as starting for period 0
    e_end = 0,                           # Ending values same as starting for period 0
    i_end = i_start_abrysvo,             # Ending values same as starting for period 0
    h_end = h_start_abrysvo,             # Ending values same as starting for period 0
    r_end = 0,                           # Ending values same as starting for period 0
    sq_end = 0,                          # Ending values same as starting for period 0
    d_end = 0,                           # Ending values same as starting for period 0
    cost_total = 0,                      # No costs incurred in initial state
    utility_total = 0                    # No utility accrued in initial state
  )
), use.names = TRUE)

# ==============================================================================
# Step 5c: Simulation Loop - Abrysvo Vaccine
# ==============================================================================

# Initialize loop control variables
period <- 1                          # Start from period 1 (period 0 is initial condition)
no_susceptible_count <- 0            # Counter for consecutive periods with no susceptible population
max_periods <- 20000                 # Safety limit to prevent infinite loops
susceptible_threshold <- 1           # Threshold below which susceptible population is considered depleted

# Run simulation until stopping condition is met
while (no_susceptible_count < 2 && period <= max_periods) {

  # Get previous period's ending values as current period's starting values
  prev_row <- results_abrysvo_dt[period]

  # Starting populations for current period
  s_vn_t <- prev_row$s_vn_end        # Susceptible, vaccine-naive
  s_ve_t <- prev_row$s_ve_end        # Susceptible, vaccine-experienced
  v_t <- prev_row$v_end              # Vaccinated
  e_t <- prev_row$e_end              # Exposed
  i_t <- prev_row$i_end              # Infectious
  h_t <- prev_row$h_end              # Hospitalized
  r_t <- prev_row$r_end              # Recovered
  sq_t <- prev_row$sq_end            # Sequelae
  d_t <- prev_row$d_end              # Deceased

  # Calculate total living population for force of infection
  n_t <- s_vn_t + s_ve_t + v_t + e_t + i_t + h_t + r_t + sq_t

  # Check if everyone is dead
  if (n_t <= 0) break

  # Calculate force of infection (transmission depends on proportion infectious)
  foi <- beta * i_t / n_t

  # ----------------------------------------------------------------------------
  # Calculate transitions from each compartment
  # For competing risks, calculate total exit rate and proportions to each destination
  # ----------------------------------------------------------------------------

  # S_VN transitions: vaccination (psi), infection (foi), or death (mu)
  total_exit_s_vn <- psi + foi + mu
  if (total_exit_s_vn > 0) {
    n_exit_s_vn <- min(s_vn_t, s_vn_t * total_exit_s_vn)
    s_vn_to_v <- n_exit_s_vn * (psi / total_exit_s_vn)
    s_vn_to_e <- n_exit_s_vn * (foi / total_exit_s_vn)
    s_vn_to_d <- n_exit_s_vn * (mu / total_exit_s_vn)
  } else {
    s_vn_to_v <- 0
    s_vn_to_e <- 0
    s_vn_to_d <- 0
  }

  # V transitions: waning immunity (gamma_v), breakthrough infection (kappa_abrysvo), or death (mu)
  total_exit_v <- gamma_v + kappa_abrysvo + mu
  if (total_exit_v > 0) {
    n_exit_v <- min(v_t, v_t * total_exit_v)
    v_to_s_ve <- n_exit_v * (gamma_v / total_exit_v)
    v_to_e <- n_exit_v * (kappa_abrysvo / total_exit_v)
    v_to_d <- n_exit_v * (mu / total_exit_v)
  } else {
    v_to_s_ve <- 0
    v_to_e <- 0
    v_to_d <- 0
  }

  # S_VE transitions: infection (foi) or death (mu)
  total_exit_s_ve <- foi + mu
  if (total_exit_s_ve > 0) {
    n_exit_s_ve <- min(s_ve_t, s_ve_t * total_exit_s_ve)
    s_ve_to_e <- n_exit_s_ve * (foi / total_exit_s_ve)
    s_ve_to_d <- n_exit_s_ve * (mu / total_exit_s_ve)
  } else {
    s_ve_to_e <- 0
    s_ve_to_d <- 0
  }

  # E transitions: progression to infectious (sigma) or death (mu)
  total_exit_e <- sigma + mu
  if (total_exit_e > 0) {
    n_exit_e <- min(e_t, e_t * total_exit_e)
    e_to_i <- n_exit_e * (sigma / total_exit_e)
    e_to_d <- n_exit_e * (mu / total_exit_e)
  } else {
    e_to_i <- 0
    e_to_d <- 0
  }

  # I transitions: hospitalization (phi), recovery (tau_1), or death (mu + delta_m)
  total_exit_i <- phi + tau_1 + mu + delta_m
  if (total_exit_i > 0) {
    n_exit_i <- min(i_t, i_t * total_exit_i)
    i_to_h <- n_exit_i * (phi / total_exit_i)
    i_to_r <- n_exit_i * (tau_1 / total_exit_i)
    i_to_d <- n_exit_i * ((mu + delta_m) / total_exit_i)
  } else {
    i_to_h <- 0
    i_to_r <- 0
    i_to_d <- 0
  }

  # H transitions: sequelae (rho), recovery (tau_2), or death (mu + delta_s)
  total_exit_h <- rho + tau_2 + mu + delta_s
  if (total_exit_h > 0) {
    n_exit_h <- min(h_t, h_t * total_exit_h)
    h_to_sq <- n_exit_h * (rho / total_exit_h)
    h_to_r <- n_exit_h * (tau_2 / total_exit_h)
    h_to_d <- n_exit_h * ((mu + delta_s) / total_exit_h)
  } else {
    h_to_sq <- 0
    h_to_r <- 0
    h_to_d <- 0
  }

  # SQ transitions: recovery (tau_3) or death (mu + delta_sq)
  total_exit_sq <- tau_3 + mu + delta_sq
  if (total_exit_sq > 0) {
    n_exit_sq <- min(sq_t, sq_t * total_exit_sq)
    sq_to_r <- n_exit_sq * (tau_3 / total_exit_sq)
    sq_to_d <- n_exit_sq * ((mu + delta_sq) / total_exit_sq)
  } else {
    sq_to_r <- 0
    sq_to_d <- 0
  }

  # R transitions: waning immunity back to susceptible (gamma_r) or death (mu)
  total_exit_r <- gamma_r + mu
  if (total_exit_r > 0) {
    n_exit_r <- min(r_t, r_t * total_exit_r)
    r_to_s_vn <- n_exit_r * (gamma_r / total_exit_r)
    r_to_d <- n_exit_r * (mu / total_exit_r)
  } else {
    r_to_s_vn <- 0
    r_to_d <- 0
  }

  # ----------------------------------------------------------------------------
  # Calculate ending populations for current period
  # ----------------------------------------------------------------------------

  s_vn_end <- s_vn_t - s_vn_to_v - s_vn_to_e - s_vn_to_d + r_to_s_vn
  s_ve_end <- s_ve_t + v_to_s_ve - s_ve_to_e - s_ve_to_d
  v_end <- v_t + s_vn_to_v - v_to_s_ve - v_to_e - v_to_d
  e_end <- e_t + s_vn_to_e + s_ve_to_e + v_to_e - e_to_i - e_to_d
  i_end <- i_t + e_to_i - i_to_h - i_to_r - i_to_d
  h_end <- h_t + i_to_h - h_to_sq - h_to_r - h_to_d
  r_end <- r_t + i_to_r + h_to_r + sq_to_r - r_to_s_vn - r_to_d
  sq_end <- sq_t + h_to_sq - sq_to_r - sq_to_d
  d_end <- d_t + s_vn_to_d + s_ve_to_d + v_to_d + e_to_d + i_to_d + h_to_d + sq_to_d + r_to_d

  # ----------------------------------------------------------------------------
  # Calculate costs for current period
  # ----------------------------------------------------------------------------

  # Vaccination costs: incurred when transitioning from S_VN to V
  cost_vaccine <- s_vn_to_v * (c_abrysvo + c_vac_admin)

  # Hospitalization costs: incurred when entering H compartment
  cost_hosp <- i_to_h * o_los_h * c_hosp

  # Indirect costs: incurred when entering I compartment (per case)
  cost_indirect <- e_to_i * (c_pl + c_care + c_hfr)

  # Death costs: incurred when entering D compartment (one-time per death)
  cost_death <- (s_vn_to_d + s_ve_to_d + v_to_d + e_to_d + i_to_d + h_to_d + sq_to_d + r_to_d) * c_death

  # Total costs for period
  cost_total <- cost_vaccine + cost_hosp + cost_indirect + cost_death

  # ----------------------------------------------------------------------------
  # Calculate utilities for current period
  # ----------------------------------------------------------------------------

  # Ongoing state utilities (people currently in each state)
  utility_states <- (s_vn_t * u_s_vn +
                     s_ve_t * u_s_ve +
                     v_t * u_v +
                     e_t * u_e +
                     i_t * u_i +
                     h_t * u_h +
                     r_t * u_r +
                     sq_t * u_sq +
                     d_t * u_d) / 52

  # One-time disutility for NEW vaccinations (adverse events)
  utility_adverse_events <- s_vn_to_v * u_adverse_events * (duration_adverse_events / 52)

  # Total utility for period
  utility_total <- utility_states + utility_adverse_events

  # ----------------------------------------------------------------------------
  # Store results for current period
  # ----------------------------------------------------------------------------

  results_abrysvo_dt <- rbindlist(list(
    results_abrysvo_dt,
    data.table(
      period = period,
      s_vn_start = s_vn_t,
      s_ve_start = s_ve_t,
      v_start = v_t,
      e_start = e_t,
      i_start = i_t,
      h_start = h_t,
      r_start = r_t,
      sq_start = sq_t,
      d_start = d_t,
      s_vn_end = s_vn_end,
      s_ve_end = s_ve_end,
      v_end = v_end,
      e_end = e_end,
      i_end = i_end,
      h_end = h_end,
      r_end = r_end,
      sq_end = sq_end,
      d_end = d_end,
      cost_total = cost_total,
      utility_total = utility_total
    )
  ), use.names = TRUE)

  # ----------------------------------------------------------------------------
  # Real-time validation check
  # ----------------------------------------------------------------------------

  validate_period("Abrysvo", period, results_abrysvo_dt[period + 1], p_60,
                  if (period > 0) results_abrysvo_dt[period]$d_end else 0)

  # ----------------------------------------------------------------------------
  # Check stopping conditions
  # ----------------------------------------------------------------------------

  # Count consecutive periods with susceptible population below threshold
  if (s_vn_end < susceptible_threshold) {
    no_susceptible_count <- no_susceptible_count + 1
  } else {
    no_susceptible_count <- 0
  }

  # Increment period counter
  period <- period + 1
}

# Print simulation summary
cat("Abrysvo simulation completed after", period - 1, "periods\n")
cat("Final population distribution:\n")
cat("  Deceased:", d_end, "\n")
cat("  Living:", n_t - (d_end - d_t), "\n")

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

# Calculate ICER for both vaccines vs SOC to determine row ordering
inc_cost_arexvy_vs_soc <- summary_arexvy$total_cost - summary_soc$total_cost
inc_qaly_arexvy_vs_soc <- summary_arexvy$total_qaly - summary_soc$total_qaly
icer_arexvy_vs_soc <- inc_cost_arexvy_vs_soc / inc_qaly_arexvy_vs_soc

inc_cost_abrysvo_vs_soc <- summary_abrysvo$total_cost - summary_soc$total_cost
inc_qaly_abrysvo_vs_soc <- summary_abrysvo$total_qaly - summary_soc$total_qaly
icer_abrysvo_vs_soc <- inc_cost_abrysvo_vs_soc / inc_qaly_abrysvo_vs_soc

# Determine which vaccine has lower ICER vs SOC (goes in row 2)
if (icer_arexvy_vs_soc < icer_abrysvo_vs_soc) {
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

# Calculate incrementals for row 3 (compared to row 2, not SOC)
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
cat("      Row 2 compared to Standard of Care\n")
cat("      Row 3 compared to Row 2\n")
cat("\n")

# ==============================================================================
# Step 7: Post-Simulation Validation Report
# ==============================================================================

# Run comprehensive validation report
validation_results <- generate_validation_report()

