#' Sample Size Calculation for Matched Logistic Regression
#' 
#' @param alpha Significance level
#' @param power Desired power (optional if n is provided)
#' @param n Sample size (optional if power is provided)
#' @param OR Odds Ratio to detect
#' @param p_discordant Proportion of discordant pairs
#' @param m Number of controls per case (default 1 for 1:1 matching)
#' @param method Calculation method: "rho" or "se" (standard error)
#' @param rho_values Vector of rho values for sample size calculation
#' @param SE Standard error of coefficient (not of OR)
#' @param CI_upper Upper confidence interval limit of OR
#' @param CI_lower Lower confidence interval limit of OR
#' @param n_previous Sample size of previous study (required for method "se")
#' @return A data frame with results according to chosen method
#' @export
SampleCCMatchedMultiENG <- function(alpha, power = NULL, n = NULL,
                                   OR = NULL,
                                   p_discordant,
                                   m = 1,
                                   method = "rho",
                                   rho_values = seq(0, 0.9, by = 0.1),
                                   SE = NULL,
                                   CI_upper = NULL, 
                                   CI_lower = NULL,
                                   n_previous = NULL) {
  
  if (method == "se") {
    if (is.null(SE) && (is.null(CI_upper) || is.null(CI_lower))) {
      stop("For 'se' method, SE or both CI limits are required")
    }
    if (is.null(n_previous)) {
      stop("For 'se' method, n_previous is required")
    }
  }
  
  # If CI provided, calculate SE
  if (!is.null(CI_upper) && !is.null(CI_lower)) {
    SE <- (log(CI_upper) - log(OR)) / qnorm(0.975)
    cat("Calculated coefficient Standard Error:", SE, "
")
  }
  
  z_alpha <- qnorm(1 - alpha/2)
  
  if (method == "rho") {
    if (!is.null(power)) {
      z_beta <- qnorm(power)
      
      calculate_n <- function(rho) {
        numerator <- (z_alpha * sqrt((m + 1) / m) + z_beta)^2
        denominator <- p_discordant * (log(OR))^2 * (1 - rho^2)
        
        n_pairs <- ceiling(numerator / denominator)
        return(n_pairs)
      }
      
      results <- data.frame(
        rho = rho_values,
        n_pairs = sapply(rho_values, calculate_n)
      )
      
      results$total_size <- results$n_pairs * (m + 1)
      
    } else {
      # Power calculation
      calculate_power <- function(rho) {
        z_beta <- sqrt(n * p_discordant * (log(OR))^2 * (1 - rho^2)) / 
          sqrt((m + 1) / m) - z_alpha
        power <- pnorm(z_beta)
        return(power)
      }
      
      results <- data.frame(
        rho = rho_values,
        power = sapply(rho_values, calculate_power)
      )
    }
  } else {
    # Standard error based method
    if (!is.null(power)) {
      z_gamma <- qnorm(power)
      n <- ceiling((z_alpha + z_gamma)^2 * n_previous * SE^2 / (log(OR))^2)
      results <- data.frame(
        n_pairs = n
      )
      results$total_size <- n * (m + 1)
    } else {
      z_gamma <- sqrt(n * (log(OR))^2 / (n_previous * SE^2)) - z_alpha
      power <- pnorm(z_gamma)
      results <- data.frame(
        power = power,
        n_pairs = n,
        total_size = n * (m + 1)
      )
    }
  }
  
  return(results)
}

#' Logistics for Matched Study
#'
#' @param n_pairs Number of required pairs
#' @param m Number of controls per case
#' @param pair_identification_rate Rate of pair identification per day
#' @param pair_rejection_rate Expected rejection rate for pairs
#' @param working_days_month Number of working days per month
#' @return A list with study logistics calculations
#' @export
matched_study_logistics <- function(n_pairs,
                                  m = 1,
                                  pair_identification_rate,
                                  pair_rejection_rate,
                                  working_days_month) {
  
  pairs_to_contact <- n_pairs / (1 - pair_rejection_rate)
  
  total_days <- ceiling(pairs_to_contact / pair_identification_rate)
  total_months <- total_days / working_days_month
  
  results <- list(
    required_pairs = n_pairs,
    total_participants = n_pairs * (m + 1),
    pairs_to_contact = ceiling(pairs_to_contact),
    recruitment_days = total_days,
    recruitment_months = round(total_months, 2)
  )
  
  # Print summary
  cat("
Matched study logistics summary:
")
  cat("----------------------------------------
")
  cat("Required pairs:", n_pairs, "
")
  cat("Total participants:", n_pairs * (m + 1),
      "(", m, ":", 1, "matching)
")
  cat("Pairs to contact:", ceiling(pairs_to_contact),
      "(", pair_rejection_rate*100, "% rejection)
")
  cat("Days needed:", total_days,
      "(", pair_identification_rate, "pairs per day)
")
  cat("Months needed:", round(total_months, 2),
      "(", working_days_month, "working days per month)
")
  cat("----------------------------------------
")
  
  return(invisible(results))
}
