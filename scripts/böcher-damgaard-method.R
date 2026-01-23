##### böcher damgaard method


# importing test data
library(tidyverse)
böcher_test_data <- read_excel("~/Library/CloudStorage/OneDrive-Aarhusuniversitet/MappingPlants/01 Vegetation changes Kobbefjord/böcher-method/böcher_test_data.xlsx")

df <- böcher_test_data |> 
  pivot_longer(cols = 2:21, names_to = "plot", values_to = "score") |> 
  select(species,plot, score)

head(df)

# --------------------------------------------------------------------
# Implementation of Damgaard (2015) Böcher-modified Raunkiær model
# Annotated by correspondence to Eqs. (1–6) in the paper
# --------------------------------------------------------------------

# ---------- PARAMETERS ----------
d <- c(1, 10, 25, 100)    # circle sizes in units of smallest subplot
n <- max(d)               # total number of smallest subplots in largest circle = 100

# ---------- (Eq. 1) Beta–binomial distribution for number of occupied subplots ----------
# Equation (1) in Damgaard (2015):
#   P(R = r | p, δ) = choose(n, r) * B(r + α, n - r + β) / B(α, β)
# where α = p(1−δ)/δ and β = (1−p)(1−δ)/δ
dbetabinom <- function(r, n, p, delta, log = FALSE) {
  if (delta <= 0) { # special case δ→0 ⇒ simple binomial (random)
    val <- dbinom(r, size = n, prob = p, log = log)
    return(val)
  }
  alpha <- (p * (1 - delta)) / delta   # Eq. (2) parameterization
  beta  <- ((1 - p) * (1 - delta)) / delta
  logpmf <- lchoose(n, r) + lbeta(r + alpha, n - r + beta) - lbeta(alpha, beta)
  if (log) return(logpmf)
  exp(logpmf)
}

# ---------- (Eqs. 3–4) Conditional probability of category y given r ----------
# Eq. (3): Probability that a species occurs in at least one of the d smallest subplots
# Eq. (4): Probabilities for the annuli (categories) derived from differences between circles
prob_y_given_r <- function(y, r, n, d) {
  if (r == 0) {
    if (y == 0) return(1)
    return(0)
  }
  
  prob_present_in_circle <- function(dsize) {
    1 - choose(n - dsize, r) / choose(n, r)
  }
  
  if (y == 4) {
    # present in smallest circle (d[1])
    return(prob_present_in_circle(d[1]))
  } else if (y == 3) {
    # present in second circle but not in smallest
    return(prob_present_in_circle(d[2]) - prob_present_in_circle(d[1]))
  } else if (y == 2) {
    # present in third circle but not smaller
    return(prob_present_in_circle(d[3]) - prob_present_in_circle(d[2]))
  } else if (y == 1) {
    # present in largest but not smaller
    return(prob_present_in_circle(d[4]) - prob_present_in_circle(d[3]))
  } else if (y == 0) {
    # absent from all circles
    return(choose(n - d[4], r) / choose(n, r))
  } else {
    stop("Invalid y. Use 0,1,2,3,4.")
  }
}

# ---------- (Eq. 5) Marginal probability of observing category y ----------
# Eq. (5): P(y | p, δ) = Σ_r [ P(y | r) * P(r | p, δ) ]
prob_y_given_pdelta <- function(y, p, delta, n, d) {
  prs <- sapply(0:n, function(r) prob_y_given_r(y, r, n, d))
  pmf_r <- sapply(0:n, function(r) dbetabinom(r, n, p, delta))
  sum(prs * pmf_r)
}

# ---------- (Eq. 6) Log-likelihood for a vector of independent plots ----------
# Eq. (6): logL(p, δ | y₁,…,yₙ) = Σ_i log P(y_i | p, δ)
logLik_pdelta <- function(params, y_obs, n, d) {
  logit_p <- params[1]; logit_delta <- params[2]
  p <- 1 / (1 + exp(-logit_p))
  delta <- 1 / (1 + exp(-logit_delta))
  delta <- pmax(delta, 1e-9)
  p <- pmin(pmax(p, 1e-9), 1-1e-9)
  logps <- sapply(y_obs, function(y) {
    pr <- prob_y_given_pdelta(y, p, delta, n, d)
    if (pr <= 0) return(-1e12)
    log(pr)
  })
  sum(logps)
}

# ---------- Numerical maximization of Eq. (6) ----------
fit_mle <- function(y_obs, n, d, start = c(0, 0)) {
  res <- optim(par = start,
               fn = function(par) -logLik_pdelta(par, y_obs, n, d),
               method = "Nelder-Mead",
               control = list(maxit = 5e4))
  p_hat <- 1 / (1 + exp(-res$par[1]))
  delta_hat <- 1 / (1 + exp(-res$par[2]))
  list(par = res$par, p = p_hat, delta = delta_hat, value = res$value, conv = res$convergence)
}

# ---------- Simple Metropolis sampler (no direct paper equation; exploratory extension) ----------
# This section is not in Damgaard (2015), but provides a Bayesian way
# to sample from posterior(p, δ | y) ∝ L(p, δ) under flat priors.
metropolis_simple <- function(y_obs, n, d, niter = 20000, start = c(0.2, 0.2), prop_sd = c(0.1, 0.1)) {
  chain <- matrix(NA, nrow = niter, ncol = 2)
  cur <- start
  cur_ll <- logLik_pdelta(cur, y_obs, n, d)
  accept <- 0
  for (i in 1:niter) {
    prop <- cur + rnorm(2, mean = 0, sd = prop_sd)
    prop_ll <- logLik_pdelta(prop, y_obs, n, d)
    if (log(runif(1)) < (prop_ll - cur_ll)) {
      cur <- prop; cur_ll <- prop_ll; accept <- accept + 1
    }
    chain[i, ] <- cur
  }
  list(chain = chain, accept_rate = accept / niter)
}

# ---------- Example usage ----------
set.seed(42)
y_obs <- sample(c(0,1,2,3,4), size = 30, replace = TRUE, prob = c(0.2, 0.2, 0.2, 0.2, 0.2))
start_logit <- qlogis(c(0.2, 0.2))
#fit_mle(y_obs, n, d, start = qlogis(c(0.7, 0.7))) to test to get app. the same delta and p
#fit_mle(y_obs, n, d, start = qlogis(c(0.05, 0.9))) to test to get app. the same delta and p
mle <- fit_mle(y_obs, n, d, start = start_logit)
mle

# ---------- Example usage ----------
results <- df |> 
  group_by(species)  |> 
  summarize(fit = list(fit_mle(score, n, d, start = qlogis(c(0.2, 0.2))))) |> 
  mutate(
    p = sapply(fit, `[[`, "p"),
    delta = sapply(fit, `[[`, "delta")
  )

joined_results <- böcher_test_data |>
  left_join(results, by = "species")

ggplot(joined_results, aes(x = p.x, ))  