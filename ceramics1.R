library(stats4)

# Data for Ceramic Phase 1
phase1 <- data.frame(
  x = c(2530, 2420, 2160, 2770, 2370, 2440, 2330, 2300, 2460, 2210, 2450),
  s = c(110, 80, 80, 70, 80, 70, 60, 100, 60, 60, 80)
)

# Data for Ceramic Phase 2
phase2 <- data.frame(
  x = c(2290, 2330, 2340, 2270, 2140, 2300, 2120, 2580, 2180),
  s = c(60, 90, 100, 90, 80, 90, 70, 80, 90)
)

# Data for Ceramic Phase 3
phase3 <- data.frame(
  x = c(2230, 2060, 2210, 2120, 2210, 2470, 2520, 2220, 2210, 2090, 2090, 2210, 2250, 2280, 2300, 2330, 2110),
  s = c(70, 80, 70, 80, 70, 90, 80, 70, 70, 90, 70, 70, 70, 80, 70, 70, 70)
)

# Data for Ceramic Phase 4
phase4 <- data.frame(
  x = c(2140, 2030, 2100, 2110, 2060, 1990, 2170, 2040, 2160, 2200, 2300, 2060, 2170, 2370, 2120, 2150, 1980, 2260, 2090, 2120, 2000, 2130, 1900),
  s = c(80, 80, 90, 80, 60, 70, 70, 70, 70, 80, 70, 70, 70, 70, 70, 70, 80, 80, 60, 70, 70, 60, 60)
)

# calibration curve data
calibration_data <- data.frame(
  calendar = c(1735, 1765, 1830, 1890, 1950, 2045, 2110, 2155, 2165, 2175, 2270, 2285, 2190, 2250, 2305, 2230, 2315, 2320, 2350, 2375, 2400, 2440, 2480, 2550, 2705, 2750, 2790, 2830, 2875, 2925, 2980, 3050, 3125, 2510, 2530, 2565, 2575, 2620, 2655, 2585, 2605, 2665, 2680),
  radiocarbon = c(1800, 1850, 1900, 1950, 2000, 2050, 2100, 2150, 2160, 2170, 2170, 2170, 2180, 2180, 2180, 2190, 2190, 2200, 2250, 2300, 2350, 2400, 2420, 2450, 2500, 2550, 2600, 2650, 2700, 2750, 2800, 2850, 2900, 2430, 2440, 2460, 2465, 2465, 2465, 2470, 2470, 2470, 2480)
)


calibration_data <- calibration_data[order(calibration_data$calendar), ]


# Define the piecewise linear function mu(theta)
mu_theta <- function(theta) {
  # get known calendar years and corresponding radiocarbon dates
  knots.cal <- calibration_data$calendar
  values.cal <- calibration_data$radiocarbon
  
  # find interval that has theta
  k_map=function(theta){
  if (theta <= min(knots.cal)) {
    k <- which.min(knots.cal) # use first interval
  } else if (theta >= max(knots.cal)) {
    k <- which.max(knots.cal) # use last interval
  } else {
    k <- max(which(knots.cal <= theta))
  }
  
  t_k <- knots.cal[k]
  t_k1 <- knots.cal[k + 1]
  m_k <- values.cal[k]
  m_k1 <- values.cal[k + 1]
  
  # slope and intercept
  b_k <- (m_k1 - m_k) / (t_k1 - t_k)
  a_k <- m_k - b_k * t_k
  
  # return interpolated value mu(theta)
  res=a_k + b_k * theta
  return(res)
  }
  sapply(theta,k_map)
}


#  likelihood p(x | s, j, alpha)
likelihood_x <- function(x, s, alpha_j, alpha_j1) {
  val=x
  std=s
  lb=alpha_j
  ub=alpha_j1
  
  p_theta <- function(theta,val.=val,std.=std,lb.=lb,ub.=ub) {
    cat("theta =", theta, "\n")
    mu <- mu_theta(theta)
    cat("mu =", mu, "\n")
    p_x_given_mu_s = dnorm(val.,mean=mu,sd=std.)
    cat("norm dens = ",p_x_given_mu_s, "\n")
    p_theta_given_alpha <- dunif(x=theta, min=lb., max=ub.)
    cat("unif dens =", p_theta_given_alpha, "\n")
    #return(p_x_given_mu_s)
    return(p_x_given_mu_s * p_theta_given_alpha)
  }
  integrate(p_theta, lower = lb, upper = ub, val=val,std=std)$value
}

# likelihoods-- all observations for a given phase
#phase 1
phase1$likelihood <- mapply(likelihood_x, phase1$x, phase1$s, 
                            MoreArgs = list(alpha_j = 2000, alpha_j1 = 2800))
sum_likelihood_phase1 <- -sum(log(phase1$likelihood))

# phase 2: 2000 - 2650
#phase 3: 1900-2600
# phase 4: 1850 - 2450


# neg log likely cause mle() 
log_likelihood <- function(alpha_j, alpha_j1, data) {
  likelihoods <- mapply(likelihood_x, data$x, data$s, 
                        MoreArgs = list(alpha_j = alpha_j, alpha_j1 = alpha_j1))
  
  -sum(log(likelihoods))
}


# use mle() for phase
mle_phase1 <- mle(log_likelihood, 
                  start = list(alpha_j = 2000, alpha_j1 = 2900), 
                  data = phase1)

mle_phase2 <- mle(minuslogl = log_likelihood, 
                  start = list(alpha_j = 2000, alpha_j1 = 2650), 
                  data = phase2)

mle_phase3 <- mle(log_likelihood, 
                  start = list(alpha_j = 1900, alpha_j1 = 2600), 
                  data = phase3)

mle_phase4 <- mle(log_likelihood,
                  start = list(alpha_j = 1850, alpha_j1 = 2450), 
                  data = phase4)


