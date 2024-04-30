foi <- function(Y, t, j, params,
                inf_state_names = c("UV0I0", "UV0IR", "LV0I0", "LV0IR","UV1I0", "UV1IR", "LV1I0", "LV1IR","UV2I0", "UV2IR", "LV2I0", "LV2IR","UV3I0", "UV3IR", "LV3I0", "LV3IR")) {
  
  
  # N[i]: total population in age i
  # phi: peak infection time
  # xi: infection amplitude
  
  n_age <- length(params$age_grps)
  pop <- matrix(Y, nrow = n_age, byrow = FALSE)
  colnames(pop) <- params$state_names
  rownames(pop) <- params$age_grps
  
  inf <- pop[, inf_state_names]
  N <- rowSums(pop)
  
  #sh=0.17
  #c_season <- params$xi*((sh/2)*(1 + cos(pi/6*(t - params$phi)))+(1-sh))
  
  shift <-6
  
  if((t+shift)>12) {
    c_season <- params$xi*cos(0.018181818*pi*(t + shift - 12 - 1))
  } else {
    c_season <- params$xi*cos(0.018181818*pi*(t + shift - 1))
  }
  
  valpha <- c(params$alpha)
  c_infect <- c(inf %*% valpha)
  res <- 0
  
  for (i in seq_len(n_age)) {
    
    res <- res + ifelse(N[i] != 0, params$contact_adj[j]*params$contact[j, i]/N[i], 0) * c_season * c_infect[i]
  }
  return(res)
}

foi3 <- function(Y, t, j, params) {
  # inf_state_names <- c("UV0I0", "UV0IR", "LV0I0", "LV0IR","UV1I0", "UV1IR", "LV1I0", "LV1IR","UV2I0", "UV2IR",
  #                     "LV2I0", "LV2IR","UV3I0", "UV3IR", "LV3I0", "LV3IR")
  
  # position inf_states <- c(3, 7, 4, 8, 12, 16, 13, 17, 20, 24, 21, 25, 28, 32, 29, 33) - 1
  
  contact <- params$contact
  contact_adj <- params$contact_adj
  age_grps <- params$age_grps
  inf_state_pos <- c(3, 7, 4, 8, 12, 16, 13, 17, 20, 24, 21, 25, 28, 32, 29, 33)
  n_age <- length(age_grps)
  
  pop <- matrix(0, n_age, 34)
  pop[] <- Y
  
  # N[i]: total population in age i
  # phi: peak infection time
  # xi: infection amplitude
  inf <- pop[, inf_state_pos]
  N <- rowSums(pop)
  
  c_season <- NULL
  xi <- params$xi
  phi <- params$phi
  ns <- params$ns
  
  # convert t into month in a year (1-12)
  tt <- t %% 12
  if (tt == 0) tt <- 12
  
  c_season <- (xi * (2^(ns - 1))) * (cos(2 * pi * (tt - phi) / 12) + 1)^ns + 1 - xi
  

  valpha <- params$alpha
  c_infect <- as.vector(inf %*% valpha)
  
  res <- 0
  
  for (i in seq_len(n_age)) {
    res <- res + ifelse(N[i] != 0, contact_adj[j] * contact[j, i] / N[i], 0) * c_season * c_infect[i]
  }
  
  return(res)
}
