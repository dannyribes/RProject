# read_data (from the original code) is modified to read 
# the contact matrix as the matrix type for using with Rcpp foi function
# i - set of parameter files
read_data <- function(i, vac_type) {
  
  if (i == 0) 
  {
    df1 <- read.csv("data/params.csv", header = TRUE)
    i <- 1
  }
  else
  {
    df1 <- read.csv(paste0("data/params_",args[3],".csv"), header = TRUE)
  }
  
  i <- as.numeric(i)
  
  start <- (i-1)*34 + 1
  end <- i*34
  prms <- lapply(df1[start:end,1:45], function(col)col[!is.na(col)])
  
  age_grps <- prms[[1]]  
  
  state_names <- prms[[2]]
  
  Nage <- length(age_grps)
  Nstate <- length(state_names)
  
  initial_period <- ifelse((vac_type==1),prms[[5]],prms[[3]])
  
  accuracy_step <- prms[[4]]
  
  time_horizon <- prms[[5]]
  startYear <-prms[[6]]
  startMonth <- prms[[7]]
  monthlyBirths <- prms[[8]]
  xi <- prms[[9]]
  phi <- prms[[10]]
  q <- prms[[11]]
  gamma0 <- prms[[12]]
  gammaR <- prms[[13]]
  omegaVac <- prms[[14]] 
  lambda <- prms[[15]] 
  sigmaR <- prms[[16]]
  
  
  pM <- prms[[17]] #TB: testing without conversion for maternal vaccination, I think ODE treats this as a proportion, not a rate: -log(1-prms[[17]])

  pVac<- matrix(c(-log(1-prms[[33]]), -log(1-prms[[34]]), -log(1-prms[[35]]),-log(1-prms[[36]]),-log(1-prms[[37]]),-log(1-prms[[38]]),
                  -log(1-prms[[39]]),-log(1-prms[[40]]),-log(1-prms[[41]]),-log(1-prms[[42]]),-log(1-prms[[43]]),-log(1-prms[[44]])), ncol = 12)


  
  
  
  rho <- prms[[18]]
  omegaM <- rep(prms[[19]], Nage) #2 converted from time in state estimate
  alpha <-matrix(prms[[20]], nrow = 2,byrow = TRUE)
  mortality<-prms[[21]]
  
  contact <- as.matrix(read.csv("data/contactMatrix.csv", header = FALSE))  # modified to matrix type
  
  sigmaL <- matrix(c(prms[[22]], prms[[23]], prms[[24]]), ncol = 3)
  sigmaU <- matrix(c(prms[[25]], prms[[26]], prms[[27]]), ncol = 3)
  omega <- matrix(c(prms[[28]], prms[[29]], prms[[30]]), ncol = 3)
  
  kaging<-prms[[31]]
  pLRTI<-prms[[32]]
  
  contact_adj<-prms[[45]]
  
  list(age_grps = age_grps,
       Nage = Nage,   #MZG
       monthlyBirths = monthlyBirths,
       startMonth = startMonth, 
       state_names = state_names,
       initial_period = initial_period,
       accuracy_step = accuracy_step,
       time_horizon = time_horizon,
       kappa = kaging,
       omegaM = omegaM,
       xi = xi,
       phi = phi,
       pLRTI = pLRTI,
       pM = pM,
       gamma0 = gamma0,
       rho = rho,
       sigmaR = sigmaR,
       gammaR = gammaR,
       alpha = alpha,
       mortality = mortality,
       contact = contact,
       sigmaL = sigmaL,
       sigmaU = sigmaU,
       omega = omega,
       pVac = pVac,
       contact_adj = contact_adj)
}


# converting deSolve output to ICON traj_output
#  => (compartments,time)
ICON.traj <- function(odeout, file = NULL) {
  odeout.size <- dim(odeout)
  
  # drop the time column and the first row and then do the transpose
  tmp <- odeout[, 2:odeout.size[2]]
  tmp <- tmp[2:odeout.size[1], ]
  colnames(tmp) <- NULL
  tmp <- t(tmp)
  
  if (is.null(file)) {
    return(tmp)
  } else{
    write.table(
      tmp,
      file = file,
      append = T, # added this line for appending trajectory
      sep = ",",
      row.names = F,
      col.names = F
    )
    message(paste("\n", file, "has been written.\n"))
  }
  
}


# convert t into month in a year (1-12)
# for using with parms$pVac i.e. parms$pVac[j, ConvMonth(month)])  
# in complete_ode function
ConvMonth <- function(t){
  tt <- (t %% 12)
  if (tt ==0) tt <- 12
  
  return(ceiling(tt))
}


read_data2 <- function(i, vac_type, parameter.file = NULL, contactmatrix.file = NULL) {
  
  if (i  == 0)
  {
    if(is.null(parameter.file)==TRUE){
      df1 <- read.csv("data/params.csv", header = TRUE)
    }else{
      df1 <- read.csv(file = parameter.file, header = TRUE)
    }
    i <- 1
  }
  else
  {
    if(is.null(parameter.file)==TRUE){
      df1 <- read.csv(paste0("data/params_",args[3],".csv"), header = TRUE)
    }else{
      df1 <- read.csv(file = parameter.file, header = TRUE)
    }
  }
  
  i <- as.numeric(i)
  
  start <- (i-1)*34 + 1
  end <- i*34
  prms <- lapply(df1[start:end,1:46], function(col)col[!is.na(col)])
  
  age_grps <- prms[[1]]
  
  state_names <- prms[[2]]
  
  Nage <- length(age_grps)
  Nstate <- length(state_names)
  
  # initial_period <- ifelse((vac_type == 1), prms[[5]], prms[[3]])
  initial_period <- prms[[3]]
  
  accuracy_step <- prms[[4]]
  
  time_horizon <- prms[[5]]
  startYear <-prms[[6]]
  startMonth <- prms[[7]]
  monthlyBirths <- prms[[8]]
  xi <- prms[[9]]
  phi <- prms[[10]]
  q <- prms[[11]]
  gamma0 <- prms[[12]]
  gammaR <- prms[[13]]
  omegaVac <- prms[[14]]
  lambda <- prms[[15]]
  sigmaR <- prms[[16]]
  
  # pM <- prms[[17]]
  # pVac<- matrix(c(-log(1-prms[[33]]), -log(1-prms[[34]]), -log(1-prms[[35]]), -log(1-prms[[36]]), -log(1-prms[[37]]), -log(1-prms[[38]]),
  #                 -log(1-prms[[39]]), -log(1-prms[[40]]), -log(1-prms[[41]]), -log(1-prms[[42]]), -log(1-prms[[43]]), -log(1-prms[[44]])), ncol = 12)
  
  # following conditional serves to provide the world without vaccination (vac_type == 1) with fixed 0 prob of vaccination in all age groups.
  if (vac_type == 1) {
    pM <- rep(0, 12)
    pVac <- matrix(0, nrow = 26, ncol = 12)
  } else {
    pM <- prms[[17]]
    pVac<- matrix(c(-log(1-prms[[33]]), -log(1-prms[[34]]), -log(1-prms[[35]]),-log(1-prms[[36]]),-log(1-prms[[37]]),-log(1-prms[[38]]),
                    -log(1-prms[[39]]), -log(1-prms[[40]]), -log(1-prms[[41]]),-log(1-prms[[42]]),-log(1-prms[[43]]),-log(1-prms[[44]])), ncol = 12)
  }
  
  rho <- prms[[18]]
  omegaM <- rep(prms[[19]], Nage) #2 converted from time in state estimate
  alpha <-matrix(prms[[20]], nrow = 2,byrow = TRUE)
  mortality<-prms[[21]]
  
  if(is.null(contactmatrix.file)==TRUE){
    contact <- as.matrix(read.csv("data/contactMatrix.csv", header = FALSE))
  }else{
    contact <- as.matrix(read.csv(contactmatrix.file, header = FALSE))
  }
  
  sigmaL <- matrix(c(prms[[22]], prms[[23]], prms[[24]]), ncol = 3)
  sigmaU <- matrix(c(prms[[25]], prms[[26]], prms[[27]]), ncol = 3)
  omega <- matrix(c(prms[[28]], prms[[29]], prms[[30]]), ncol = 3)
  
  kaging<-prms[[31]]
  pLRTI<-prms[[32]]
  
  contact_adj<-prms[[45]]
  
  ns <- prms[[46]]
  
  list(age_grps = age_grps,
       Nage = Nage,   #MZG
       monthlyBirths = monthlyBirths,
       startMonth = startMonth,
       state_names = state_names,
       initial_period = initial_period,
       accuracy_step = accuracy_step,
       time_horizon = time_horizon,
       kappa = kaging,
       omegaM = omegaM,
       xi = xi,
       phi = phi,
       pLRTI = pLRTI,
       pM = pM,
       gamma0 = gamma0,
       rho = rho,
       sigmaR = sigmaR,
       gammaR = gammaR,
       alpha = alpha,
       mortality = mortality,
       contact = contact,
       sigmaL = sigmaL,
       sigmaU = sigmaU,
       omega = omega,
       pVac = pVac,
       contact_adj = contact_adj,
       ns = ns)
}


# clear all incidence data
# for reset the values in each month
ClearIncidenceVars <- function(Nage){
  
  # vars for keeping incidence values
  incidence <<- numeric(Nage)
  incidence_U <<- numeric(Nage)
  
  incidence_U0O <<- numeric(Nage)
  incidence_U0R <<- numeric(Nage)
  incidence_U1O <<- numeric(Nage)
  incidence_U1R <<- numeric(Nage)
  incidence_U2O <<- numeric(Nage)
  incidence_U2R <<- numeric(Nage)
  incidence_U3O <<- numeric(Nage)
  incidence_U3R <<- numeric(Nage)
  
  incidence_L <<- numeric(Nage)
  
  incidence_L0O <<- numeric(Nage)
  incidence_L0R <<- numeric(Nage)
  incidence_L1O <<- numeric(Nage)
  incidence_L1R <<- numeric(Nage)
  incidence_L2O <<- numeric(Nage)
  incidence_L2R <<- numeric(Nage)
  incidence_L3O <<- numeric(Nage)
  incidence_L3R <<- numeric(Nage)
  
  vaccination <<- numeric(Nage)
  
}


# print information
# phi: peak infection time
# xi: infection amplitude
printInfo <- function(params){
  message(paste("\ninitial period (years) :", params$initial_period))
  message(paste("time horizon (years) :", params$time_horizon))
  message(paste("accuracy step :", params$accuracy_step))
  message(paste("infection amplitude (xi) :", params$xi))
  message(paste("peak infection time (phi) :", params$phi))
  message(paste("recovery rate from temporary immunity to susceptible (rho) :", params$rho))
  message(paste("waning of maternal immunity (omegaM) :", params$omegaM[1]))
  message(paste("pointiness parameter (ns) :", params$ns))
}

