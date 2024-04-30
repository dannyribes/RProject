
## combine two system functions into one system function  
ode_system.deSolve <- function(month, Y, parms){
  if(month <= 12*parms$initial_period) 
    return(complete_ode_initial.deSolve(month, Y, parms))
  else 
    return(complete_ode.deSolve(month, Y, parms))
    
}


# for using with incidence module from the original code
ode_system.deSolve2 <- function(month, Y, parms) {
  
  if (exists("pb.DTM", envir = .GlobalEnv))
    setTxtProgressBar(pb.DTM, month / max(times))
  
  if (month <= 12 * parms$initial_period)
    return(complete_ode_initial.deSolve2(month, Y, parms))
  else
    return(complete_ode.deSolve2(month, Y, parms))

}


# run the model using Euler method from deSolve 
# complete_ode_initial.list is from complete_ode_initial but step = 1 and dY is returned as a list
#
# args - c(WORKINGPATH, paramsset, vac_type)
#
# this deSolve2 uses the incidence model from the original code
runDTMX.deSolve2 <- function(args, vac_type, init.pop = "data/startpop1.csv",
  APPEND = FALSE, method = "euler",  ProgressBar = TRUE, Write.Output = FALSE) {

  # load the input parameters
  params <- read_data2(args[2],vac_type)
  
  # load population structure data
  start_pop1 <- read.csv(init.pop)
  
  time_horizon <- params$time_horizon  # the unit is in years
  state_names <- params$state_names
  age_grps <- params$age_grps
  Nage <- length(age_grps)
  Nstate <- length(state_names)
  
  # intial population structure
  state_ini1 <- matrix(data = 0, nrow = Nage*34, ncol = 1)
  state_ini1[,1] <- start_pop1[, "start_pop"]
  # print information
  printInfo(params)
  
  # Progress Bar
  if(ProgressBar==TRUE){
    assign("pb.DTM",NULL,envir = .GlobalEnv)
    #ss for showing the progress bar
    pb.DTM <<- txtProgressBar(style = 3)
    #ss
    
  }else{
    rm("pb.DTM", envir = .GlobalEnv) 
  }
  
  # vars for keeping incidence values
  # these vars will be exported to .GlobalEnv 
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
        
  incidence_monthly <<- matrix(NA, nrow = Nage, ncol = 12*time_horizon)
  incidence_monthly_U <<- matrix(NA, nrow = Nage, ncol = 12*time_horizon)
  
  incidence_monthly_U0O <<- matrix(NA, nrow = Nage, ncol = 12*time_horizon)
  incidence_monthly_U0R <<- matrix(NA, nrow = Nage, ncol = 12*time_horizon)
  incidence_monthly_U1O <<- matrix(NA, nrow = Nage, ncol = 12*time_horizon)
  incidence_monthly_U1R <<- matrix(NA, nrow = Nage, ncol = 12*time_horizon)
  incidence_monthly_U2O <<- matrix(NA, nrow = Nage, ncol = 12*time_horizon)
  incidence_monthly_U2R <<- matrix(NA, nrow = Nage, ncol = 12*time_horizon)
  incidence_monthly_U3O <<- matrix(NA, nrow = Nage, ncol = 12*time_horizon)
  incidence_monthly_U3R <<- matrix(NA, nrow = Nage, ncol = 12*time_horizon)
  
  incidence_monthly_L <<- matrix(NA, nrow = Nage, ncol = 12*time_horizon)
  
  incidence_monthly_L0O <<- matrix(NA, nrow = Nage, ncol = 12*time_horizon)
  incidence_monthly_L0R <<- matrix(NA, nrow = Nage, ncol = 12*time_horizon)
  incidence_monthly_L1O <<- matrix(NA, nrow = Nage, ncol = 12*time_horizon)
  incidence_monthly_L1R <<- matrix(NA, nrow = Nage, ncol = 12*time_horizon)
  incidence_monthly_L2O <<- matrix(NA, nrow = Nage, ncol = 12*time_horizon)
  incidence_monthly_L2R <<- matrix(NA, nrow = Nage, ncol = 12*time_horizon)
  incidence_monthly_L3O <<- matrix(NA, nrow = Nage, ncol = 12*time_horizon)
  incidence_monthly_L3R <<- matrix(NA, nrow = Nage, ncol = 12*time_horizon)
  
  vaccination_monthly <<- matrix(NA, nrow = Nage, ncol = 12*time_horizon)
 

  # time in months
  times <<- seq(from = 1, to = (12 * time_horizon) + 1, by = 1)
  
  # global index for counting the months
  gb.indx <<- 1
  
  # solve the system using the ode solver from the deSolve package
  out <- ode(y = state_ini1, times = times,func = ode_system.deSolve2, 
    parms = params, method = method, hini = 1/params$accuracy_step)
  
  if(ProgressBar==TRUE){
    rm("pb.DTM",envir = .GlobalEnv)
  }
  
  # write the monthly outputs
  if(Write.Output == TRUE){

    ICON.traj(out, file = paste0('results/Monthly_trajectories', '_', vac_type,'_', args[3], '.csv'))
    
    # write.table(traj_monthly, paste0('results/Monthly_trajectories','_',vac_type,'_',args[3],'.csv'), append = APPEND, sep = ",", col.names = FALSE, row.names = FALSE, quote = FALSE)
    
    # write.table(incidence_monthly,paste0('results/Incidence_count','_',vac_type,'_',args[3],'.csv'), append = APPEND, sep = ",", col.names = FALSE, row.names = FALSE, quote = FALSE)
    
    # write.table(incidence_monthly_U,paste0('results/Incidence_count_U','_',vac_type,'_',args[3],'.csv'), append = APPEND, sep = ",", col.names = FALSE, row.names = FALSE, quote = FALSE)
    # write.table(incidence_monthly_L,paste0('results/Incidence_count_L','_',vac_type,'_',args[3],'.csv'), append = APPEND, sep = ",", col.names = FALSE, row.names = FALSE, quote = FALSE)
    
    
    write.table(incidence_monthly_U0O,paste0('results/Incidence_count_UVI','_',vac_type,'_',args[3],'.csv'), append = APPEND, sep = ",", col.names = FALSE, row.names = FALSE, quote = FALSE)
    write.table(incidence_monthly_U0R,paste0('results/Incidence_count_UVI','_',vac_type,'_',args[3],'.csv'), append = TRUE, sep = ",", col.names = FALSE, row.names = FALSE, quote = FALSE)
    write.table(incidence_monthly_U1O,paste0('results/Incidence_count_UVI','_',vac_type,'_',args[3],'.csv'), append = TRUE, sep = ",", col.names = FALSE, row.names = FALSE, quote = FALSE)
    write.table(incidence_monthly_U1R,paste0('results/Incidence_count_UVI','_',vac_type,'_',args[3],'.csv'), append = TRUE, sep = ",", col.names = FALSE, row.names = FALSE, quote = FALSE)
    write.table(incidence_monthly_U2O,paste0('results/Incidence_count_UVI','_',vac_type,'_',args[3],'.csv'), append = TRUE, sep = ",", col.names = FALSE, row.names = FALSE, quote = FALSE)
    write.table(incidence_monthly_U2R,paste0('results/Incidence_count_UVI','_',vac_type,'_',args[3],'.csv'), append = TRUE, sep = ",", col.names = FALSE, row.names = FALSE, quote = FALSE)
    write.table(incidence_monthly_U3O,paste0('results/Incidence_count_UVI','_',vac_type,'_',args[3],'.csv'), append = TRUE, sep = ",", col.names = FALSE, row.names = FALSE, quote = FALSE)
    write.table(incidence_monthly_U3R,paste0('results/Incidence_count_UVI','_',vac_type,'_',args[3],'.csv'), append = TRUE, sep = ",", col.names = FALSE, row.names = FALSE, quote = FALSE)
    
    write.table(incidence_monthly_L0O,paste0('results/Incidence_count_LVI','_',vac_type,'_',args[3],'.csv'), append = APPEND, sep = ",", col.names = FALSE, row.names = FALSE, quote = FALSE)
    write.table(incidence_monthly_L0R,paste0('results/Incidence_count_LVI','_',vac_type,'_',args[3],'.csv'), append = TRUE, sep = ",", col.names = FALSE, row.names = FALSE, quote = FALSE)
    write.table(incidence_monthly_L1O,paste0('results/Incidence_count_LVI','_',vac_type,'_',args[3],'.csv'), append = TRUE, sep = ",", col.names = FALSE, row.names = FALSE, quote = FALSE)
    write.table(incidence_monthly_L1R,paste0('results/Incidence_count_LVI','_',vac_type,'_',args[3],'.csv'), append = TRUE, sep = ",", col.names = FALSE, row.names = FALSE, quote = FALSE)
    write.table(incidence_monthly_L2O,paste0('results/Incidence_count_LVI','_',vac_type,'_',args[3],'.csv'), append = TRUE, sep = ",", col.names = FALSE, row.names = FALSE, quote = FALSE)
    write.table(incidence_monthly_L2R,paste0('results/Incidence_count_LVI','_',vac_type,'_',args[3],'.csv'), append = TRUE, sep = ",", col.names = FALSE, row.names = FALSE, quote = FALSE)
    write.table(incidence_monthly_L3O,paste0('results/Incidence_count_LVI','_',vac_type,'_',args[3],'.csv'), append = TRUE, sep = ",", col.names = FALSE, row.names = FALSE, quote = FALSE)
    write.table(incidence_monthly_L3R,paste0('results/Incidence_count_LVI','_',vac_type,'_',args[3],'.csv'), append = TRUE, sep = ",", col.names = FALSE, row.names = FALSE, quote = FALSE)
    
    write.table(vaccination_monthly,paste0('results/Vaccine_count','_',vac_type,'_',args[3],'.csv'), append = APPEND, sep = ",", col.names = FALSE, row.names = FALSE, quote = FALSE)

  }

  # return monthly incidence
  output <- incidence_monthly
  return(output)
}

