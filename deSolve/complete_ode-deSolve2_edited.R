# slphyx@SHIFT-ENTER
# the modified version of complete_ode for using with deSolve
# and foi from the Rcpp version
# parms  must be from read_data 
complete_ode.deSolve2 <- function(month, Y, parms) {
  
  step <- 1/parms$accuracy_step
  
  dY <- numeric(length(Y))
  
  loc_incidence <- numeric(parms$Nage)
  
  loc_incidence_U <- numeric(parms$Nage)
  
  loc_incidence_U0O <- numeric(parms$Nage)
  loc_incidence_U0R <- numeric(parms$Nage)
  loc_incidence_U1O <- numeric(parms$Nage)
  loc_incidence_U1R <- numeric(parms$Nage)
  loc_incidence_U2O <- numeric(parms$Nage)
  loc_incidence_U2R <- numeric(parms$Nage)
  loc_incidence_U3O <- numeric(parms$Nage)
  loc_incidence_U3R <- numeric(parms$Nage)
  
  loc_incidence_L <- numeric(parms$Nage)
  
  loc_incidence_L0O <- numeric(parms$Nage)
  loc_incidence_L0R <- numeric(parms$Nage)
  loc_incidence_L1O <- numeric(parms$Nage)
  loc_incidence_L1R <- numeric(parms$Nage)
  loc_incidence_L2O <- numeric(parms$Nage)
  loc_incidence_L2R <- numeric(parms$Nage)
  loc_incidence_L3O <- numeric(parms$Nage)
  loc_incidence_L3R <- numeric(parms$Nage)
  
  loc_vaccination <- numeric(parms$Nage)
   
  
  pURTI <- 1 - parms$pLRTI
  
  for (j in seq_along(parms$age_grps)) {
    
    loc_incidence[j] <- step*(Y[parms$Nage + j]*foi3(Y, month, j, parms) +
                                 Y[5*parms$Nage + j]*parms$sigmaR*foi3(Y, month, j, parms) +
                                 Y[10*parms$Nage + j]*(parms$sigmaL[j,1]*parms$pLRTI[j] + parms$sigmaU[j,1]*pURTI[j])*foi3(Y, month, j, parms) +
                                 Y[14*parms$Nage + j]*(min(parms$sigmaL[j,1], parms$sigmaR)*parms$pLRTI[j] + min(parms$sigmaU[j,1], parms$sigmaR)*pURTI[j])*foi3(Y, month, j, parms) +
                                 Y[18*parms$Nage + j]*(parms$sigmaL[j,2]*parms$pLRTI[j] + parms$sigmaU[j,2]*pURTI[j])*foi3(Y, month, j, parms) +
                                 Y[22*parms$Nage + j]*(min(parms$sigmaL[j,2],parms$sigmaR)*parms$pLRTI[j] + min(parms$sigmaU[j,2], parms$sigmaR)*pURTI[j])*foi3(Y, month, j, parms) +
                                 Y[26*parms$Nage + j]*(parms$sigmaL[j,3]*parms$pLRTI[j] + parms$sigmaU[j,3]*pURTI[j])*foi3(Y, month, j, parms) +
                                 Y[30*parms$Nage + j]*(min(parms$sigmaL[j,3],parms$sigmaR)*parms$pLRTI[j] + min(parms$sigmaU[j,3], parms$sigmaR)*pURTI[j])*foi3(Y, month, j, parms))
    
    loc_incidence_U[j] <- step*(Y[parms$Nage + j]*foi3(Y, month, j, parms)*pURTI[j] +
                                 Y[5*parms$Nage + j]*parms$sigmaR*foi3(Y, month, j, parms)*pURTI[j] +                             # infection
                                 Y[10*parms$Nage + j]*parms$sigmaU[j,1]*pURTI[j]*foi3(Y, month, j, parms) +
                                 Y[14*parms$Nage + j]*min(parms$sigmaU[j,1], parms$sigmaR)*pURTI[j]*foi3(Y, month, j, parms) +
                                 Y[18*parms$Nage + j]*parms$sigmaU[j,2]*pURTI[j]*foi3(Y, month, j, parms) +                    
                                 Y[22*parms$Nage + j]*min(parms$sigmaU[j,2], parms$sigmaR)*pURTI[j]*foi3(Y, month, j, parms) +
                                 Y[26*parms$Nage + j]*parms$sigmaU[j,3]*pURTI[j]*foi3(Y, month, j, parms) +                    
                                 Y[30*parms$Nage + j]*min(parms$sigmaU[j,3], parms$sigmaR)*pURTI[j]*foi3(Y, month, j, parms))
    
    loc_incidence_U0O[j] <- step*(Y[parms$Nage + j]*foi3(Y, month, j, parms)*pURTI[j])
    loc_incidence_U0R[j] <- step*(Y[5*parms$Nage + j]*parms$sigmaR*foi3(Y, month, j, parms)*pURTI[j])                             
    loc_incidence_U1O[j] <- step*(Y[10*parms$Nage + j]*parms$sigmaU[j,1]*pURTI[j]*foi3(Y, month, j, parms))
    loc_incidence_U1R[j] <- step*(Y[14*parms$Nage + j]*min(parms$sigmaU[j,1], parms$sigmaR)*pURTI[j]*foi3(Y, month, j, parms))
    loc_incidence_U2O[j] <- step*(Y[18*parms$Nage + j]*parms$sigmaU[j,2]*pURTI[j]*foi3(Y, month, j, parms))                    
    loc_incidence_U2R[j] <- step*(Y[22*parms$Nage + j]*min(parms$sigmaU[j,2], parms$sigmaR)*pURTI[j]*foi3(Y, month, j, parms))
    loc_incidence_U3O[j] <- step*(Y[26*parms$Nage + j]*parms$sigmaU[j,3]*pURTI[j]*foi3(Y, month, j, parms))                    
    loc_incidence_U3R[j] <- step*(Y[30*parms$Nage + j]*min(parms$sigmaU[j,3], parms$sigmaR)*pURTI[j]*foi3(Y, month, j, parms))
    
    loc_incidence_L[j] <- step*(Y[parms$Nage + j]*foi3(Y, month, j, parms)*parms$pLRTI[j] +
                                Y[5*parms$Nage + j]*parms$sigmaR*parms$pLRTI[j]*foi3(Y, month, j, parms) +
                                Y[10*parms$Nage + j]*parms$sigmaL[j,1]*parms$pLRTI[j]*foi3(Y, month, j, parms) +
                                Y[14*parms$Nage + j]*min(parms$sigmaL[j,1], parms$sigmaR)*parms$pLRTI[j]*foi3(Y, month, j, parms) +
                                Y[18*parms$Nage + j]*parms$sigmaL[j,2]*parms$pLRTI[j]*foi3(Y, month, j, parms) +
                                Y[22*parms$Nage + j]*min(parms$sigmaL[j,2], parms$sigmaR)*parms$pLRTI[j]*foi3(Y, month, j, parms) +
                                Y[26*parms$Nage + j]*parms$sigmaL[j,3]*parms$pLRTI[j]*foi3(Y, month, j, parms) +
                                Y[30*parms$Nage + j]*min(parms$sigmaL[j,3], parms$sigmaR)*parms$pLRTI[j]*foi3(Y, month, j, parms)) 
                       
    loc_incidence_L0O[j] <- step*(Y[parms$Nage + j]*foi3(Y, month, j, parms)*parms$pLRTI[j])
    loc_incidence_L0R[j] <- step*(Y[5*parms$Nage + j]*parms$sigmaR*parms$pLRTI[j]*foi3(Y, month, j, parms))
    loc_incidence_L1O[j] <- step*(Y[10*parms$Nage + j]*parms$sigmaL[j,1]*parms$pLRTI[j]*foi3(Y, month, j, parms))
    loc_incidence_L1R[j] <- step*(Y[14*parms$Nage + j]*min(parms$sigmaL[j,1], parms$sigmaR)*parms$pLRTI[j]*foi3(Y, month, j, parms))
    loc_incidence_L2O[j] <- step*(Y[18*parms$Nage + j]*parms$sigmaL[j,2]*parms$pLRTI[j]*foi3(Y, month, j, parms))
    loc_incidence_L2R[j] <- step*(Y[22*parms$Nage + j]*min(parms$sigmaL[j,2], parms$sigmaR)*parms$pLRTI[j]*foi3(Y, month, j, parms))
    loc_incidence_L3O[j] <- step*(Y[26*parms$Nage + j]*parms$sigmaL[j,3]*parms$pLRTI[j]*foi3(Y, month, j, parms))
    loc_incidence_L3R[j] <- step*(Y[30*parms$Nage + j]*min(parms$sigmaL[j,3], parms$sigmaR)*parms$pLRTI[j]*foi3(Y, month, j, parms)) 
    
    loc_vaccination[j] <- ifelse(j > 1, 
                                 step*(Y[parms$Nage + j]*parms$pVac[j, ConvMonth(month)] + Y[5*parms$Nage + j]*parms$pVac[j, ConvMonth(month)]),
                                 step*(parms$monthlyBirths*parms$pM[ConvMonth(month)] + Y[parms$Nage + j]*parms$pVac[j, ConvMonth(month)] + Y[5*parms$Nage + j]*parms$pVac[j, ConvMonth(month)]))  
   
    # MV0I0: 1
    dY[j] <- 
      (ifelse(j > 1, parms$kappa[j-1]*Y[j-1], parms$monthlyBirths*(1  - parms$pM[ConvMonth(month)])) - ifelse(j <parms$Nage ,parms$kappa[j]*Y[j],0) -  # ageing
              Y[j]*(parms$mortality[j] + parms$omegaM[j]))                                              # dead, waning immunity
    
    # SV0I0: 2
    dY[parms$Nage + j] <-
      (ifelse(j > 1, parms$kappa[j-1]*Y[parms$Nage + j-1], 0) - ifelse(j <parms$Nage ,parms$kappa[j]*Y[parms$Nage + j],0) +
              Y[j]*parms$omegaM[j] -                                                         # waning from unvaccinated maternal immunity
              Y[parms$Nage + j]*parms$mortality[j] -
              Y[parms$Nage + j]*foi3(Y, month, j, parms) -
              Y[parms$Nage + j]*parms$pVac[j, ConvMonth(month)])    # dead, infected, vaccinated
    
    # UV0I0: 3
    dY[2*parms$Nage + j] <-
      (ifelse(j > 1, parms$kappa[j-1]*Y[2*parms$Nage + j-1], 0) - ifelse(j <parms$Nage ,parms$kappa[j]*Y[2*parms$Nage + j],0) +
         Y[parms$Nage + j]*foi3(Y, month, j, parms)*pURTI[j] -                    
         Y[2*parms$Nage + j]*(parms$mortality[j] + parms$gamma0))                                    # dead, recovering 
    
    # LV0I0: 4
    dY[3*parms$Nage + j] <-
      (ifelse(j > 1, parms$kappa[j-1]*Y[3*parms$Nage + j-1], 0) - ifelse(j <parms$Nage ,parms$kappa[j]*Y[3*parms$Nage + j],0) +
         Y[parms$Nage + j]*foi3(Y, month, j, parms)*parms$pLRTI[j] -
         Y[3*parms$Nage + j]*(parms$mortality[j] + parms$gamma0))                                    # dead, recovering
    
    # RV0I0: 5
    dY[4*parms$Nage + j] <-
      (ifelse(j > 1, parms$kappa[j-1]*Y[4*parms$Nage + j-1], 0) - ifelse(j <parms$Nage ,parms$kappa[j]*Y[4*parms$Nage + j],0) +
         parms$gamma0*(Y[2*parms$Nage + j] + Y[3*parms$Nage + j]) -
         Y[4*parms$Nage + j]*(parms$mortality[j] + parms$rho))                                       # dead, recovered
    
    # SV0IR: 6
    dY[5*parms$Nage + j] <-
      (ifelse(j > 1, parms$kappa[j-1]*Y[5*parms$Nage + j-1], 0) - ifelse(j <parms$Nage ,parms$kappa[j]*Y[5*parms$Nage + j],0) +
         parms$rho*(Y[4*parms$Nage + j] + Y[8*parms$Nage + j]) +                                               # from recovery
         (Y[26*parms$Nage + j] + Y[30*parms$Nage + j])*parms$omega[j,3] -                                      # waning from post vaccine
         Y[5*parms$Nage + j]*parms$sigmaR*foi3(Y, month, j, parms) -
         Y[5*parms$Nage + j]*(parms$mortality[j] + parms$pVac[j, ConvMonth(month)]))      # risk adjusted infection, parms$mortality, vaccination
    
    # UV0IR: 7
    dY[6*parms$Nage + j] <-
      (ifelse(j > 1, parms$kappa[j-1]*Y[6*parms$Nage + j-1], 0) - ifelse(j <parms$Nage ,parms$kappa[j]*Y[6*parms$Nage + j],0) +
         (Y[27*parms$Nage + j] + Y[31*parms$Nage + j])*parms$omega[j,3] +                                      # waning from post vaccine
         Y[5*parms$Nage + j]*parms$sigmaR*foi3(Y, month, j, parms)*pURTI[j] -                             # infection
         Y[6*parms$Nage + j]*(parms$mortality[j] + parms$gammaR))
    
    # LV0IR: 8
    dY[7*parms$Nage + j] <-
      (ifelse(j > 1, parms$kappa[j-1]*Y[7*parms$Nage + j-1], 0) - ifelse(j <parms$Nage ,parms$kappa[j]*Y[7*parms$Nage + j],0) +
         (Y[28*parms$Nage + j] + Y[32*parms$Nage + j])*parms$omega[j,3] +                                      # waning from post vaccine
         Y[5*parms$Nage + j]*parms$sigmaR*parms$pLRTI[j]*foi3(Y, month, j, parms) -
         Y[7*parms$Nage + j]*(parms$mortality[j] + parms$gammaR))
    
    # RV0IR: 9
    dY[8*parms$Nage + j] <-
      (ifelse(j > 1, parms$kappa[j-1]*Y[8*parms$Nage + j-1], 0) - ifelse(j <parms$Nage ,parms$kappa[j]*Y[8*parms$Nage + j],0) +
         (Y[29*parms$Nage + j] + Y[33*parms$Nage + j])*parms$omega[j,3] +                                      # waning from post vaccine
         parms$gammaR*(Y[6*parms$Nage + j] + Y[7*parms$Nage + j]) -
         Y[8*parms$Nage + j]*(parms$mortality[j] + parms$rho))
    
    # -------------------------------------------------------------------------
    
    # MV1I0: 10
    dY[9*parms$Nage + j] <-
      (ifelse(j > 1, parms$kappa[j-1]*Y[9*parms$Nage + j-1], parms$monthlyBirths*parms$pM[ConvMonth(month)]) - ifelse(j <parms$Nage , parms$kappa[j]*Y[9*parms$Nage + j],0) -  # maternal vaccination
         Y[9*parms$Nage + j]*(parms$mortality[j] + parms$omegaM[j]))                                                # dead, waning immunity
    
    # SV1I0: 11
    dY[10*parms$Nage + j] <-
      (ifelse(j > 1, parms$kappa[j-1]*Y[10*parms$Nage + j-1], 0) - ifelse(j <parms$Nage , parms$kappa[j]*Y[10*parms$Nage + j],0) +
         Y[parms$Nage + j]*parms$pVac[j, ConvMonth(month)] +                                                                                                         # vaccinated
         Y[9*parms$Nage + j]*parms$omegaM[j] -                                                                                                            # waning from vaccinated maternal immunity
         Y[10*parms$Nage + j]*(parms$sigmaL[j,1]*parms$pLRTI[j] + parms$sigmaU[j,1]*pURTI[j])*foi3(Y, month, j, parms) -
         Y[10*parms$Nage + j]*(parms$mortality[j] + parms$omega[j,1]))  # dead, infected, vaccinated
    
    # UV1I0: 12
    dY[11*parms$Nage + j] <-
      (ifelse(j > 1, parms$kappa[j-1]*Y[11*parms$Nage + j-1], 0) - ifelse(j <parms$Nage , parms$kappa[j]*Y[11*parms$Nage + j],0) +
         Y[10*parms$Nage + j]*parms$sigmaU[j,1]*pURTI[j]*foi3(Y, month, j, parms) -              # adjusted infection
         Y[11*parms$Nage + j]*(parms$mortality[j] + parms$gamma0 + parms$omega[j,1]))                        # dead, recovered, waning vaccine
    
    # LV1I0: 13
    dY[12*parms$Nage + j] <-
      (ifelse(j > 1, parms$kappa[j-1]*Y[12*parms$Nage + j-1], 0) - ifelse(j <parms$Nage ,parms$kappa[j]*Y[12*parms$Nage + j],0) +
         Y[10*parms$Nage + j]*parms$sigmaL[j,1]*parms$pLRTI[j]*foi3(Y, month, j, parms) -              # adjusted infection
         Y[12*parms$Nage + j]*(parms$mortality[j] + parms$gamma0 + parms$omega[j,1]))                        # dead, recovered, waning vaccine
    
    # RV1I0: 14
    dY[13*parms$Nage + j] <-
      (ifelse(j > 1, parms$kappa[j-1]*Y[13*parms$Nage + j-1], 0) - ifelse(j <parms$Nage ,parms$kappa[j]*Y[13*parms$Nage + j],0) +
         parms$gamma0*(Y[11*parms$Nage + j] + Y[12*parms$Nage + j]) -
         Y[13*parms$Nage + j]*(parms$mortality[j] + parms$rho + parms$omega[j,1]))                           # dead, recovered, waning vaccine
    
    # SV1IR: 15
    dY[14*parms$Nage + j] <-
      (ifelse(j > 1, parms$kappa[j-1]*Y[14*parms$Nage + j-1], 0) - ifelse(j <parms$Nage ,parms$kappa[j]*Y[14*parms$Nage + j],0) +
         Y[5*parms$Nage + j]*parms$pVac[j, ConvMonth(month)] +                                                                                                               # vaccinated
         parms$rho*(Y[13*parms$Nage + j] + Y[17*parms$Nage + j]) -
         Y[14*parms$Nage + j]*(min(parms$sigmaL[j,1], parms$sigmaR)*parms$pLRTI[j] + min(parms$sigmaU[j,1], parms$sigmaR)*pURTI[j])*foi3(Y, month, j, parms) -  # recovered
         Y[14*parms$Nage + j]*(parms$mortality[j] + parms$omega[j,1]))
    
    # UV1IR: 16
    dY[15*parms$Nage + j] <-
      (ifelse(j > 1, parms$kappa[j-1]*Y[15*parms$Nage + j-1], 0) - ifelse(j <parms$Nage ,parms$kappa[j]*Y[15*parms$Nage + j],0) +
         Y[14*parms$Nage + j]*min(parms$sigmaU[j,1], parms$sigmaR)*pURTI[j]*foi3(Y, month, j, parms) -
         Y[15*parms$Nage + j]*(parms$mortality[j] + parms$gammaR + parms$omega[j,1]))
    
    # LV1IR: 17
    dY[16*parms$Nage + j] <-
      (ifelse(j > 1, parms$kappa[j-1]*Y[16*parms$Nage + j-1], 0) - ifelse(j <parms$Nage ,parms$kappa[j]*Y[16*parms$Nage + j],0) +
         Y[14*parms$Nage + j]*min(parms$sigmaL[j,1], parms$sigmaR)*parms$pLRTI[j]*foi3(Y, month, j, parms) -
         Y[16*parms$Nage + j]*(parms$mortality[j] + parms$gammaR + parms$omega[j,1]))
    
    # RV1IR: 18
    dY[17*parms$Nage + j] <-
      (ifelse(j > 1, parms$kappa[j-1]*Y[17*parms$Nage + j-1], 0) - ifelse(j <parms$Nage ,parms$kappa[j]*Y[17*parms$Nage + j],0) +
         parms$gammaR*(Y[15*parms$Nage + j] + Y[16*parms$Nage + j]) -
         Y[17*parms$Nage + j]*(parms$mortality[j] + parms$rho + parms$omega[j,1]))
    
    # -------------------------------------------------------------------------
    
    # SV2I0: 19
    dY[18*parms$Nage + j] <-
      (ifelse(j > 1, parms$kappa[j-1]*Y[18*parms$Nage + j-1], 0) - ifelse(j <parms$Nage ,parms$kappa[j]*Y[18*parms$Nage + j],0) +
         Y[10*parms$Nage + j]*parms$omega[j,1] -
         Y[18*parms$Nage + j]*parms$mortality[j] - 
         Y[18*parms$Nage + j]*(parms$sigmaL[j,2]*parms$pLRTI[j] + parms$sigmaU[j,2]*pURTI[j])*foi3(Y, month, j, parms) -
         Y[18*parms$Nage + j]*parms$omega[j,2])  # dead, infected, vaccinated, waning
    
    # UV2I0: 20
    dY[19*parms$Nage + j] <-
      (ifelse(j > 1, parms$kappa[j-1]*Y[19*parms$Nage + j-1], 0) - ifelse(j <parms$Nage ,parms$kappa[j]*Y[19*parms$Nage + j],0) +
         Y[11*parms$Nage + j]*parms$omega[j,1] +
         Y[18*parms$Nage + j]*parms$sigmaU[j,2]*pURTI[j]*foi3(Y, month, j, parms) -                    
         Y[19*parms$Nage + j]*(parms$mortality[j] + parms$gamma0 + parms$omega[j,2]))
    
    # LV2I0: 21
    dY[20*parms$Nage + j] <-
      (ifelse(j > 1, parms$kappa[j-1]*Y[20*parms$Nage + j-1], 0) - ifelse(j <parms$Nage ,parms$kappa[j]*Y[20*parms$Nage + j],0) +
         Y[12*parms$Nage + j]*parms$omega[j,1] +
         Y[18*parms$Nage + j]*parms$sigmaL[j,2]*parms$pLRTI[j]*foi3(Y, month, j, parms) -
         Y[20*parms$Nage + j]*(parms$mortality[j] + parms$gamma0 + parms$omega[j,2]))                        # dead, recovered, waning
    
    # RV2I0: 22
    dY[21*parms$Nage + j] <-
      (ifelse(j > 1, parms$kappa[j-1]*Y[21*parms$Nage + j-1], 0) - ifelse(j <parms$Nage ,parms$kappa[j]*Y[21*parms$Nage + j],0) +
         Y[13*parms$Nage + j]*parms$omega[j,1] +
         parms$gamma0*(Y[19*parms$Nage + j] + Y[20*parms$Nage + j]) -
         Y[21*parms$Nage + j]*(parms$mortality[j] + parms$rho + parms$omega[j,2]))                           # dead, recovered
    
    # SV2IR: 23
    dY[22*parms$Nage + j] <-
      (ifelse(j > 1, parms$kappa[j-1]*Y[22*parms$Nage + j-1], 0) - ifelse(j <parms$Nage ,parms$kappa[j]*Y[22*parms$Nage + j],0) +
         parms$rho*(Y[21*parms$Nage + j] + Y[25*parms$Nage + j]) +
         Y[14*parms$Nage + j]*parms$omega[j, 1] - 
         Y[22*parms$Nage + j]*parms$mortality[j] -
         Y[22*parms$Nage + j]*(min(parms$sigmaL[j,2],parms$sigmaR)*parms$pLRTI[j] + min(parms$sigmaU[j,2], parms$sigmaR)*pURTI[j])*foi3(Y, month, j, parms) -  # risk adjusted infection
         Y[22*parms$Nage + j]*parms$omega[j,2])  # risk adjusted infection
    
    # UV2IR: 24
    dY[23*parms$Nage + j] <-
      (ifelse(j > 1, parms$kappa[j-1]*Y[23*parms$Nage + j-1], 0) - ifelse(j <parms$Nage ,parms$kappa[j]*Y[23*parms$Nage + j],0) +
         Y[15*parms$Nage + j]*parms$omega[j,1] +
         Y[22*parms$Nage + j]*min(parms$sigmaU[j,2], parms$sigmaR)*pURTI[j]*foi3(Y, month, j, parms) -
         Y[23*parms$Nage + j]*(parms$mortality[j] + parms$gammaR + parms$omega[j,2]))
    
    # LV2IR: 25
    dY[24*parms$Nage + j] <-
      (ifelse(j > 1, parms$kappa[j-1]*Y[24*parms$Nage + j-1], 0) - ifelse(j <parms$Nage ,parms$kappa[j]*Y[24*parms$Nage + j],0) +
         Y[16*parms$Nage + j]*parms$omega[j,1] +
         Y[22*parms$Nage + j]*min(parms$sigmaL[j,2], parms$sigmaR)*parms$pLRTI[j]*foi3(Y, month, j, parms) -
         Y[24*parms$Nage + j]*(parms$mortality[j] + parms$gammaR + parms$omega[j,2]))
    
    # RV2IR: 26
    dY[25*parms$Nage + j] <-
      (ifelse(j > 1, parms$kappa[j-1]*Y[25*parms$Nage + j-1], 0) - ifelse(j <parms$Nage ,parms$kappa[j]*Y[25*parms$Nage + j],0) +
         Y[17*parms$Nage + j]*parms$omega[j,1] +
         parms$gammaR*(Y[23*parms$Nage + j] + Y[24*parms$Nage + j]) -
         Y[25*parms$Nage + j]*(parms$mortality[j] + parms$rho + parms$omega[j,2]))
    
    # -------------------------------------------------------------------------
    
    # SV3I0: 27
    dY[26*parms$Nage + j] <-
      (ifelse(j > 1, parms$kappa[j-1]*Y[26*parms$Nage + j-1], 0) - ifelse(j <parms$Nage ,parms$kappa[j]*Y[26*parms$Nage + j],0) +
         Y[18*parms$Nage + j]*parms$omega[j,2] -
         Y[26*parms$Nage + j]*parms$mortality[j] -
         Y[26*parms$Nage + j]*(parms$sigmaL[j,3]*parms$pLRTI[j] + parms$sigmaU[j,3]*pURTI[j])*foi3(Y, month, j, parms) -   # dead, infected, vaccinated
         Y[26*parms$Nage + j]*parms$omega[j,3])   # dead, infected, vaccinated
    
    # UV3I0: 28
    dY[27*parms$Nage + j] <-
      (ifelse(j > 1, parms$kappa[j-1]*Y[27*parms$Nage + j-1], 0) - ifelse(j <parms$Nage ,parms$kappa[j]*Y[27*parms$Nage + j],0) +
         Y[19*parms$Nage + j]*parms$omega[j,2] +
         Y[26*parms$Nage + j]*parms$sigmaU[j,3]*pURTI[j]*foi3(Y, month, j, parms) -                    
         Y[27*parms$Nage + j]*(parms$mortality[j] + parms$gamma0 + parms$omega[j,3]))
    
    # LV3I0: 29
    dY[28*parms$Nage + j] <-
      (ifelse(j > 1, parms$kappa[j-1]*Y[28*parms$Nage + j-1], 0) - ifelse(j <parms$Nage ,parms$kappa[j]*Y[28*parms$Nage + j],0) +
         Y[20*parms$Nage + j]*parms$omega[j,2] +
         Y[26*parms$Nage + j]*parms$sigmaL[j,3]*parms$pLRTI[j]*foi3(Y, month, j, parms) -
         Y[28*parms$Nage + j]*(parms$mortality[j] + parms$gamma0 + parms$omega[j,3]))                        # dead, recovered 
    
    # RV3I0: 30
    dY[29*parms$Nage + j] <-
      (ifelse(j > 1, parms$kappa[j-1]*Y[29*parms$Nage + j-1], 0) - ifelse(j <parms$Nage ,parms$kappa[j]*Y[29*parms$Nage + j],0) +
         Y[21*parms$Nage + j]*parms$omega[j,2] +
         parms$gamma0*(Y[27*parms$Nage + j] + Y[28*parms$Nage + j]) -
         Y[29*parms$Nage + j]*(parms$mortality[j] + parms$rho + parms$omega[j,3]))                          # dead, recovered
    
    # SV3IR: 31
    dY[30*parms$Nage + j] <-
      (ifelse(j > 1, parms$kappa[j-1]*Y[30*parms$Nage + j-1], 0) - ifelse(j <parms$Nage ,parms$kappa[j]*Y[30*parms$Nage + j],0) +
         parms$rho*(Y[29*parms$Nage + j] + Y[33*parms$Nage + j]) +
         Y[22*parms$Nage + j]*parms$omega[j, 2] - 
         Y[30*parms$Nage + j]*parms$mortality[j] -
         Y[30*parms$Nage + j]*(min(parms$sigmaL[j,3],parms$sigmaR)*parms$pLRTI[j] + min(parms$sigmaU[j,3], parms$sigmaR)*pURTI[j])*foi3(Y, month, j, parms) -
         Y[30*parms$Nage + j]*parms$omega[j,3])  # risk adjusted infection
    
    # UV3IR: 32
    dY[31*parms$Nage + j] <-
      (ifelse(j > 1, parms$kappa[j-1]*Y[31*parms$Nage + j-1], 0) - ifelse(j <parms$Nage ,parms$kappa[j]*Y[31*parms$Nage + j],0) +
         Y[23*parms$Nage + j]*parms$omega[j,2] +
         Y[30*parms$Nage + j]*min(parms$sigmaU[j,3], parms$sigmaR)*pURTI[j]*foi3(Y, month, j, parms) -
         Y[31*parms$Nage + j]*(parms$mortality[j] + parms$gammaR + parms$omega[j,3]))
    
    # LV3IR: 33
    dY[32*parms$Nage + j] <-
      (ifelse(j > 1, parms$kappa[j-1]*Y[32*parms$Nage + j-1], 0) - ifelse(j <parms$Nage ,parms$kappa[j]*Y[32*parms$Nage + j],0) +
         Y[24*parms$Nage + j]*parms$omega[j,2] +
         Y[30*parms$Nage + j]*min(parms$sigmaL[j,3], parms$sigmaR)*parms$pLRTI[j]*foi3(Y, month, j, parms) -
         Y[32*parms$Nage + j]*(parms$mortality[j] + parms$gammaR + parms$omega[j,3]))
    
    # RV3IR: 34
    dY[33*parms$Nage + j] <-
      (ifelse(j > 1, parms$kappa[j-1]*Y[33*parms$Nage + j-1], 0) - ifelse(j <parms$Nage ,parms$kappa[j]*Y[33*parms$Nage + j],0) +
         Y[25*parms$Nage + j]*parms$omega[j,2] +
         parms$gammaR*(Y[31*parms$Nage + j] + Y[32*parms$Nage + j]) -
         Y[33*parms$Nage + j]*(parms$mortality[j] + parms$rho + parms$omega[j,3]))
    
  }  #close 'for' loop defining the transitions between states by age group.
  
  # collect all the incidence data
  incidence <<- incidence + loc_incidence 
  incidence_U <<- incidence_U + loc_incidence_U 
        
  incidence_U0O <<- incidence_U0O + loc_incidence_U0O 
  incidence_U0R <<- incidence_U0R + loc_incidence_U0R 
  incidence_U1O <<- incidence_U1O + loc_incidence_U1O 
  incidence_U1R <<- incidence_U1R + loc_incidence_U1R 
  incidence_U2O <<- incidence_U2O + loc_incidence_U2O 
  incidence_U2R <<- incidence_U2R + loc_incidence_U2R 
  incidence_U3O <<- incidence_U3O + loc_incidence_U3O 
  incidence_U3R <<- incidence_U3R + loc_incidence_U3R 
        
  incidence_L <<- incidence_L + loc_incidence_L 
        
  incidence_L0O <<- incidence_L0O + loc_incidence_L0O 
  incidence_L0R <<- incidence_L0R + loc_incidence_L0R 
  incidence_L1O <<- incidence_L1O + loc_incidence_L1O 
  incidence_L1R <<- incidence_L1R + loc_incidence_L1R 
  incidence_L2O <<- incidence_L2O + loc_incidence_L2O 
  incidence_L2R <<- incidence_L2R + loc_incidence_L2R 
  incidence_L3O <<- incidence_L3O + loc_incidence_L3O 
  incidence_L3R <<- incidence_L3R + loc_incidence_L3R 
        
  vaccination <<- vaccination + loc_vaccination 

  # export monthly incidence
  if(month >= 1){

    # gb.indx <<- gb.indx + 1
    if(trunc(month) > gb.indx){
      # cat('\n month : ',month,' gb.indx ',gb.indx ,' incidence_monthly ', incidence[1:2])
      ClearIncidenceVars(parms$Nage)
      gb.indx <<- gb.indx + 1
    }
    
    incidence_monthly[, trunc(month)] <<- incidence   
    
    incidence_monthly_U[,trunc(month)] <<-incidence_U
      
    incidence_monthly_U0O[,trunc(month)] <<-incidence_U0O
    incidence_monthly_U0R[,trunc(month)] <<-incidence_U0R
    incidence_monthly_U1O[,trunc(month)] <<-incidence_U1O
    incidence_monthly_U1R[,trunc(month)] <<-incidence_U1R
    incidence_monthly_U2O[,trunc(month)] <<-incidence_U2O
    incidence_monthly_U2R[,trunc(month)] <<-incidence_U2R
    incidence_monthly_U3O[,trunc(month)] <<-incidence_U3O
    incidence_monthly_U3R[,trunc(month)] <<-incidence_U3R
      
    incidence_monthly_L[,trunc(month)] <<-incidence_L
      
    incidence_monthly_L0O[,trunc(month)] <<-incidence_L0O
    incidence_monthly_L0R[,trunc(month)] <<-incidence_L0R
    incidence_monthly_L1O[,trunc(month)] <<-incidence_L1O
    incidence_monthly_L1R[,trunc(month)] <<-incidence_L1R
    incidence_monthly_L2O[,trunc(month)] <<-incidence_L2O
    incidence_monthly_L2R[,trunc(month)] <<-incidence_L2R
    incidence_monthly_L3O[,trunc(month)] <<-incidence_L3O
    incidence_monthly_L3R[,trunc(month)] <<-incidence_L3R
      
    vaccination_monthly[,trunc(month)] <<-vaccination
    
  }

  # return as a list for using with deSolve
  list(c(dY))
} 


