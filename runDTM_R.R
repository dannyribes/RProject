
###################################################################
# runDTM.R: top level script file for RSV DTM control		 	#
# developed by: Mohammadreza Amiri: mohammadreza.amiri@iconplc.com#
# editted by: Tmiothy Baker: timothy.baker@iconplc.com		#
# version for model shell 3.0, June23-July3 2023			#
###################################################################

###################################################################
# initialize R environment, load libraries and set options		#
###################################################################

rm(list = ls())

# load the required libraries
library(deSolve)
library(Rcpp)
library(RcppArmadillo)

# Set path
# the args statement passes specifics of the conditions under which the model should be run;
# use the following specific assignment to test the R code outside the model shell
# the sattement must mimic the command line call that the shell uses, specifying the path and options for deterministic or probabilistic analysis 
# args <- c("C:/Users/bakerti/OneDrive - ICON CLINICAL RESEARCH LTD - IRELAND/Desktop/copy_20230620","2","2")

# the following statement sets the options for R according to those passed by the command line call, letting the R model run as part of the shell operation
args = commandArgs(trailingOnly=TRUE)

args<-c("C:/Users/danri/Downloads/RSV_DTM version 4.0 final/","0","0")
#C:\Users\danri\Downloads\RSV_DTM_v3.3 Current shell and model code
#cmd.exe /c C:\"Program Files"\R\R-4.3.1\bin\Rscript.exe "C:\Users\danri\Downloads\RSV_DTM_v3.3 Current shell and model code\runDTM.R" "C:/Users/danri/Downloads/RSV_DTM_v3.3 Current shell and model code/" 0 0
print(args)
setwd(args[1])

###################################################################"
# load the component functions of the DTM by source location	#
###################################################################

source("deSolve/utils_edited.R")
source("deSolve/complete_ode_initial-deSolve2_edited.R")
source("deSolve/complete_ode-deSolve2_edited.R")
source("deSolve/ode_system_edited.R")
sourceCpp("Rcpp/foi.cpp")


# clear any stored results from R overhead				#
times <- NULL   

# confirm starting population file for read-in
path.init.pop <- "data/startpop1.csv"

# diagnostic code, not to be activated
# Monthly_incidence <- colSums(odeOutput)
# plot(times[-1],Monthly_incidence,type = "l",xlab = "time (months)",ylab = "incidence",main = "Total incidence")


###################################################################
# call DTM and subsidiary functions to execute requested interations	#
###################################################################

s <- as.numeric(args[2])

for (vac_type in 1:2) {
    if(args[2] == 0){
        odeOutput <- runDTMX.deSolve2(args = args, vac_type = vac_type, 
                                      init.pop = path.init.pop,
                                      APPEND = F, method = "euler",
                                      ProgressBar = T, Write.Output = T)
        Monthly_incidence <- colSums(odeOutput)
    } else {
        for (x in 1:s) {
            args[2] <- x
            odeOutput <- runDTMX.deSolve2(args = args, vac_type = vac_type, 
                                          init.pop = path.init.pop,
                                          APPEND = T, method = "euler",
                                          ProgressBar = T, Write.Output = T)
            Monthly_incidence <- colSums(odeOutput)
        }
    }
}

###################################################################
# DTM has completed, write notice file to drive to notify shell	#
# that import can proceed							#
###################################################################

tmp_file <- paste0("temp_file_", args[3], ".txt")
writeLines("R analyses completed.", tmp_file)

###################################################################
# end of file									#
###################################################################