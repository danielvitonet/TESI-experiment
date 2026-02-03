
# Forecasting natural gas prices in real time

# Define the list of required packages
required_packages <- c("here","R.matlab","dplyr")
# Check each package; if it's not installed, install it
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    renv::install(pkg)
    #renv::snapshot()
    renv::status()
  }
}
# Load the libraries
lapply(required_packages, library, character.only = TRUE)
# Source several helper functions 
set.seed(123) #set a seed to allow reproducibility

#############################################################################
# MF-BAVART SET-UP (ENVIROMENT and DEPENDENCES)                             #                                                                           
# 1) Path alla repo "vendorizzata" dentro ext/                              #
mf_dir <- here("ext", "mf-bavart")                                          #
#
# 2) Crea un environment dedicato per non sporcare il Global Env            #
#mf <- new.env(parent = baseenv())                                          #
mf <- new.env(parent = globalenv())                                         #
#
# 3) Carica i file R nell’environment mf                                    #
# sys.source() è preferibile quando vuoi specificare l’envir esplicitamente #
sys.source(file.path(mf_dir, "aux_func.R"),      envir = mf, chdir = TRUE)  #
sys.source(file.path(mf_dir, "mfbavart_func.R"), envir = mf, chdir = TRUE)  #
#############################################################################

###################################################################################
# FUNCTIONS NEEDED:                                                               #
# 1) lagn function                                                                #
lagn <- function(data, m) {                                                       #
  # input: data matrix and lag m                                                  #
  data[(m+1):nrow(data), , drop = FALSE] - data[1:(nrow(data)-m), , drop = FALSE] #
}                                                                                 #
###################################################################################

# Load real-time datasets
# Rows: T (starting in 1973M1+36)
# Columns: point in real time (1991M1 to 2024M2)

NG_HENRY <- as.matrix(read.table("C:/Users/HP/Downloads/DATA - LNG/DATA - models/NG_HENRY.txt", sep="\t", header=FALSE))	         
# nominal gas price, Henry Hub
CPI_AC <- as.matrix(read.table("C:/Users/HP/Downloads/DATA - LNG/DATA - models/CPI_AC.txt", sep="\t", header=FALSE))		           
# average change nowcast of US CPI

# Economic predictor variables:
IPALF <- as.matrix(read.table("C:/Users/HP/Downloads/DATA - LNG/DATA - models/IPALF.txt", sep="\t", header=FALSE))		             
# industrial production index (ALFRED vintages)
CAPUALF <- as.matrix(read.table("C:/Users/HP/Downloads/DATA - LNG/DATA - models/CAPUALF.txt", sep="\t", header=FALSE))		         
# US capacity utilization rate (ALFRED vintages)

PROD_DRY_SA <- as.matrix(read.table("C:/Users/HP/Downloads/DATA - LNG/DATA - models/PROD_DRY_SA.txt", sep="\t", header=FALSE))	   
# marketed NG production
STORE_WORK_SA <- as.matrix(read.table("C:/Users/HP/Downloads/DATA - LNG/DATA - models/STORE_WORK_SA.txt", sep="\t", header=FALSE)) 
# underground storage of working gas
CONS_TOT <- as.matrix(read.table("C:/Users/HP/Downloads/DATA - LNG/DATA - models/CONS_TOT.txt", sep="\t", header=FALSE))	         
# total NG consumption
RIGS <- as.matrix(read.table("C:/Users/HP/Downloads/DATA - LNG/DATA - models/RIGS.txt", sep="\t", header=FALSE))		               
# rig count

CDDdev <- as.matrix(read.table("C:/Users/HP/Downloads/DATA - LNG/DATA - models/CDDdev.txt", sep="\t", header=FALSE))		           
# cooling degree days in deviation from historical average
HDDdev <- as.matrix(read.table("C:/Users/HP/Downloads/DATA - LNG/DATA - models/HDDdev.txt", sep="\t", header=FALSE))		           
# heating degree days in deviation from historical average

# GAS Futures start in April 1990 (1990M4) but next year gas futures start in June 1990 (1990M6)
gas_futures <- as.matrix(read.table("C:/Users/HP/Downloads/DATA - LNG/DATA - models/gas_futures.txt", header=FALSE))
fut <- rbind(matrix(NA, nrow = 207, ncol = 9), gas_futures)

# May 2024 vintage of Henry Hub spot price and CPI (final release data)
# rows: T (1973M1-2024M5)
# columns: most recent vintage of 2024M5
# non-NaN rows start on 1997M1 (row 289)
# C:/Users/HP/Downloads/DATA - LNG/HH_CPI_May2024vintage.mat
HH_CPI <- readMat("C:/Users/HP/Downloads/DATA - LNG/DATA - models/HH_CPI_May2024vintage.mat",
                  fixNames = FALSE,          # conserva nomi MATLAB (underscore)
                  drop = "singletonLists",
                  verbose = TRUE)

CPI_May24 <- as.matrix(HH_CPI$CPI_May24)
NG_May24  <- as.matrix(HH_CPI$NG_May24)


# Basic parameters
#ind <- 288+1   # indicates length of initial real-time sample (up to 1997.1)
ind <- 437+1   # indicates length of initial real-time sample (up to 2009.6)
h <- 1         # forecast horizon
p <- 6         # fixed lag order

#jx = 73   
# 12*6+1 = 73 (6 years and 1 month), is the column representing the 1997M1 data vintage for when we do not use future contracts data.
#jx = 222  
# 222 is the column representing the 2009M6 data vintage for which at least 228 observations for future contracts data is available.
# change 37 with 210 , and change 36 with 209 , when using future contracts data

for (jx in (12*18+6):(ncol(CPI_AC)-h)) {  # adjust for evaluation period 2009M6 to 2024M2
  print(jx)
  
  # Create VAR data in column format 
  rpg <- log(100 * NG_HENRY[210:ind, jx] / CPI_AC[210:ind, jx])  # log Real gas price (nominal deflated by US CPI)
  
  # variables in levels 
  capu <- CAPUALF[210:ind, jx]
  hd <- HDDdev[210:ind, jx]
  cd <- CDDdev[210:ind, jx]
  
  # variables in logs
  conslog <- log(CONS_TOT[210:ind, jx])
  
  # variables in growth rates
  ipg <- 100 * lagn(log(IPALF[1:ind, jx, drop = FALSE]), 1)
  dryprod <- 100 * lagn(log(PROD_DRY_SA[1:ind, jx, drop = FALSE]), 1)
  inventories <- 100 * lagn(log(STORE_WORK_SA[1:ind, jx, drop = FALSE]), 1)
  consg <- 100 * lagn(log(CONS_TOT[1:ind, jx, drop = FALSE]), 1)
  rigcount <- 100 * lagn(log(RIGS[1:ind, jx, drop = FALSE]), 1)
 
  
  # one observation lost due to differencing
  ip <- ipg[209:nrow(ipg), 1]
  prod <- dryprod[209:nrow(dryprod), 1]
  store <- inventories[209:nrow(inventories), 1]
  cons <- consg[209:nrow(consg), 1]
  rigs <- rigcount[209:nrow(rigcount), 1]

  
  # futures contract variable
  futgasrt <- log(fut[210:ind, ])
  spotgasrt <- log(NG_HENRY[210:ind, jx])
  #inflrt <- log(CPI_AC[(ind - 120 + 2):ind, jx]) - log(CPI_AC[(ind - 120 + 1):(ind - 1), jx]) # average inflation over the past 10 years
  inflrt <- log(CPI_AC[210:ind, jx]) - log(CPI_AC[209:(ind - 1), jx]) # monthly inflation
  
  jj <- switch(as.character(h),
               "1" = 1,
               "3" = 2,
               "6" = 3,
               "9" = 4,
               "12" = 5,
               "15" = 6,
               "18" = 7,
               "21" = 8,
               "24" = 9)
  
  #futs <- futgasrt[, jj] - spotgasrt - ((1 + mean(inflrt))^h - 1)
  infl_ma <- stats::filter(inflrt, rep(1 / 120, 120), sides = 1) # 10-year rolling average
  infl_exp <- (1 + infl_ma)^h - 1
  #futs <- ts(futgasrt[, jj] - spotgasrt - as.numeric(infl_exp), frequency = 12)
  futs <- futgasrt[, jj] - spotgasrt - as.numeric(infl_exp)
  
  # Create revised real price of natural gas (most recent vintage)
  x <- 100 * NG_May24[210:(ind + h), 1] / CPI_May24[210:(ind + h), 1]  # Real gas price (nominal deflated by US CPI)
  
  # Estimate BAVART
  # (Estimation code needed here)
  data <- list(
    log_GAS_Price = ts(rpg, frequency = 12),             # monthly
    mgr_DRY_Production = ts(prod, frequency = 12),       # monthly
    mgr_Working_Inventories = ts(store, frequency = 12), # monthly
    #CONS = ts(cons, frequency = 12),                     # monthly
    mgr_Rig_counts = ts(rigs, frequency = 12),           # monthly
    mgr_INDUSTRIAL_PRODUCTION = ts(ip, frequency = 12),  # monthly
    #CAP_UT = ts(capu, frequency = 12),                   # monthly
    HDD_dev = ts(hd, frequency = 12),                    # monthly
    #CDD_dev = ts(cd, frequency = 12),                    # monthly
    Fut_sp = ts(futs, frequency = 12)                    # monthly
  )
  
  # qgr_Real_GDP = ts(..., frequency = 4)                # quarterly


  # Estimate BVAR from fatBVARS package BVAR.novol(...)
  # (Estimation code needed here)

  # Evaluate h-step ahead BAVAR forecast
  # (Forecast evaluation code needed here)
  # save results

  # Estimate MF-BAVART-SV from mf$mfbavart(...)
  # (Estimation code needed here)
  
  # Evaluate h-step ahead MF-BAVART-SV forecast
  # (Forecast evaluation code needed here)
  # save results
  
  # Estimate mixBART from flexBART(...)
  # (Estimation code needed here)
  
  # Evaluate h-step ahead mixBART forecast
  # (Forecast evaluation code needed here)
  # save results
  
  # Evaluate h-step ahead RW forecast
  sfe[jx - 12 * 6, 2] <- (exp(rpg[length(rpg)]) - x[t + h])^2  # Squared forecast error
  
  # Keep track of forecasts
  # (Code needed here)
  
  # Update index for recursive estimation
  ind <- ind + 1
}

# Evaluate real-time recursive forecast accuracy
# (Code needed here)

# Compute Average Quantile Scores
# (Code needed here)

# Compute quantile-weighted Continuous Ranked Probability Scores (qw-CRPS)
# (Code needed here)

# Compute Directional Symmetry (DS) statistic
# (Code needed here)
