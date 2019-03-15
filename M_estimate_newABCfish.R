#' R code for M (natural mortality) estimation for fish stocks for which new ABC rules are applied 
#'
#' This R code was used to estimate M using emirical estimators and
#' "FishLife" (Thorson et al 2017) for fish stocks for which new ABC rules are
#' applied.
#' 
#' @author norio takahashi
#' 
#' @param arg1 description
#' @importFrom pkg func
#' @return description
#' @examples code
#' @export
#'


#setwd("E:/norio_top/国内_資源評価2018/M推定2018/M_estimate")

#### load tidyverse package
library(tidyverse)

#### Install and load FishLife package
#devtools::install_github("james-thorson/FishLife")
library(FishLife)

#### load TMB package for FishLife (when re-run or update)
library(TMB)

remove(list=ls())

#### function M estimator function
source("M_estimator.R")

#### define function for von Bertalanffy growth curve
vonBertalanffy <- function(
  Linf = 0,
  KK = 0,
  t0 = 0,
  age = 0
){
  
  La <- Linf * ( 1 - exp(-KK*(age-t0)) )
  
  return(La)
  
}


#### biological parameters values

guess_encoding("bio_param.csv")

bparam <- read_csv("bio_param.csv", locale = locale(encoding = "Shift_JIS"))
head(bparam)
str(bparam)
glimpse(bparam)


df <- tibble(
  
  jname = "dummy",
  stock = "dummy",
  genus = "dummy",
  spp = "dummy",
  Pauly = NA,
  Pauly_update = NA,
  Jensen = NA,
  Hoenig = NA,
  Hoenig_update = NA,
  Gislason1 = NA,
  Gislason2 = NA,
  FishLife = NA,
  FishLife_lowSD = NA,
  FishLife_uppSD = NA,
  FishLife_estimate_taxon_level = NA,
  median = NA
  
)

for (ii in 1:nrow(bparam)) {
  
  M_pau <- list(M = NA, M_1SD = c(NA, NA), FishLife_estimate_taxon_level = NA)
  M_pau_up <- list(M = NA, M_1SD = c(NA, NA), FishLife_estimate_taxon_level = NA)
  M_jen <- list(M = NA, M_1SD = c(NA, NA), FishLife_estimate_taxon_level = NA)
  M_hoe <- list(M = NA, M_1SD = c(NA, NA), FishLife_estimate_taxon_level = NA)
  M_hoe_up <- list(M = NA, M_1SD = c(NA, NA), FishLife_estimate_taxon_level = NA)
  M_gis1 <- list(M = NA, M_1SD = c(NA, NA), FishLife_estimate_taxon_level = NA)
  M_gis2 <- list(M = NA, M_1SD = c(NA, NA), FishLife_estimate_taxon_level = NA)
  M_tho <- list(M = NA, M_1SD = c(NA, NA), FishLife_estimate_taxon_level = NA)
  median_val <- NA

  if( !is.na(bparam$Linf[ii]) & !is.na(bparam$K[ii]) & !is.na(bparam$Temp[ii]) ){
    
    M_pau <- Calc_M(whose = "pauly1_L", Linf = bparam$Linf[ii], KK = bparam$K[ii], temp = bparam$Temp[ii])

  }
  
  if( !is.na(bparam$Linf[ii]) & !is.na(bparam$K[ii]) ){
    
    M_pau_up <- Calc_M(whose="pauly_L_update", Linf = bparam$Linf[ii], KK = bparam$K[ii])
    
  }
  
  if( !is.na(bparam$K[ii]) ){
    
    M_jen <- Calc_M(whose = "jensen2", KK = bparam$K[ii])
    
  }
  
  if( !is.na(bparam$Amax[ii]) ){
    
    M_hoe <- Calc_M(whose = "hoenig", max_age = bparam$Amax[ii])
    M_hoe_up <- Calc_M(whose = "hoenig_update", max_age = bparam$Amax[ii])
    
  }
  
  if( !is.na(bparam$Linf[ii]) & !is.na(bparam$K[ii]) & !is.na(bparam$Lsize[ii]) ){
    
    M_gis1 <- Calc_M(whose = "gislason1", Linf = bparam$Linf[ii], KK = bparam$K[ii], size_len = bparam$Lsize[ii])
    M_gis2 <- Calc_M(whose = "gislason2", Linf = bparam$Linf[ii], KK = bparam$K[ii], size_len = bparam$Lsize[ii])
    
  }
  
  if( !is.na(bparam$genus[ii]) & !is.na(bparam$spp[ii]) ){
    
    M_tho <- Calc_M(whose = "thorson", genus = bparam$genus[ii], species = bparam$spp[ii])
    
  }
  
  median_val <- median(na.omit(c(
    as.numeric(M_pau["M"]),
    as.numeric(M_pau_up["M"]),
    as.numeric(M_jen["M"]),
    as.numeric(M_hoe["M"]),
    as.numeric(M_hoe_up["M"]),
    as.numeric(M_gis1["M"]),
    as.numeric(M_gis2["M"]),
    as.numeric(M_tho["M"])
  )))
  
  wdf <- tibble(
    
    jname = bparam$jname[ii],
    stock = bparam$stock[ii],
    genus = bparam$genus[ii],
    spp = bparam$spp[ii],
    Pauly = as.numeric(M_pau["M"]),
    Pauly_update = as.numeric(M_pau_up["M"]),
    Jensen = as.numeric(M_jen["M"]),
    Hoenig = as.numeric(M_hoe["M"]),
    Hoenig_update = as.numeric(M_hoe_up["M"]),
    Gislason1 = as.numeric(M_gis1["M"]),
    Gislason2 = as.numeric(M_gis2["M"]),
    FishLife = as.numeric(M_tho["M"]),
    FishLife_lowSD = as.numeric(M_tho[[2]][1]),
    FishLife_uppSD = as.numeric(M_tho[[2]][2]),
    FishLife_estimate_taxon_level = as.character(M_tho["FishLife_estimate_taxon_level"]),
    median = median_val
    
  )
  
  df <- rbind(df, wdf)
  
  df <- df %>%
    filter(jname != "dummy")
  
}

df %>% 
  write_excel_csv("M_estimate_newABCfish.csv") # note!!: encoding for Japanese characters




Predict.org <- FishLife::Plot_taxa( 
  FishLife::Search_species(Genus = "Scomber", Species = "japonicus")$match_taxonomy,
  mfrow=c(2,2) )

Predict.org[[1]]$Mean_pred
log_M <- as.numeric(Predict.org[[1]]$Mean_pred["M"])

Predict.org [[1]]$Cov_pred  #分散共分散 (temperature以外はlogscale)
sqrt(diag(Predict.org[[1]]$Cov_pred))  #SD (CV)
log_SD <- as.numeric(sqrt(diag(Predict.org[[1]]$Cov_pred))["M"])

# mean M +- 1SD
low_M <- exp(log_M-log_SD)
upp_M <- exp(log_M+log_SD)


Ynew_ij <- matrix( c("Loo" = log(44.6),
                     "K" = log(0.20),
                     "Winfinity" = NA,
                     "tmax" = log(10),
                     "tm" = NA,
                     "M" = NA,
                     "Lm" = NA,
                     "Temperature" = 16.7),
                     nrow=1)

Update.mean <- FishLife::Update_prediction( 
  Taxon = FishLife::Search_species(Genus = "Scomber", Species = "japonicus", add_ancestors = FALSE)$match_taxonomy,
  Ynew_ij = Ynew_ij)
