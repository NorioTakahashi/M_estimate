#' R code for M (natural mortality) estimation for north Pacicfic chub mackerel 
#'
#' This R code was used to estimate M using emirical estimators and
#' "FishLife" (Thorson et al 2017) for north Pacific chub mackerel.
#' This work was done for NPFC 2nd Meeting of the Technical Working
#' Group on Chub Mackerel Stock Assessment.
#' 
#' @author norio takahashi
#' 
#' @param arg1 description
#' @importFrom pkg func
#' @return description
#' @examples code
#' @export
#'


setwd("E:/norio_top/国内_資源評価2018/M推定2018/NPFC2019_TWG_CMSA02対応作業")

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


#### parameters values

max_age_2006_18_all <- 10  # Kamimura et al
max_age_1960_70     <- 11  # Iiduka (2002)

temp_mean <- 16.7  # f/ data provided by Kamimura
temp_min  <- 15.1  # f/ data provided by Kamimura
temp_max  <- 18.9  # f/ data provided by Kamimura

size_len_mod <- 29.  # in cm
size_len_min <- 5.   # in cm
size_len_max <- 48.  # in cm


################ results for paper draft 1 (start) ################

# read VB Linf and K for each cohort (Kamimura et al)
dVB <- read.csv("VB_Linf_K_cohort.csv")
dVB$Linf <- dVB$Linf/10. # must convert mm -> cm

df <- data.frame(estimator="dummy", M=NA, Linf=NA, K=NA, cohort=NA, max_age=NA, length=NA, temp=NA)

#### Pauly1 w/ Kamimura VB parameters and mean, min, and max water temp
for(ii in seq(dVB$cohort)){
  
  Mvalue <- Calc_M(whose="pauly1_L", Linf=dVB$Linf[ii], KK=dVB$K[ii], temp=temp_mean)
  
  wdf <- data.frame(estimator="Pauly", M=Mvalue, Linf=dVB$Linf[ii], K=dVB$K[ii], cohort=dVB$cohort[ii],
                    max_age=NA, length=NA, temp=temp_mean)
  df <- rbind(df, wdf)
  
}

for(ii in seq(dVB$cohort)){
  
  Mvalue <- Calc_M(whose="pauly1_L", Linf=dVB$Linf[ii], KK=dVB$K[ii], temp=temp_min)
  
  wdf <- data.frame(estimator="Pauly", M=Mvalue, Linf=dVB$Linf[ii], K=dVB$K[ii], cohort=dVB$cohort[ii],
                    max_age=NA, length=NA, temp=temp_min)
  df <- rbind(df, wdf)
  
}

for(ii in seq(dVB$cohort)){
  
  Mvalue <- Calc_M(whose="pauly1_L", Linf=dVB$Linf[ii], KK=dVB$K[ii], temp=temp_max)
  
  wdf <- data.frame(estimator="Pauly", M=Mvalue, Linf=dVB$Linf[ii], K=dVB$K[ii], cohort=dVB$cohort[ii],
                    max_age=NA, length=NA, temp=temp_max)
  df <- rbind(df, wdf)
  
}


#### Pauly (updated) w/ Kamimura VB parameters
for(ii in seq(dVB$cohort)){
  
  Mvalue <- Calc_M(whose="pauly_L_update", Linf=dVB$Linf[ii], KK=dVB$K[ii])
  
  wdf <- data.frame(estimator="Pauly_update", M=Mvalue, Linf=dVB$Linf[ii], K=dVB$K[ii], cohort=dVB$cohort[ii],
                    max_age=NA, length=NA, temp=NA)
  df <- rbind(df, wdf)
  
}


#### Jensen2 w/ Kamimura VB parameters
for(ii in seq(dVB$cohort)){
  
  Mvalue <- Calc_M(whose="jensen2", KK=dVB$K[ii])
  
  wdf <- data.frame(estimator="Jensen", M=Mvalue, Linf=NA, K=dVB$K[ii], cohort=dVB$cohort[ii],
                    max_age=NA, length=NA, temp=NA)
  df <- rbind(df, wdf)
  
}


#### Hoenig w/ Kamimura max_age (for all cohort combined)
Mvalue <- Calc_M(whose="hoenig", max_age=max_age_2006_18_all)

wdf <- data.frame(estimator="Hoenig", M=Mvalue, Linf=NA, K=NA, cohort=NA,
                  max_age=max_age_2006_18_all, length=NA, temp=NA)
df <- rbind(df, wdf)


#### Hoenig w/ Iiduka (2002) max_age
Mvalue <- Calc_M(whose="hoenig", max_age=max_age_1960_70)

wdf <- data.frame(estimator="Hoenig", M=Mvalue, Linf=NA, K=NA, cohort=NA,
                  max_age=max_age_1960_70, length=NA, temp=NA)
df <- rbind(df, wdf)


#### Hoenig (updated) w/ Kamimura max_age (for all cohort combined)
Mvalue <- Calc_M(whose="hoenig_update", max_age=max_age_2006_18_all)

wdf <- data.frame(estimator="Hoenig_update", M=Mvalue, Linf=NA, K=NA, cohort=NA,
                  max_age=max_age_2006_18_all, length=NA, temp=NA)
df <- rbind(df, wdf)


#### Hoenig (updated) w/ Iiduka (2002) max_age
Mvalue <- Calc_M(whose="hoenig_update", max_age=max_age_1960_70)

wdf <- data.frame(estimator="Hoenig_update", M=Mvalue, Linf=NA, K=NA, cohort=NA,
                  max_age=max_age_1960_70, length=NA, temp=NA)
df <- rbind(df, wdf)


#### Gislason1 w/ Kamimura VB parameters and modal, min, and Linf length
for(ii in seq(dVB$cohort)){
  
  Mvalue <- Calc_M(whose="gislason1", Linf=dVB$Linf[ii], KK=dVB$K[ii], size_len=size_len_mod)
  
  wdf <- data.frame(estimator="Gislason1", M=Mvalue, Linf=dVB$Linf[ii], K=dVB$K[ii], cohort=dVB$cohort[ii],
                    max_age=NA, length=size_len_mod, temp=NA)
  df <- rbind(df, wdf)
  
}

for(ii in seq(dVB$cohort)){
  
  Mvalue <- Calc_M(whose="gislason1", Linf=dVB$Linf[ii], KK=dVB$K[ii], size_len=size_len_min)
  
  wdf <- data.frame(estimator="Gislason1", M=Mvalue, Linf=dVB$Linf[ii], K=dVB$K[ii], cohort=dVB$cohort[ii],
                    max_age=NA, length=size_len_min, temp=NA)
  df <- rbind(df, wdf)
  
}

for(ii in seq(dVB$cohort)){
  
  #Mvalue <- Calc_M(whose="gislason1", Linf=dVB$Linf[ii], KK=dVB$K[ii], size_len=size_len_max)
  Mvalue <- Calc_M(whose="gislason1", Linf=dVB$Linf[ii], KK=dVB$K[ii], size_len=dVB$Linf[ii])
  
  #wdf <- data.frame(estimator="Gislason1", M=Mvalue, Linf=dVB$Linf[ii], K=dVB$K[ii], cohort=dVB$cohort[ii],
  #                  max_age=NA, length=size_len_max, temp=NA)
  wdf <- data.frame(estimator="Gislason1", M=Mvalue, Linf=dVB$Linf[ii], K=dVB$K[ii], cohort=dVB$cohort[ii],
                    max_age=NA, length=dVB$Linf[ii], temp=NA)
  df <- rbind(df, wdf)
  
}


#### Gislason2 w/ Kamimura VB parameters and modal, min, and Linf length
for(ii in seq(dVB$cohort)){
  
  Mvalue <- Calc_M(whose="gislason2", Linf=dVB$Linf[ii], KK=dVB$K[ii], size_len=size_len_mod)
  
  wdf <- data.frame(estimator="Gislason2", M=Mvalue, Linf=dVB$Linf[ii], K=dVB$K[ii], cohort=dVB$cohort[ii],
                    max_age=NA, length=size_len_mod, temp=NA)
  df <- rbind(df, wdf)
  
}

for(ii in seq(dVB$cohort)){
  
  Mvalue <- Calc_M(whose="gislason2", Linf=dVB$Linf[ii], KK=dVB$K[ii], size_len=size_len_min)
  
  wdf <- data.frame(estimator="Gislason2", M=Mvalue, Linf=dVB$Linf[ii], K=dVB$K[ii], cohort=dVB$cohort[ii],
                    max_age=NA, length=size_len_min, temp=NA)
  df <- rbind(df, wdf)
  
}

for(ii in seq(dVB$cohort)){
  
  #Mvalue <- Calc_M(whose="gislason2", Linf=dVB$Linf[ii], KK=dVB$K[ii], size_len=size_len_max)
  Mvalue <- Calc_M(whose="gislason2", Linf=dVB$Linf[ii], KK=dVB$K[ii], size_len=dVB$Linf[ii])
  
  #wdf <- data.frame(estimator="Gislason2", M=Mvalue, Linf=dVB$Linf[ii], K=dVB$K[ii], cohort=dVB$cohort[ii],
  #                  max_age=NA, length=size_len_max, temp=NA)
  wdf <- data.frame(estimator="Gislason2", M=Mvalue, Linf=dVB$Linf[ii], K=dVB$K[ii], cohort=dVB$cohort[ii],
                    max_age=NA, length=dVB$Linf[ii], temp=NA)
  df <- rbind(df, wdf)
  
}


#### Thorson (original, not updated w/ Kamimura et al. data)
Mvalue <- Calc_M(whose="thorson")

wdf <- data.frame(estimator="Thorson", M=Mvalue, Linf=NA, K=NA, cohort=NA,
                  max_age=NA, length=NA, temp=NA)
df <- rbind(df, wdf)


#### write results to csv file
df <- df[-1,]
write.csv(df, file="M_estimate.csv", row.names=FALSE)

################ results for paper draft 1 (end) ################



################ results for paper draft 2 (start) ################

#### parameter values

Linf_all <- 446.25/10.  # must be converted mm to cm
K_all <- 0.20
t0_all <- -3.05

Linf_200612 <- 394.99/10.  # must be converted mm to cm
K_200612 <- 0.42

Linf_201315 <- 327.82/10.  # must be converted mm to cm
K_201315 <- 0.51

df2 <- data.frame(estimator="dummy", M=NA, Linf=NA, K=NA, cohort=NA, max_age=NA, length=NA, temp=NA)

#### Pauly1 w/ Kamimura VB parameters (all, 2006-12, 2013-15) and mean water temp
Mvalue <- Calc_M(whose="pauly1_L", Linf=Linf_all, KK=K_all, temp=temp_mean)

wdf <- data.frame(estimator="Pauly", M=Mvalue, Linf=Linf_all, K=K_all, cohort=200615,
                  max_age=NA, length=NA, temp=temp_mean)
df2 <- rbind(df2, wdf)

Mvalue <- Calc_M(whose="pauly1_L", Linf=Linf_200612, KK=K_200612, temp=temp_mean)

wdf <- data.frame(estimator="Pauly", M=Mvalue, Linf=Linf_200612, K=K_200612, cohort=200612,
                  max_age=NA, length=NA, temp=temp_mean)
df2 <- rbind(df2, wdf)

Mvalue <- Calc_M(whose="pauly1_L", Linf=Linf_201315, KK=K_201315, temp=temp_mean)

wdf <- data.frame(estimator="Pauly", M=Mvalue, Linf=Linf_201315, K=K_201315, cohort=201315,
                  max_age=NA, length=NA, temp=temp_mean)
df2 <- rbind(df2, wdf)


#### Pauly (updated) w/ Kamimura VB parameters (all, 2006-12, 2013-15)
Mvalue <- Calc_M(whose="pauly_L_update", Linf=Linf_all, KK=K_all)

wdf <- data.frame(estimator="Pauly_update", M=Mvalue, Linf=Linf_all, K=K_all, cohort=200615,
                  max_age=NA, length=NA, temp=NA)
df2 <- rbind(df2, wdf)

Mvalue <- Calc_M(whose="pauly_L_update", Linf=Linf_200612, KK=K_200612)

wdf <- data.frame(estimator="Pauly_update", M=Mvalue, Linf=Linf_200612, K=K_200612, cohort=200612,
                  max_age=NA, length=NA, temp=NA)
df2 <- rbind(df2, wdf)

Mvalue <- Calc_M(whose="pauly_L_update", Linf=Linf_201315, KK=K_201315)

wdf <- data.frame(estimator="Pauly_update", M=Mvalue, Linf=Linf_201315, K=K_201315, cohort=201315,
                  max_age=NA, length=NA, temp=NA)
df2 <- rbind(df2, wdf)


#### Jensen2 w/ Kamimura VB parameters (all, 2006-12, 2013-15)
Mvalue <- Calc_M(whose="jensen2", KK=K_all)

wdf <- data.frame(estimator="Jensen", M=Mvalue, Linf=NA, K=K_all, cohort=200615,
                  max_age=NA, length=NA, temp=NA)
df2 <- rbind(df2, wdf)

Mvalue <- Calc_M(whose="jensen2", KK=K_200612)

wdf <- data.frame(estimator="Jensen", M=Mvalue, Linf=NA, K=K_200612, cohort=200612,
                  max_age=NA, length=NA, temp=NA)
df2 <- rbind(df2, wdf)

Mvalue <- Calc_M(whose="jensen2", KK=K_201315)

wdf <- data.frame(estimator="Jensen", M=Mvalue, Linf=NA, K=K_201315, cohort=201315,
                  max_age=NA, length=NA, temp=NA)
df2 <- rbind(df2, wdf)


#### Hoenig w/ Kamimura max_age (for all cohort combined)
Mvalue <- Calc_M(whose="hoenig", max_age=max_age_2006_18_all)

wdf <- data.frame(estimator="Hoenig", M=Mvalue, Linf=NA, K=NA, cohort=NA,
                  max_age=max_age_2006_18_all, length=NA, temp=NA)
df2 <- rbind(df2, wdf)


#### Hoenig w/ Iiduka (2002) max_age
Mvalue <- Calc_M(whose="hoenig", max_age=max_age_1960_70)

wdf <- data.frame(estimator="Hoenig", M=Mvalue, Linf=NA, K=NA, cohort=NA,
                  max_age=max_age_1960_70, length=NA, temp=NA)
df2 <- rbind(df2, wdf)


#### Hoenig (updated) w/ Kamimura max_age (for all cohort combined)
Mvalue <- Calc_M(whose="hoenig_update", max_age=max_age_2006_18_all)

wdf <- data.frame(estimator="Hoenig_update", M=Mvalue, Linf=NA, K=NA, cohort=NA,
                  max_age=max_age_2006_18_all, length=NA, temp=NA)
df2 <- rbind(df2, wdf)


#### Hoenig (updated) w/ Iiduka (2002) max_age
Mvalue <- Calc_M(whose="hoenig_update", max_age=max_age_1960_70)

wdf <- data.frame(estimator="Hoenig_update", M=Mvalue, Linf=NA, K=NA, cohort=NA,
                  max_age=max_age_1960_70, length=NA, temp=NA)
df2 <- rbind(df2, wdf)


#### Gislason1 w/ Kamimura VB parameters (all, 2006-12, 2013-15) and modal length
Mvalue <- Calc_M(whose="gislason1", Linf=Linf_all, KK=K_all, size_len=size_len_mod)

wdf <- data.frame(estimator="Gislason1", M=Mvalue, Linf=Linf_all, K=K_all, cohort=200615,
                  max_age=NA, length=size_len_mod, temp=NA)
df2 <- rbind(df2, wdf)

Mvalue <- Calc_M(whose="gislason1", Linf=Linf_200612, KK=K_200612, size_len=size_len_mod)

wdf <- data.frame(estimator="Gislason1", M=Mvalue, Linf=Linf_200612, K=K_200612, cohort=200612,
                  max_age=NA, length=size_len_mod, temp=NA)
df2 <- rbind(df2, wdf)

Mvalue <- Calc_M(whose="gislason1", Linf=Linf_201315, KK=K_201315, size_len=size_len_mod)

wdf <- data.frame(estimator="Gislason1", M=Mvalue, Linf=Linf_201315, K=K_201315, cohort=201315,
                  max_age=NA, length=size_len_mod, temp=NA)
df2 <- rbind(df2, wdf)


#### Gislason2 w/ Kamimura VB parameters (all, 2006-12, 2013-15) and modal length
Mvalue <- Calc_M(whose="gislason2", Linf=Linf_all, KK=K_all, size_len=size_len_mod)

wdf <- data.frame(estimator="Gislason2", M=Mvalue, Linf=Linf_all, K=K_all, cohort=200615,
                  max_age=NA, length=size_len_mod, temp=NA)
df2 <- rbind(df2, wdf)

Mvalue <- Calc_M(whose="gislason2", Linf=Linf_200612, KK=K_200612, size_len=size_len_mod)

wdf <- data.frame(estimator="Gislason2", M=Mvalue, Linf=Linf_200612, K=K_200612, cohort=200612,
                  max_age=NA, length=size_len_mod, temp=NA)
df2 <- rbind(df2, wdf)

Mvalue <- Calc_M(whose="gislason2", Linf=Linf_201315, KK=K_201315, size_len=size_len_mod)

wdf <- data.frame(estimator="Gislason2", M=Mvalue, Linf=Linf_201315, K=K_201315, cohort=201315,
                  max_age=NA, length=size_len_mod, temp=NA)
df2 <- rbind(df2, wdf)


#### Thorson (original, not updated w/ Kamimura et al. data)
Mvalue <- Calc_M(whose="thorson")

wdf <- data.frame(estimator="Thorson", M=Mvalue, Linf=NA, K=NA, cohort=NA,
                  max_age=NA, length=NA, temp=NA)
df2 <- rbind(df2, wdf)


#### write results to csv file
df2 <- df2[-1,]
write.csv(df2, file="M_estimate_rev1.csv", row.names=FALSE)


#### M by age using gislason1 and gislason2 and Kamimura (all)
df3 <- data.frame(estimator="dummy", M=NA, Linf=NA, K=NA, cohort=NA, age=NA, length=NA)

vonBertalanffy(Linf=Linf_all, KK=K_all, t0=t0_all, age=0)  # length at age 0?

# gislason1
for(iage in 1:6){
  
  La[iage] <- vonBertalanffy(Linf=Linf_all, KK=K_all, t0=t0_all, age=iage)
  
  Mvalue <- Calc_M(whose="gislason1", Linf=Linf_all, KK=K_all, size_len=La[iage])
  
  wdf <- data.frame(estimator="Gislason1", M=Mvalue, Linf=Linf_all, K=K_all, cohort=200615,
                    age=iage, length=La[iage])
  df3 <- rbind(df3, wdf)
  
}

# gislason2
for(iage in 1:6){
  
  Mvalue <- Calc_M(whose="gislason2", Linf=Linf_all, KK=K_all, size_len=La[iage])
  
  wdf <- data.frame(estimator="Gislason2", M=Mvalue, Linf=Linf_all, K=K_all, cohort=200615,
                    age=iage, length=La[iage])
  df3 <- rbind(df3, wdf)
  
}

# write to csv file
df3 <- df3[-1,]
write.csv(df3, file="M_estimate_byage_gislason_all.csv", row.names=FALSE)

par(mfrow=c(1,1))
plot(rep(1:6,2), df3$M, xlab="Age", ylab="M", type="n")
lines(1:6, df3$M[df3$estimator=="Gislason1"], type="b", col="blue", pch="1", lwd=2)
lines(1:6, df3$M[df3$estimator=="Gislason2"], type="b", col="green", pch="2", lwd=2)
legend(4.0, 0.45, unique(df3$estimator), pch=c("1","2"), lty=rep(1,2), col=c("blue","green"), lwd=2, cex=1.5)



Predict.org <- Plot_taxa( 
  Search_species(Genus="Scomber",Species="japonicus")$match_taxonomy,
  mfrow=c(2,2) )

Predict.org[[1]]$Mean_pred
log_M <- as.numeric(Predict.org[[1]]$Mean_pred["M"])

Predict.org [[1]]$Cov_pred  #分散共分散 (temperature以外はlogscale)
sqrt(diag(Predict.org[[1]]$Cov_pred))  #SD (CV)
log_SD <- as.numeric(sqrt(diag(Predict.org[[1]]$Cov_pred))["M"])

# mean M +- 1SD
low_M <- exp(log_M-log_SD)
upp_M <- exp(log_M+log_SD)


Ynew_ij <- matrix( c("Loo"=log(Linf_all),
                     "K"=log(K_all),
                     "Winfinity"=NA,
                     "tmax"=log(max_age_2006_18_all),
                     "tm"=NA,
                     "M"=NA,
                     "Lm"=NA,
                     "Temperature"=temp_mean),
                   nrow=1)

Update.mean <- Update_prediction( 
  Taxon=Search_species(Genus="Scomber",Species="japonicus",add_ancestors=FALSE)$match_taxonomy,
  Ynew_ij=Ynew_ij)





################ results for paper draft 2 (end) ################



#### how to use Thorson's package
vignette("tutorial","FishLife")

# original (not updated)
Predict.org <- Plot_taxa( 
  Search_species(Genus="Scomber",Species="japonicus")$match_taxonomy,
  mfrow=c(2,2) )

Predict.org[[1]]$Mean_pred

# update with Kamimura et al data
Ynew_ij <- matrix( c("Loo"=log(Linf_2006_18_all),
                     "K"=log(KK_2006_18_all),
                     "Winfinity"=NA,
                     "tmax"=log(max_age_2006_18_all),
                     "tm"=NA,
                     "M"=NA,
                     "Lm"=NA,
                     "Temperature"=temp_mean),
                   nrow=1)

Update.mean <- Update_prediction( 
  Taxon=Search_species(Genus="Scomber",Species="japonicus",add_ancestors=FALSE)$match_taxonomy,
  Ynew_ij=Ynew_ij)



#### for check
#Linf <- 12.93
#KK <- 0.23
#temp <- 17
#max_age <- 3
#size_len1 <- 3
#size_len2 <- 10
#
#Calc_M(whose="norio")
#Calc_M(whose="tanaka", max_age=max_age)
#Calc_M(whose="pauly1_L", Linf=Linf, KK=KK, temp=temp)
#Calc_M(whose="pauly_L_update", Linf=Linf, KK=KK)
#Calc_M(whose="jensen2", KK=KK)
#Calc_M(whose="hoenig", max_age=max_age)
#Calc_M(whose="hoenig_update", max_age=max_age)
#Calc_M(whose="gislason1", Linf=Linf, KK=KK, size_len=size_len1)
#Calc_M(whose="gislason1", Linf=Linf, KK=KK, size_len=size_len2)
#Calc_M(whose="gislason2", Linf=Linf, KK=KK, size_len=size_len1)
#Calc_M(whose="gislason2", Linf=Linf, KK=KK, size_len=size_len2)
#Calc_M(whose="thorson")
####

####
