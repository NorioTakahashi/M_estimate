#' function for M (natural mortality rate) estimator
#'
#' this function uses a subset of simple empirical M estimators
#' and "FishLife" (Thorson et al 2017)
#' 
#' @author norio takahashi
#' 
#' @param whose M estimator ID c("tanaka","pauly1_L","pauly_L_update","jensen2","hoenig","hoenig_update","gislason1","gislason2","thorson")
#' @param max_age observed maximum age
#' @param Linf L(length) infinity in cm (growth parameter)
#' @param KK K (growth parameter)
#' @param temp ambient water temperature in degree C
#' @param size_len length in cm
#' 
#' @importFrom FishLife Plot_taxa
#' 
#' @return estimated M value
#'  
#' @examples Calc_M(whose = "pauly1_L", Linf = 200, KK = 0.5, temp = 12)
#' 
#' @export
#'

Calc_M <- function(
  
  whose = "not_specified",
  max_age = 0,
  Linf = 0,  # in cm
  KK = 0,
  temp = 0, # in degree C
  size_len = 0 # in cm
  
  # -- available M estimators --
  # tanaka: Tanaka (1960)
  # pauly1_L: Pauly (1980)
  # pauly_L_update: Then et al. (2015)
  # jensen2: Jensen (1996)
  # hoenig: Hoenig (1983)
  # hoenig_update: Then et al. (2015)
  # gislason1: Gislason (2010)
  # gislason2: Gislason (2010) and Charnov et al. (2013)
  # thorson: Thorson et al. (2017)
  
){
  
  nat_mort <- 0
  
  if(whose=="tanaka"){
    
    if(max_age==0) warning("max_age given as 0 for tanaka")
    
    nat_mort <- 2.5/max_age
    
  }else if(whose=="pauly1_L"){
    
    if(Linf==0) warning("Linf given as 0 for pauly1_L")
    if(KK==0) warning("KK given as 0 for pauly1_L")
    if(temp==0) warning("temp given as 0 for pauly1_L")
    
    nat_mort <- 0.9849 * Linf^(-0.279) * KK^0.6543 * temp^0.4634
    
  }else if(whose=="pauly_L_update"){
    
    if(Linf==0) warning("Linf given as 0 for pauly_L_update")
    if(KK==0) warning("KK given as 0 for pauly_L_update")
    
    nat_mort <- 4.118 * Linf^(-0.33) * KK^0.73
    
  }else if(whose=="jensen2"){
    
    if(KK==0) warning("KK given as 0 for jensen2")
    
    nat_mort <- 1.5*KK
    
  }else if(whose=="hoenig"){
    
    if(max_age==0) warning("max_age given as 0 for hoenig")
    
    nat_mort <- 4.3/max_age
    
  }else if(whose=="hoenig_update"){
    
    if(max_age==0) warning("max_age given as 0 for hoenig_update")
    
    nat_mort <- 4.899 * max_age^(-0.916)
    
  }else if(whose=="gislason1"){
    
    if(size_len==0) warning("size_len given as 0 for gislason1")
    if(Linf==0) warning("Linf given as 0 for gislason1")
    if(KK==0) warning("KK given as 0 for gislason1")
    
    nat_mort <- 1.73 * size_len^(-1.61) *  Linf^(1.44) * KK
    
  }else if(whose=="gislason2"){
    
    if(size_len==0) warning("size_len given as 0 for gislason2")
    if(Linf==0) warning("Linf given as 0 for gislason2")
    if(KK==0) warning("KK given as 0 for gislason2")
    
    nat_mort <- KK * (size_len/Linf)^(-1.5)
    
  }else if(whose=="thorson"){
    
    Predict.org <- Plot_taxa( 
      Search_species(Genus="Scomber",Species="japonicus")$match_taxonomy,
      mfrow=c(2,2) )
    
    nat_mort <- exp(as.numeric(Predict.org[[1]]$Mean_pred["M"]))
    
  }else{
    
    warning("name of M estimator was not specified or was mis-specified")
    
  }
  
  return(nat_mort)
  
}

