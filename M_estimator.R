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
#' @param genus
#' @param species
#' 
#' @importFrom FishLife Plot_taxa
#' 
#' @return estimated M and M+-1SD values 
#'  
#' @examples Calc_M(whose = "pauly1_L", Linf = 200, KK = 0.5, temp = 12)
#' 
#' @export
#'

Calc_M <- function(
  
  whose = NULL,
  max_age = NULL,
  Linf = NULL,  # in cm
  KK = NULL,
  temp = NULL, # in degree C
  size_len = NULL, # in cm
  genus = NULL,
  species = NULL
  
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
  
  nat_mort <- NULL
  low_nat_mort <- NULL
  upp_nat_mort <- NULL
  FishLife_estimate_taxon_level <- NULL
  
  if(is.null(whose)){
  
    stop("name of M estimator was not specified")
  
  }else if(whose=="tanaka"){
    
    if(is.null(max_age)) stop("max_age was not given for tanaka")

    nat_mort <- 2.5/max_age

  }else if(whose=="pauly1_L"){
    
    if(is.null(Linf)) stop("Linf was not given for pauly1_L")
    if(is.null(KK)) stop("KK was not given for pauly1_L")
    if(is.null(temp)) stop("temp was not given for pauly1_L")
    
    nat_mort <- 0.9849 * Linf^(-0.279) * KK^0.6543 * temp^0.4634
    
  }else if(whose=="pauly_L_update"){
    
    if(is.null(Linf)) stop("Linf was not given for pauly_L_update")
    if(is.null(KK)) stop("KK was not given for pauly_L_update")
    
    nat_mort <- 4.118 * Linf^(-0.33) * KK^0.73
    
  }else if(whose=="jensen2"){
    
    if(is.null(KK)) stop("KK was not given for jensen2")
    
    nat_mort <- 1.5*KK
    
  }else if(whose=="hoenig"){
    
    if(is.null(max_age)) stop("max_age was not given for hoenig")
    
    nat_mort <- 4.3/max_age
    
  }else if(whose=="hoenig_update"){
    
    if(is.null(max_age)) stop("max_age was not given for hoenig_update")
    
    nat_mort <- 4.899 * max_age^(-0.916)
    
  }else if(whose=="gislason1"){
    
    if(is.null(size_len)) stop("size_len was not given for gislason1")
    if(is.null(Linf)) stop("Linf was not given for gislason1")
    if(is.null(KK)) stop("KK was not given for gislason1")
    
    nat_mort <- 1.73 * size_len^(-1.61) *  Linf^(1.44) * KK
    
  }else if(whose=="gislason2"){
    
    if(is.null(size_len)) stop("size_len was not given for gislason2")
    if(is.null(Linf)) stop("Linf was not given for gislason2")
    if(is.null(KK)) stop("KK was not given for gislason2")
    
    nat_mort <- KK * (size_len/Linf)^(-1.5)
    
  }else if(whose=="thorson"){
    
    if( !is.null(genus) & !is.null(species) ){
      
      Predict.org <- try(
          FishLife::Plot_taxa( 
            FishLife::Search_species(Genus = genus, Species = species)$match_taxonomy,
            mfrow=c(2,2)
          )
      )
      
      if(class(Predict.org) != "try-error"){ # check whether M estimate at spp level is available
        
        nat_mort <- exp(as.numeric(Predict.org[[1]]$Mean_pred["M"]))
        
        log_M <- as.numeric(Predict.org[[1]]$Mean_pred["M"])
        log_SD <- as.numeric(sqrt(diag(Predict.org[[1]]$Cov_pred))["M"])
        
        # mean M +- 1SD
        low_nat_mort <- exp(log_M-log_SD)
        upp_nat_mort <- exp(log_M+log_SD)
        
        FishLife_estimate_taxon_level <- "Species"
        
      }else{
        
        warning("M estimate at species level was not available, only at genus level for thorson")
        
        Predict.org <- try(
          FishLife::Plot_taxa( 
            FishLife::Search_species(Genus = genus)$match_taxonomy,
            mfrow=c(2,2)
          )
        )
        
        if(class(Predict.org) != "try-error"){ # check whether M estimate at genus level is available
          
          nat_mort <- exp(as.numeric(Predict.org[[1]]$Mean_pred["M"]))
          
          log_M <- as.numeric(Predict.org[[1]]$Mean_pred["M"])
          log_SD <- as.numeric(sqrt(diag(Predict.org[[1]]$Cov_pred))["M"])
          
          # mean M +- 1SD
          low_nat_mort <- exp(log_M-log_SD)
          upp_nat_mort <- exp(log_M+log_SD)
          
          FishLife_estimate_taxon_level <- "Genus"
          
        }else{
          
          stop("M estimate at genus level was not available for thorson")
          
        }
        
      }
      
    }else{
      
      stop("genus and species names were not given for thorson")
      
    }
    
  }else{
    
    stop("name of M estimator was mis-specified")
    
  }
  
  M_estimate <- list(
    M = nat_mort,
    M_1SD = c(low_nat_mort, upp_nat_mort),
    FishLife_estimate_taxon_level = FishLife_estimate_taxon_level
  )
  
  return( M_estimate )
  
}

