#' Cohen's kappa between BCMs
#'
#' Thi function calculates the Cohen's kappa index
#' between n-couples of Boolean Classification Matrices
#'
#' @param nrt Number of raters.
#' @param ni Number of items.
#' @param na Number of symptoms/attributes investigated by items.
#' @param data Dataframe containing the original attribution.
#'     It can be a `as_BCM` object.
#' @return A dataframe containing the couples compared (first column),
#'     and the kappa value (second column).
#' @export
#' @seealso `as_BCM`, `toyEx`
#' @examples
#' data(toyEx)
#'
#' ## Compute SPA only between rater 1 and 3
#' kappa_BCMs(3,12,8,data=toyEx,whichRaters=c(1,3))
kappa_BCMs<-function(nrt,ni,na,data,whichRaters){
  N<-ni*na
  data=create_BCM_list(nrt,ni,na,data)
  data_kappa<-data.frame(Couple=character(),kappa=numeric(),
                       stringsAsFactors = FALSE)
  if(is.list(data) & whichRaters=="ALL"){
    combinations<-combn(1:nrt,2)
    for (i in 1:ncol(combinations)) {
      couple<-combinations[,i]
      BCM.r.1 <- c(unlist(data[[couple[1]]]))
      BCM.r.2 <- c(unlist(data[[couple[2]]]))
      BCM<- c(BCM.r.1,BCM.r.2)
      data_kappa[i,1]<-as.character(paste0(couple,collapse = "-"))
      data_kappa[i,2]<- round(irr::kappa2(cbind(BCM.r.1,BCM.r.2))$value,2)
    }
  } else if(is.list(data) & whichRaters!="ALL"){
    combinations<-combn(whichRaters,2)
    for (i in 1:ncol(combinations)) {
      couple<-combinations[,i]
      BCM.r.1 <- c(unlist(data[[couple[1]]]))
      BCM.r.2 <- c(unlist(data[[couple[2]]]))
      BCM<-BCM <- c(BCM.r.1,BCM.r.2)
      data_kappa[i,1]<-as.character(paste0(couple,collapse = "-"))
      data_kappa[i,2]<- round(irr::kappa2(cbind(BCM.r.1,BCM.r.2))$value,2)
    }
  } else if(is.list(data)==FALSE & whichRaters=="ALL"){
    df<-create_BCM_list(nrt,ni,na,data)
    combinations<-combn(1:nrt,2)
    for (i in 1:ncol(combinations)) {
      couple<-combinations[,i]
      BCM.r.1 <- c(unlist(df[[couple[1]]]))
      BCM.r.2 <- c(unlist(df[[couple[2]]]))
      BCM<-BCM <- c(BCM.r.1,BCM.r.2)
      data_kappa[i,1]<-as.character(paste0(couple,collapse = "-"))
      data_kappa[i,2]<- round(irr::kappa2(cbind(BCM.r.1,BCM.r.2))$value,2)
    }
  } else {
    df<-create_BCM_list(nrt,ni,na,data)
    combinations<-combn(whichRaters,2)
    for (i in 1:ncol(combinations)) {
      couple<-combinations[,i]
      BCM.r.1 <- c(unlist(df[[couple[1]]]))
      BCM.r.2 <- c(unlist(df[[couple[2]]]))
      BCM<-BCM <- c(BCM.r.1,BCM.r.2)
      data_kappa[i,1]<-as.character(paste0(couple,collapse = "-"))
      data_kappa[i,2]<- round(irr::kappa2(cbind(BCM.r.1,BCM.r.2))$value,2)
    }

  }

  return(data_kappa)

}
