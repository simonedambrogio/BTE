#' Belonging Threshold (BT) estimation-weakly informative prior
#'
#' This function estimates theBelonging Threshold for each rater
#' in case of a weakly informative prior
#'
#' @param nrt Number of raters.
#' @param ni Number of items.
#' @param na Number of symptoms/attributes investigated by items.
#' @param data Dataframe containing the original attribution.
#'     It can be a `as_BCM` object.
#' @param NC Number of chains. Default is equal to four.
#' @param NI Number of iteration per chain. Default is equal to 100
#' @return TO BE DEFINED...
#' @export
#' @seealso `as_BCM`, `toyEx`
#' @examples
#' data(toyEx)
#'
#' ## Estimate BTs for three rater with a weakly informative prior
#' BT_WIP(3,12,8,data=toyEx)
BT_WIP<-function(nrt,ni,na,data,NC=4,NI=100){
  N<-ni*na
  nr=2
  df<-create_BCM_list(nrt,ni,na,data)
  mybiglist_WIP<-list()
  combinations<-combn(1:nrt,2)
    for (i in 1:ncol(combinations)) {
      couple<-combinations[,i]
      BCM.r.1 <- c(unlist(df[[couple[1]]]))
      BCM.r.2 <- c(unlist(df[[couple[2]]]))
      BCM<- c(BCM.r.1,BCM.r.2)

      weakly.informative.prior <- ifelse(xor(BCM.r.1, BCM.r.2),.5,ifelse(BCM.r.1==1,.75,.25)) # Mean Raters' BMM prior
      sigma <- rep(.5, N*nr) # Standard Deviation Raters' BMM prior

      dataList <- list(N=N, nr=nr, BCM=BCM, PUS_BMM=weakly.informative.prior, sigma=sigma) # data for Stan model

      BT.model.wi <- stan(file="stan\\BT_estimation.stan",
                          data=dataList, chains=NC, iter=NI) # run Stan model

      BT.post.wi <- rstan::extract(BT.model.wi,"BT")$BT   # BT posterior distribution
      BT.e.wi <- round(apply(BT.post.wi, 2, mean),3)  # BT estimated value


      BMM.e.wi <- rstan::extract(BT.model.wi,"R_BMM")$R_BMM # BMM posterior distribution

      post.pred.wi <- c(ifelse(apply(BMM.e.wi[,    1:N    ], 2, mean) > BMM.e.wi[1],1,0),
                        ifelse(apply(BMM.e.wi[,(N+1):(N*2)], 2, mean) > BMM.e.wi[2],1,0)) # posterior predictive

      badfitted.wi <- xor(post.pred.wi, BCM) # difference between posterior predictive and observed data

      name <- paste(as.character(couple),collapse="-")
      tmp <- list(weakly.informative.prior,
                  BT.e.wi, badfitted.wi)
      mybiglist_WIP[[name]] <- tmp

    }


  return(mybiglist_WIP)

}
