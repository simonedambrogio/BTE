#' Belonging Threshold (BT) estimation- informative prior
#'
#' This function estimates theBelonging Threshold for each rater
#' in case of an informative prior
#'
#' @param nrt Number of raters.
#' @param ni Number of items.
#' @param na Number of symptoms/attributes investigated by items.
#' @param data Dataframe containing the original attribution.
#'     It can be a `as_BCM` object.
#' @param NC Number of chains. Default is equal to four.
#' @param NI Number of iteration per chain. Default is equal to 100
#' @param data_prior A dataframe containing data representing the informative prior
#' @return TO BE DEFINED...
#' @export
#' @seealso `as_BCM`, `toyEx`
#' @examples
#' data(toyEx)
#'
#' ## Estimate BTs for three rater with a weakly informative prior
#' expertBMM<-expertBMM[,2:9]
#' BT_IP(3,12,8,data=toyEx,data_prior=expertBMM)
BT_IP<-function(nrt,ni,na,data,NC=4,NI=100,data_prior){
  N<-ni*na
  nr=2
  df<-create_BCM_list(nrt,ni,na,data)
  mybiglist_IP<-list()
  combinations<-combn(1:nrt,2)
  for (i in 1:ncol(combinations)) {
    couple<-combinations[,i]
    BCM.r.1 <- c(unlist(df[[couple[1]]]))
    BCM.r.2 <- c(unlist(df[[couple[2]]]))
    BCM<- c(BCM.r.1,BCM.r.2)

    informative.prior <- c(as.matrix(data_prior)) # Mean Raters' BMM prior from US_BMM
    sigma <- rep(.05,N*nr)      # Standard Deviation Raters' BMM prior

    dataList <- list(N=N, nr=nr, BCM=BCM, PUS_BMM=informative.prior, sigma=sigma) #  data for Stan model

    BT.model.i <- stan(file="inst\\stan\\BT_estimation.stan",
                       data=dataList, chains=NC, iter=NI) # run Stan model


    BT.post.i <- rstan::extract(BT.model.i,"BT")$BT      # BT posterior distribution
    BT.e.i <- round(apply(BT.post.i, 2, mean),3)  # BT estimated value


    if(BT.e.i[1]<.11) BT.e.i [1] <- .11 # check
    if(BT.e.i[2]<.11) BT.e.i [2] <- .11 # check

    BMM.e.i <- rstan::extract(BT.model.i,"R_BMM")$R_BMM # BMM posterior distribution

    post.pred.i <- c(ifelse(apply(BMM.e.i[,    1:N    ], 2, mean) > BT.e.i[1],1,0),
                     ifelse(apply(BMM.e.i[,(N+1):(N*2)], 2, mean) > BT.e.i[2],1,0)) # posterior predictive

    badfitted.i <- xor(post.pred.i, BCM) # bifference between posterior predictive and observed data


    name <- paste(as.character(couple),collapse="-")
    tmp <- list(informative.prior,
                BT.e.i, badfitted.i)
    mybiglist_IP[[name]] <- tmp

  }


  return(mybiglist_IP)

}
