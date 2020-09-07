#' Belonging Measure Matrix (BMM) estimation-weakly informative prior
#'
#' This function estimates the Belonging Measure Matrix, starting from
#' the Belonging Threshold in case of a weakly informative prior
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
#' BMM_WIP(3,12,8,data=toyEx)
BMM_WIP<-function(nrt,ni,na,data,NC=4,NI=100){
  N<-ni*na
  nr=2
  df<-create_BCM_list(nrt,ni,na,data)
  sv<-BT_WIP(3,12,8,data) #BT values
  myBMM_WIP<-list()
  combinations<-combn(1:nrt,2)
  for (i in 1:ncol(combinations)) {
    couple<-combinations[,i]
    BCM.r.1 <- c(unlist(df[[couple[1]]]))
    BCM.r.2 <- c(unlist(df[[couple[2]]]))
    BCM<- c(BCM.r.1,BCM.r.2)
    #a[[3]][[2]][1]
    sv_i<-sv[[paste(as.character(couple),collapse="-")]]

    init_values <- c(ifelse(BCM.r.1==0,sv_i[[2]][1]/2, (sv_i[[2]][1] + (1-sv_i[[2]][1])/2)),
                     ifelse(BCM.r.2==0,sv_i[[2]][2]/2, (sv_i[[2]][2] + (1-sv_i[[2]][2])/2)))


    # init_values <- c(ifelse(BCM.r.1==0,BT.e.wi[1]/2, (BT.e.wi[1] + (1-BT.e.wi[1])/2)),
    #                  ifelse(BCM.r.2==0,BT.e.wi[2]/2, (BT.e.wi[2] + (1-BT.e.wi[2])/2)))
    #
    init_list <- list(c1=list(theta_r=init_values),
                      c2=list(theta_r=init_values),
                      c3=list(theta_r=init_values),
                      c4=list(theta_r=init_values))


    sigma <- rep(.2,N*nr)  # Standard Deviation Raters' BMM prior
    sigma[sv_i[[3]]] <- 1
    dim(sv_i[[2]]) <- nr ### VERIFICARE QUI


    dataList <- list(N=N, nr=nr, BCM=BCM, PUS_BMM= sv_i[[1]], BT=sv_i[[2]], sigma=sigma) # data for Stan model


    BMM.model.wi <- rstan::sampling(stanmodels$BMM_estimation.stan, data=dataList, chains=NC, iter=NI, init = init_list) # run Stan model


    BMM.post.wi <- rstan::extract(BMM.model.wi,"R_BMM")$R_BMM # BMM posterior distribution
    BMM.e.wi    <- apply(BMM.post.wi,2,mean)           # BMM estimated value
    BCM.e.wi.5  <- ifelse(BMM.e.wi>.5,1,0)             # BCM estimated from BMM and standard BT (.5)



    name <- paste(as.character(couple),collapse="-")
    tmp <- list(sv_i[[1]],sv_i[[2]],sv_i[[3]],
                BCM.e.wi.5)
    myBMM_WIP[[name]] <- tmp

  }


  return(myBMM_WIP)

}
