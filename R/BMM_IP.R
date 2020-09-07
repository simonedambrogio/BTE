#' Belonging Measure Matrix (BMM) estimation- informative prior
#'
#' This function estimates the Belonging Measure Matrix, starting from
#' the Belonging Threshold in case of an informative prior
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
#' BMM_IP(3,12,8,data=toyEx,data_prior=expertBMM)
BMM_IP<-function(nrt,ni,na,data,NC=4,NI=100,data_prior){
  N<-ni*na
  nr=2
  df<-create_BCM_list(nrt,ni,na,data)
  sv<-BT_IP(nrt=nrt,ni=ni,na=na,data=data,data_prior=data_prior) #BT values
  myBMM_IP<-list()
  combinations<-combn(1:nrt,2)
  for (i in 1:ncol(combinations)) {
    couple<-combinations[,i]
    BCM.r.1 <- c(unlist(df[[couple[1]]]))
    BCM.r.2 <- c(unlist(df[[couple[2]]]))
    BCM<- c(BCM.r.1,BCM.r.2)
    #a[[3]][[2]][1]
    sv_ip<-sv[[paste(as.character(couple),collapse="-")]]

    init_values <- c(ifelse(BCM.r.1==0,sv_ip[[2]][1]/2, (sv_ip[[2]][1] + (1-sv_ip[[2]][1])/2)),
                     ifelse(BCM.r.2==0,sv_ip[[2]][2]/2, (sv_ip[[2]][2] + (1-sv_ip[[2]][2])/2)))


    # init_values <- c(ifelse(BCM.r.1==0,BT.e.wi[1]/2, (BT.e.wi[1] + (1-BT.e.wi[1])/2)),
    #                  ifelse(BCM.r.2==0,BT.e.wi[2]/2, (BT.e.wi[2] + (1-BT.e.wi[2])/2)))
    #
    init_list <- list(c1=list(theta_r=init_values),
                      c2=list(theta_r=init_values),
                      c3=list(theta_r=init_values),
                      c4=list(theta_r=init_values))


    sigma <- rep(.2,N*nr)  # Standard Deviation Raters' BMM prior
    sigma[sv_ip[[3]]] <- 1
    dim(sv_ip[[2]]) <- nr ### VERIFICARE QUI
    #print(dim(sv_ip[[2]]))


    dataList <- list(N=N, nr=nr, BCM=BCM, PUS_BMM= sv_ip[[1]], BT=sv_ip[[2]], sigma=sigma) # data for Stan model


    BMM.model.i <- rstan::sampling(stanmodels$BMM_estimation.stan, data=dataList, chains=NC, iter=NI, init = init_list) # run Stan model


    BMM.post.i <- rstan::extract(BMM.model.i,"R_BMM")$R_BMM # BMM posterior distribution
    BMM.e.i    <- apply(BMM.post.i,2,mean)           # BMM estimated value
    BCM.e.i.5  <- ifelse(BMM.e.i>.5,1,0)             # BCM estimated from BMM and standard BT (.5)



    name <- paste(as.character(couple),collapse="-")
    tmp <- list(sv_ip[[1]],sv_ip[[2]],sv_ip[[3]],
                BCM.e.i.5)
    myBMM_IP[[name]] <- tmp

  }


  return(myBMM_IP)

}
