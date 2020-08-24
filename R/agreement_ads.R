#' Agreement's coefficients adjusted
#'
#' Thi function calculates different agreement's coefficients
#' between n-couples of raters
#'
#' @param nrt Number of raters.
#' @param ni Number of items.
#' @param na Number of symptoms/attributes investigated by items.
#' @param data Dataframe containing the original attribution.
#'     It can be a `as_BCM` object.
#' @return A dataframe containing the couples compared (first column),
#'     and the SPA value (second column).
#' @export
#' @seealso `as_BCM`, `toyEx`
#' @examples
#' data(toyEx)
#'
#' ## Compute SPA only between rater 1 and 3
#' agreement_ads(3,12,8,data=toyEx,prior="WI",index="Kappa")
agreement_ads<-function(nrt,ni,na,data,NC=4,NI=100,
                        data_prior=NULL,prior,index){
  N<-ni*na
  nr=2
  df<-create_BCM_list(nrt,ni,na,data)
  nof_row<-combn(1:nrt,2)
  namesCol<-c("Couple","BCM_R1","BCM_R2","SPA",
              "Kappa","BT_wi_R1","BT_wi_R2",
              "BT_i_R1","BT_i_R2","SPA_adj_wi",
              "Kappa_adj_wi","SPA_adj_i","Kappa_adj_i")
  Results<-matrix(NA,nrow = ncol(nof_row),ncol = length(namesCol))

  if(prior=="WI"){
    a<-BMM_WIP(nrt=nrt,ni=ni,na=na,data=data,NC=NC,NI=NI)
    combinations<-combn(1:nrt,2)
    for (i in 1:ncol(combinations)) {
      couple<-combinations[,i]
      BCM.r.1 <- c(unlist(df[[couple[1]]]))
      BCM.r.2 <- c(unlist(df[[couple[2]]]))
      BCM<- c(BCM.r.1,BCM.r.2)
      #a[[3]][[2]][1]
      a_i<-a[[paste(as.character(couple),collapse="-")]]

      BCM.r.1.e.wi.5 <- unlist(a_i[[4]])[    1:N    ] # Rater 1 - BCM estimated from BMM and standard BT (.5)
      BCM.r.2.e.wi.5 <- unlist(a_i[[4]])[(N+1):(N*2)] # Rater 2 - BCM estimated from BMM and standard BT (.5)

      if(index=="SPA"){
        SPA.wi.5  <- sum(BCM.r.1.e.wi.5==BCM.r.2.e.wi.5)/N # Estimated (wi) simple Percentage Agreement whit standard threshold (.5)
        Results[i,1]<-paste(as.character(couple),collapse="-")
        Results[i,2]<-mean(BCM.r.1)
        Results[i,3]<-mean(BCM.r.2)
        Results[i,4]<-sum(BCM.r.1==BCM.r.2)/N
        Results[i,5]<-NA
        Results[i,6]<-a_i[[2]][1]
        Results[i,7]<-a_i[[2]][2]
        Results[i,8]<-NA
        Results[i,9]<-NA
        Results[i,10]<-SPA.wi.5
        Results[i,11]<-NA
        Results[i,12]<-NA
        Results[i,13]<-NA

      } else {
        CK.wi.5<-irr:: kappa2(cbind(BCM.r.1.e.wi.5,BCM.r.2.e.wi.5))$value
        Results[i,1]<-paste(as.character(couple),collapse="-")
        Results[i,2]<-mean(BCM.r.1)
        Results[i,3]<-mean(BCM.r.2)
        Results[i,4]<-NA
        Results[i,5]<- irr::kappa2(cbind(BCM.r.1,BCM.r.2))$value
        Results[i,6]<-a_i[[2]][1]
        Results[i,7]<-a_i[[2]][2]
        Results[i,8]<-NA
        Results[i,9]<-NA
        Results[i,10]<-NA
        Results[i,11]<-CK.wi.5
        Results[i,12]<-NA
        Results[i,13]<-NA

      }

    }
  } else {
    a<-BMM_IP(nrt=nrt,ni=ni,na=na,data=data,NC=NC,NI=NI,data_prior = data_prior)
    combinations<-combn(1:nrt,2)
    for (i in 1:ncol(combinations)) {
      couple<-combinations[,i]
      BCM.r.1 <- c(unlist(df[[couple[1]]]))
      BCM.r.2 <- c(unlist(df[[couple[2]]]))
      BCM<- c(BCM.r.1,BCM.r.2)
      #a[[3]][[2]][1]
      a_i<-a[[paste(as.character(couple),collapse="-")]]

      BCM.r.1.e.i.5 <- unlist(a_i[[4]])[    1:N    ] # Rater 1 - BCM estimated from BMM and standard BT (.5)
      BCM.r.2.e.i.5 <- unlist(a_i[[4]])[(N+1):(N*2)] # Rater 2 - BCM estimated from BMM and standard BT (.5)

      if(index=="SPA"){
        SPA.i.5  <- sum(BCM.r.1.e.i.5==BCM.r.2.e.i.5)/N # Estimated (wi) simple Percentage Agreement whit standard threshold (.5)
        Results[i,1]<-paste(as.character(couple),collapse="-")
        Results[i,2]<-mean(BCM.r.1)
        Results[i,3]<-mean(BCM.r.2)
        Results[i,4]<-sum(BCM.r.1==BCM.r.2)/N
        Results[i,5]<-NA
        Results[i,6]<-NA
        Results[i,7]<-NA
        Results[i,8]<-a_i[[2]][1]
        Results[i,9]<-a_i[[2]][2]
        Results[i,10]<-NA
        Results[i,11]<-NA
        Results[i,12]<-SPA.i.5
        Results[i,13]<-NA

      } else {
        CK.i.5<-irr:: kappa2(cbind(BCM.r.1.e.i.5,BCM.r.2.e.i.5))$value
        Results[i,1]<-paste(as.character(couple),collapse="-")
        Results[i,2]<-mean(BCM.r.1)
        Results[i,3]<-mean(BCM.r.2)
        Results[i,4]<-NA
        Results[i,5]<- irr::kappa2(cbind(BCM.r.1,BCM.r.2))$value
        Results[i,6]<-NA
        Results[i,7]<-NA
        Results[i,8]<-a_i[[2]][1]
        Results[i,9]<-a_i[[2]][2]
        Results[i,10]<-NA
        Results[i,11]<-NA
        Results[i,12]<-NA
        Results[i,13]<-CK.i.5

      }

    }
  }
  Results<-data.frame(Results)
  colnames(Results)<-namesCol
  Results<-Filter(function(x) !all(is.na(x)),Results)
  return(Results)
}
