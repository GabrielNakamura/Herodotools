#' Auxiliary function used to define the occurrence threshold used to define node presences
#'
#' @param x A matri containing the continuous probabilities of occurrences of 
#'     nodes accross the sites 
#' @param threshold.steps A numeric vector defining the threshold values that will
#'     be compared in the procedure
#'
#' @returns A list object containing the following information:
#' 
#'     \itemize{
#'         \item Maximum_site_probability: a vector containing the maximum probability
#'             of occurrence of nodes for ech site
#'         \item Removed_sites
#'         \item Correlation_between_steps
#'         \item Select_threshold: a scalar indicating the selected threshold
#'     }
#'     
#' @export
#'
#' @examples
find_threshold <- function(x, 
                           threshold.steps = c(0.4, 0.45, 0.5,
                                             0.55, 0.6, 0.65, 0.7, 
                                             0.75, 0.8, 0.85, 0.9)){
  max<-numeric()
  for (m in 1:nrow(x)){
    max[m]<-max(x[m,])
  }
  
  sites_off<-which(max<=min(threshold.steps))
  if (length(sites_off)>0){
    x=x[-sites_off,]
  } else { x=x
  }
  steps_off<-which(threshold.steps>max(max))
  if (length(steps_off)>0){
    threshold.steps=threshold.steps[-steps_off]
  } else { threshold.steps=threshold.steps
  }
  
  dist_threshold<-matrix(NA,nrow=(nrow(x)*(nrow(x)-1))/2,ncol=length(threshold.steps))
  colnames(dist_threshold)<-threshold.steps
  
  for (i in 1:length(threshold.steps)){
    vec<-ifelse(x>=threshold.steps[i],1,0)
    #dist_threshold[,i]<-as.vector(vegan::vegdist(vec,"bray"))
    dist_threshold[,i]<-betapart::beta.pair(vec,"jaccard")$beta.jtu
  }
  dist_threshold<-ifelse(is.nan(dist_threshold)==TRUE,1,dist_threshold)
  cor_steps<-cor(dist_threshold)
  names_obj<-c("full",threshold.steps[-c(length(threshold.steps)-1,
                                         length(threshold.steps))])
  eig_1<-eig_others<-eig_diff<-matrix(NA,nrow=length(names_obj), ncol=1,
                                      dimnames=list(names_obj,"value"))
  dist_threshold_test<-dist_threshold
  for (j in 1:length(names_obj)){
    pc_dist<-ade4::dudi.pca(dist_threshold_test,scannf = F,nf=2)
    eig_1[j,]<-pc_dist$eig[1]/sum(pc_dist$eig)
    eig_others[j,]<-sum(pc_dist$eig[2:length(pc_dist$eig)])/sum(pc_dist$eig)
    eig_diff[j,]<-eig_1[j,]-eig_others[j,]
    ade4::s.corcircle(pc_dist$co)
    dist_threshold_test<-dist_threshold_test[,-1]
  }
  eigen<-cbind(eig_1,eig_others)
  colnames(eigen)<-c("First_eigenvalue","All_other_eingenvalues")
  set_diff<-matrix(NA,nrow=length(names_obj)-1, ncol=1,
                   dimnames=list(names_obj[-1],"value"))
  for (k in 1:(length(eig_diff)-1)){
    set_diff[k,]<-eig_diff[k+1]-eig_diff[k]
  }
  #plot(x=rownames(set_diff),y=set_diff,type="l",ylab="Eigenvalue increment",xlab="Excluded threshold")
  d_crit<-hclust(dist(set_diff),"ward.D2")
  plot(d_crit,xlab="Excluded threshold")
  sel=set_diff
  return(list(Maximum_site_probability=max,Removed_sites=sites_off,Correlation_between_steps=cor_steps,Correlation_between_steps=eigen,Select_threshold=sel))
}