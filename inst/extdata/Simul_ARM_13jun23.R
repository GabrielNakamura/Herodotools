simul.comm<-function(Ncomm,Nspp,envir=FALSE,analyze.E=TRUE,replace.E=FALSE,diag=FALSE,u=5,power=5,binary=FALSE,runs=30,test=FALSE,nperm=999,parallel=8){
# Compute matrix P:
    matrix.p <- function(comm, phylodist){
        matrix.w <- as.matrix(comm)
        phylodist <- as.matrix(phylodist)
        similar.phy <- 1 - (phylodist/max(phylodist))
        matrix.phy <- 1/colSums(similar.phy)
        matrix.q <- sweep(similar.phy, 1, matrix.phy, "*")
          if(diag==FALSE){
            diag(matrix.q)<-0
          } else { matrix.q=matrix.q
            }
        matrix.P <- matrix.w %*% matrix.q
        return(list(matrix.w = matrix.w, matrix.q = matrix.q, matrix.P = matrix.P))
      }
    
# Double center matrix transformation: 
    matrix.double.center <- function(mat){
      mean_row_partial <- apply(mat, MARGIN = 2, 
                                function(x){
                                  x - rowMeans(mat)
                                }
      )
      mean_col_partial <- t(apply(mean_row_partial, MARGIN = 1, 
                                  function(x){
                                    x - colMeans(mat)
                                  }))
      double_center_matrix <- mean_col_partial + mean(as.matrix(mat))
      return(double_center_matrix)
    }

  #Res<-matrix(NA,nrow=runs,ncol=11,dimnames=list(1:runs,paste(c("R2.W","Coef.a.W","Coef.b.W","Beta.W","p.value","Value.PCPS.1","Value.PCPS.2","R2.Res","Coef.a.Res","Coef.b.Res","Beta.Res"))))
  if (analyze.E==TRUE){
    Res<-matrix(NA,nrow=runs,ncol=8,dimnames=list(1:runs,paste(c("Adj.R2","Coef.b.P","Beta.P","p.value.P","Coef.b.E","Beta.E","p.value.E","R2.Res"))))
  } else {Res<-matrix(NA,nrow=runs,ncol=5,dimnames=list(1:runs,paste(c("Adj.R2","Coef.b.P","Beta.P","p.value.P","R2.Res"))))
    }
      
  tree.list<-list()
  N<-matrix(NA,nrow=Nspp,ncol=runs,dimnames=list(1:Nspp,1:runs))
  K.N<-matrix(NA,nrow=runs,ncol=1,dimnames=list(1:runs,"K"))
  E<-matrix(NA,nrow=Ncomm,ncol=runs,dimnames=list(1:Ncomm,1:runs))
  EL<-matrix(NA,nrow=Ncomm,ncol=Nspp,dimnames=list(1:Ncomm,1:Nspp))
  L.list<-P.list<-EL.list<-pred.L.list<-resid.L.list<-list()
    for(k in 1:runs){
      phylo<-geiger::sim.bdtree(b=0.1,d=0,stop="taxa",n=Nspp,extinct=FALSE)
      tree.list[[k]]<-phylo
      a<-rnorm(Nspp,u,10)
      h<-runif(Nspp,0,30)
      L<-matrix(NA,Ncomm,Nspp,dimnames=list(paste("C", 1:Ncomm, sep=""),phylo$tip.label))
      rownames(N)<-phylo$tip.label
 	 
# Simulate niche N:    
      niche<-ape::rTraitCont(ape::compute.brlen(phylo,power=power),model="BM")
      N[,k]<-scales::rescale(niche,c(-1,101))
 	    K.N[k,]<-picante::Kcalc(N[,k],phylo)
# Simulate environment E:
 	    	 if(envir==TRUE){
	          E[,k]<-runif(Ncomm,0,100)
	       } else {E[,k]=N[,k]
	        }
# Simulate L:	    
		    for(i in 1:Ncomm){
			    for (j in 1:Nspp){
				    L[i,j]<-rpois(1,h[j]*exp((-((E[,k][i]-N[,k][j])^2))/(2*(a[j]^2))))
			    }
		    }
	    colnames(L)<-phylo$tip.label
	    rownames(L)<-rownames(as.matrix(L),FALSE,prefix="comm")
	      if(binary==TRUE){
	        L<-vegan::decostand(L,"pa")
	      } else{L=L}
	    L.list[[k]]<-L
	    L.cent<-matrix.double.center(L)

# Compute PCPS (silenced):
	    #pcps.L<-PCPS::pcps(L, cophenetic(phylo),checkdata = FALSE)
	    #Res[k,6]<-pcps.L$values[1,2]
	    #Res[k,7]<-pcps.L$values[2,2]

# Compute matrix P:	    
	    P<-matrix.p(comm=L,phylodist=cophenetic(phylo))$matrix.P
	    P.list[[k]]<-P
      P.cent<-matrix.double.center(P)
      
# run ARM:
      if(analyze.E==TRUE){
        mod.LP<-lm(as.numeric(L.cent)~as.numeric(P.cent))
        pred.LP<-matrix(predict(mod.LP),nrow(L),ncol(L),byrow=TRUE,dimnames=list(rownames(L),colnames(L)))
        resid.LP<-L.cent-pred.LP
          if(replace.E==TRUE){
            E[,k]<-runif(Ncomm,0,100)
          } else {E[,k]=E[,k]
            }
        
          for(i in 1:ncol(L)){
            EL[,i]<-as.matrix(E[,k])*resid.LP[,i]
          }
        EL.cent<-matrix.double.center(EL)
        colnames(EL.cent)<-colnames(L)
        rownames(EL.cent)<-rownames(L)
        EL.list[[k]]<-EL.cent	   
        
        mod.L<-lm(as.numeric(L.cent)~as.numeric(P.cent)+as.numeric(EL.cent))
        Res[k,1]<-summary(mod.L)$adj.r.squared
	      Res[k,2]<-summary(mod.L)$coefficients[2]
	      Res[k,3]<-QuantPsyc::lm.beta(mod.L)[1]
	      Res[k,5]<-summary(mod.L)$coefficients[3]
	      Res[k,6]<-QuantPsyc::lm.beta(mod.L)[2]
      } else {
        mod.L<-lm(as.numeric(L.cent)~as.numeric(P.cent))
        Res[k,1]<-summary(mod.L)$adj.r.squared
        Res[k,2]<-summary(mod.L)$coefficients[2]
        Res[k,3]<-QuantPsyc::lm.beta(mod.L)
        }
# Compute predicted and residual W matrices:
      pred.L<-matrix(predict(mod.L),nrow(L),ncol(L),byrow=TRUE,dimnames=list(rownames(L),colnames(L)))
      pred.L.list[[k]]<-pred.L
      resid.L<-L.cent-pred.L
      resid.L.list[[k]]<-resid.L
      
# power/type I error test####     
	      if(test==TRUE){
	        seqpermutation.E<-SYNCSA::permut.vector(nrow(as.matrix(E[,k])), nset = nperm)
	        seqpermutation.E <- lapply(seq_len(nrow(seqpermutation.E)), function(i) seqpermutation.E[i,])
	        seqpermutation.taxa <- SYNCSA::permut.vector(ncol(L), nset = nperm)
	        seqpermutation.taxa <- lapply(seq_len(nrow(seqpermutation.taxa)), function(i) seqpermutation.taxa[i,])
	        phylodist<- cophenetic(phylo)
	        p.n.taxa <- function(samp, comm, phylodist){
	          MP.null <- matrix.p(comm, phylodist[samp, samp])$matrix.P
	          return(MP.null)
	        }
	        
	        P.null <- lapply(seqpermutation.taxa, p.n.taxa, comm = L, phylodist = phylodist)
	        
	        P.null.cent_list<- lapply(P.null, FUN=function(x){
	          P.null.cent<-matrix.double.center(x)
	        })
	        
	        el.n.taxa<-function(samp, comm, E){
	          EL.null<-EL
	            for(i in 1:ncol(comm)){
	            EL.null[,i]<-as.matrix(E[samp,k])*L[,i]
	          }
	          return(EL.null)
	        }
	         
	        EL.null <- lapply(seqpermutation.E, el.n.taxa, comm = L, E=E)
	        EL.null.cent_list<- lapply(EL.null, FUN=function(x){
	          EL.null.cent<-matrix.double.center(x)
	        })

# start parallel:
	        newClusters <- FALSE
	        if (is.numeric(parallel)) {
	          parallel <- parallel::makeCluster(parallel, type = "PSOCK")
	          newClusters <- TRUE
	        }
	        if (!inherits(parallel, "cluster")) {
	          if(analyze.E==TRUE){
	            mod.L.null_list<- lapply(1:nperm, function(x){lm(as.numeric(L.cent) ~ as.numeric(P.null.cent_list[[x]]) + as.numeric(EL.null.cent_list[[x]]))})
	          } else {
	            mod.L.null_list<- lapply(P.null.cent_list, function(x){lm(as.numeric(L.cent) ~ as.numeric(x))})
	            }
	          
	        } else {
	          if(analyze.E==TRUE){
	             mod.L.null_list<- parallel::parLapply(cl = parallel,
	                                                X = 1:nperm, 
	                                                fun = function(i){lm(as.numeric(L.cent) ~ as.numeric(P.null.cent_list[[i]]) + as.numeric(EL.null.cent_list[[i]]))})
	          } else {mod.L.null_list<- parallel::parLapply(cl = parallel,
	                                                        X = P.null.cent_list, 
	                                                        fun = function(i){
	                                                          lm(as.numeric(L.cent) ~ as.numeric(i))
	                                                              }
	                                    )
	            }
	          }
	        
	        if (newClusters){
	          parallel::stopCluster(parallel)
	        }
	        
	     if(analyze.E==TRUE){
	            mod.beta.null <- matrix(unlist(lapply(mod.L.null_list, function(x){
	            QuantPsyc::lm.beta(x)})),
	            nrow = nperm, ncol= 2)
	       	    Res[k,4]<- (sum(ifelse(mod.beta.null[,1] > Res[k,3], 1, 0)) +1)/(nperm + 1)
	            Res[k,7]<- (sum(ifelse(mod.beta.null[,2] > Res[k,6], 1, 0)) +1)/(nperm + 1)
	     
	     } else {
	       mod.beta.null <- matrix(unlist(lapply(mod.L.null_list, function(x){
	       QuantPsyc::lm.beta(x)})),
	       nrow = nperm, ncol= 1)
	       Res[k,4]<- (sum(ifelse(mod.beta.null > Res[k,3], 1, 0)) +1)/(nperm + 1)
	     }
	  } else {p.value=NA
	   }

# Test residuals
	    resid.L.cent<-matrix.double.center(resid.L)
	      if(analyze.E==TRUE){
	          mod.resid.L<-lm(as.numeric(resid.L.cent)~as.numeric(P.cent)+as.numeric(EL.cent))
	           Res[k,8]<-summary(mod.resid.L)$adj.r.squared
	      } else {
	          mod.resid.L<-lm(as.numeric(resid.L.cent)~as.numeric(P.cent))
	          Res[k,5]<-summary(mod.resid.L)$adj.r.squared
	        }
      print(k)
    }
  return(list(Trees=tree.list,Niche=N,K.niche=K.N,Environment=E,EL.matrices=EL.list,
              W.matrices=L.list,P.matrices=P.list,Predicted.W.matrices=pred.L.list,
              Residual.W.matrices=resid.L.list,Model.Results=Res))
  }