

#' Ancestral Diversity Analysis for assemblage data
#'
#' @param x 
#' @param phy 
#' @param node.names 
#' @param type 
#' @param marginal 
#' @param sp.bin 
#' @param lik.threshold 
#' @param threshold 
#' @param compute.node.by.sites 
#' @param compute.fields 
#' @param allowSplit 
#' @param credMass 
#' @param compute.metrics 
#'
#' @return
#' @export
#'
#' @examples
ada <- 
  function(x, 
           phy, 
           node.names = FALSE, 
           type = c("discrete", "continuous"),
           marginal = FALSE,
           sp.bin = 15,
           lik.threshold = FALSE,
           threshold = 0.20, 
           compute.node.by.sites = FALSE, 
           compute.fields = F, 
           allowSplit = T, 
           credMass = 0.95, 
           compute.metrics = FALSE){
    # Enter and organize data:
    match <- picante::match.phylo.comm(phy, x)
    x <- match$comm
    if(node.names == FALSE){
      phy <- ape::makeNodeLabel(match$phy, prefix = NULL)
    } else {phy = phy
    }
    root.age <- max(cophenetic(phy))/2
    
    # Extract species by nodes matrix:
    spp.nodes<-adephylo::treePart(phy, result="dummy")
    colnames(spp.nodes)<-2:length(phy$node.label)
    Node1<-matrix(1,nrow(spp.nodes),1,dimnames=list(rownames(spp.nodes),1))
    spp.nodes<-cbind(Node1,spp.nodes)
    
    # Run Ancestral Area Reconstruction:
    node.list<-list()
    node.anc.area<-node.samp.mat<-matrix(NA,nrow=phy$Nnode,ncol=nrow(x),dimnames=list(phy$node.label,rownames(x)))
    
    if(type=="continuous"){
      for(i in 1:nrow(x)){
        node.list[[i]]<-phytools::fastAnc(phy,x[i,])
        node.anc.area[,i]<-node.list[[i]]
        print(i)
      }
      if(lik.threshold==FALSE){
        m.node.anc.area<-rowMeans(node.anc.area)
        sd.node.anc.area<-numeric()
        for(i in 1:nrow(node.anc.area)){
          sd.node.anc.area[i]<-sd(as.numeric(node.anc.area[i,]))
        }
        for(i in 1:nrow(node.anc.area)){
          for(p in 1:ncol(node.anc.area)){
            node.anc.area[i,p]<-pnorm(q=(node.anc.area[i,p]-m.node.anc.area[i])/sd.node.anc.area[i],
                                      mean=0,sd=1,lower.tail=TRUE)
          }
        }
      } else {
        for(i in 1:nrow(x)){
          for (j in 1:phy$Nnode){
            ifelse(node.list[[i]][j]>=threshold,1,0)
          }
        }
        node.anc.area[,i]<-node.list[[i]]
      }
      
    } else {
      node.list.single<-list()
      for(i in 1:nrow(x)){
        node.list.single[[i]]<-apply(ape::ace(t(x)[,i], phy,type=type)$lik.anc,1,which.max)
        node.anc.area[,i]<-ifelse(node.list.single[[i]]==1,0,1)
      }
    }
    
    # Compute a matrix of nodes by sites
    if(compute.node.by.sites==TRUE){
      node.names.mat<-rownames(node.samp.mat)
      for (p in 1:nrow(node.samp.mat)){
        for (i in 1:ncol(node.samp.mat)){
          comm.samp<- x[,which(x[i,]==1)]
          samp.nodes <- picante::prune.sample(samp=comm.samp, phylo=phy)$node.label
          node.samp.mat[p,i]<-ifelse(node.names.mat[p]%in%samp.nodes,1,0)
        }
        print(p)
      }
    } else {node.samp.mat<-"Nops..."}
    
    # Define joint occurrence matrix:
    joint.x<-as.matrix(cbind(x,t(node.anc.area)))
    
    # Compute  metrics:
    if(compute.metrics==TRUE){
      
      # Compute node ages:
      ages<-round(as.matrix(root.age-ape::node.depth.edgelength(phy)),5)
      ages[1:length(phy$tip.label),]<-0.00001
      colnames(ages)<-"age"
      spp.names.nodes<-c(phy$tip.label,colnames(spp.nodes))
      rownames(ages)<-spp.names.nodes
      
      age.anc.area<-matrix(0,nrow(joint.x),ncol(joint.x),
                           dimnames=list(rownames(joint.x),colnames(joint.x)))
      for (p in 1:nrow(joint.x)){
        age.anc.area[p,]<-ages*t(joint.x)[,p]
      }
      
      # Compute "LTT":
      h<-hist(ages,breaks=round(length(ages)/sp.bin,0),plot=F)
      breaks<-h$breaks[-1]
      n.breaks<-length(breaks)
      mid.time<-as.vector(h$mids)
      
      # Compute weights for species and internal nodes given "LTT":
      weights<-1/h$counts
      isna.weights<-which(is.na(ifelse(weights==Inf,NA,weights)))
      if(length(isna.weights)==0){
        weights<-weights
      } else{weights<-weights[-isna.weights]
      }
      spp.weight.class<-list()
      spp.weight.class[[1]]<-which(ages<=breaks[1])
      for(b in 1:n.breaks){
        spp.weight.class[[b+1]]<-which(ages>breaks[b]&ages<=breaks[b+1])
      }
      spp.weight.class<-spp.weight.class[-length(spp.weight.class)]
      spp.weight.class[lengths(spp.weight.class)==0]<-NA
      isna.spp.weight.class<-which(is.na(spp.weight.class))
      if(length(isna.spp.weight.class)==0){
        spp.weight.class.clean<-spp.weight.class
      } else{spp.weight.class.clean<-rlist::list.remove(spp.weight.class,
                                                        range = which(is.na(spp.weight.class)))
      }
      
      weight.mat<-matrix(NA,nrow=nrow(ages),ncol=length(weights),dimnames=list(rownames(ages),round(weights,3)))
      for(w in 1:nrow(weight.mat)){
        for(y in 1:length(weights)){
          weight.mat[w,y]<-round(ifelse(w%in%spp.weight.class.clean[[y]],yes=weights[y],no=0),3)
        }
      }
      weight.vec<-vector()
      for(w in 1:nrow(weight.mat)){
        weight.vec[w]<-sum(weight.mat[w,])
      }
      # Compute results:
      res.dens<-matrix(NA,nrow=nrow(age.anc.area),ncol=8,
                       dimnames=list(rownames(age.anc.area),c("N_species","N_ancestral_nodes",
                                                              "Min.Density.Peak.Myr","Max.Density.Peak.Myr","N.Diversity.Peaks",
                                                              "Assemblage.Age","Stability", "Mean.Diversification.Rate")))
      res.dens[,1]<-rowSums(x)
      res.dens[,2]<-rowSums(joint.x[,(length(phy$tip.label)+1):ncol(joint.x)])
      res.dens[,8]<-as.numeric(SYNCSA::matrix.t(x,
                                                traits=as.data.frame(1/picante::evol.distinct(phy,type="equal.splits",
                                                                                              scale = FALSE,use.branch.lengths = TRUE)[,2]), ranks=FALSE,
                                                notification=FALSE)$matrix.T)
      
      
      for(l in 1:nrow(age.anc.area)){
        nz.node.ages<-which(age.anc.area[l,]>0)
        weights.nz<-weight.vec[nz.node.ages]/sum(weight.vec[nz.node.ages])
        age.anc.area.nz<-age.anc.area[l,nz.node.ages]
        if(length(age.anc.area.nz)>=2){
          dens.area<-density(age.anc.area.nz,from=0,to=root.age,weights = weights.nz)
          res.dens[l,3]<-as.numeric(HDInterval::hdi(dens.area,allowSplit=F,credMass=0.05)[1])
          res.dens[l,4]<-as.numeric(HDInterval::hdi(dens.area,allowSplit=F,credMass=0.05)[2])                     
          hpd.res.dens<-HDInterval::hdi(dens.area,allowSplit=allowSplit,credMass=credMass)
          dim.hpd<-dim(hpd.res.dens)  
          occ<-vector() 
          if(dim(hpd.res.dens)[1]>1){
            for (z in 1:dim.hpd[1]){
              occ[z]<-hpd.res.dens[z,2]-hpd.res.dens[z,1]
            }
            occ.time<-sum(occ)
            res.dens[l,5]<-dim.hpd[1]
            res.dens[l,6]<-max(hpd.res.dens[,2])
            res.dens[l,7]<-as.numeric(occ.time/(hpd.res.dens[dim.hpd[1],dim.hpd[2]]-hpd.res.dens[1,1]))
          } else {res.dens[l,6]<-hpd.res.dens[2]
          res.dens[l,7]<-1
          }
        } else {res.dens[l,3]<-0
        res.dens[l,4]<-0
        res.dens[,5]<-1
        res.dens[l,6]<-0
        res.dens[l,7]<-1
        }
      }
      
      # Compute site number of nodes x time period:
      div.nodes<-matrix(NA,nrow(joint.x),ncol=length(spp.weight.class.clean),dimnames=list(rownames(joint.x),
                                                                                           c(mid.time[-which(is.na(spp.weight.class))])))
      for(d in 1:length(spp.weight.class.clean)){
        if(length(spp.weight.class.clean[[d]])>1){
          div.nodes[,d]<-rowSums(joint.x[,spp.weight.class.clean[[d]]])
        } else{div.nodes[,d]<-joint.x[,spp.weight.class.clean[[d]]]
        }
      }
      
      # Compute species fields:
      if(compute.fields==TRUE){
        hist.fields.div<-list()
        hist.fields.anc.nodes<-list()
        hist.fields.peak<-list()
        for(f in 1:ncol(x)){
          sites<-which(x[,f]>0)
          hist.fields.div[[f]]<-res.dens[sites,1]
          hist.fields.anc.nodes[[f]]<-res.dens[sites,2]
          hist.fields.peak[[f]]<-res.dens[sites,3]
        }
        Species.Fields<-list(Diversity.Fields=hist.fields.div,
                             Ancient.Node.Diversity.Fields=hist.fields.anc.nodes,
                             Diversity.Peak.Time.Fields=hist.fields.peak)
      } else {Species.Fields<-"Tell me to compute them, dude!!!"
      }   
    } else { res.dens<-"You've lost your mind..."
    Species.Fields<-"Tell me to compute them, dude!!!"
    div.nodes<-"WTF?!"
    ages<-"Not sure..."
    }   
    # Define outputs: 
    Res<-list(Phylogeny=phy,
              Root.Age=root.age,
              Species.per.node.matrix=spp.nodes,
              Node.ages=ages,
              Site.by.nodes=node.samp.mat,
              Per.node.ancestral.area=node.anc.area,
              Joint.species.nodes.matrix=joint.x,
              Diversity.Through.Time=div.nodes,
              Cell.Metrics=res.dens,
              Species.Fields,
              node.list=node.list)
    return(Res)
  }



