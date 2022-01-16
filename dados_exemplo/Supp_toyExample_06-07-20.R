####Supplementary material - Toy example########

toy_treeEx<- ape::read.tree(here::here("data", "toy_tree.new"))
toy_treeEx<- ape::makeNodeLabel(toy_treeEx, method= "user", nodeList= list(N6ABC= c("s"), 
                                                                           N7AB= c("s1", "s2", "s3"),
                                                                           N9C= c("s4", "s5"),
                                                                           N8A= c("s1", "s2")
)
)
quartz()
plot(toy_treeEx,  show.node.label = T)
ape::axisPhylo()

### current occurrence area ####
W_toy<- matrix(c(0, 1, 1,
             0, 1, 1,
             0, 1, 0,
             1, 0, 0,
             1, 0, 0
             ), nrow= 3, ncol= 5, dimnames=list(c("Comm 1", "Comm 2", "Comm 3"),
                                                        c(paste("s", 1:5, sep=""))
                                                )
           )
W_toy<- matrix(c(1, 1, 0, 0, 0,
               1, 1, 0, 0, 0,
               0, 0, 1, 1, 1,
               0, 0, 0, 1, 1, 
               0, 0, 0, 1, 1), nrow= 5, ncol= 5,
               dimnames=list(c("Comm 1", "Comm 2", "Comm 3", "Comm 4", "Comm 5"),
                             c(paste("s", 1:5, sep=""))
                             ), 
               byrow= TRUE)

biogeo_toy<- data.frame(Ecoregion= c("A", "A", "B", "C", "D"))
ancestral_area_toy<- data.frame(state= c("ABC", "AB", "A", "C"))

####calculating age arrival with toy example#####
age_arrival_toy<- diversification.assembly(W = W_toy, tree = toy_treeEx, ancestral.area = ancestral_area_toy, biogeo = biogeo_toy)$age_arrival
apply(age_arrival_toy, MARGIN = 1, function(x) mean(x[which(x != 0)]))
abs(node.depth.edgelength(phy = toy_treeEx)
    -max(node.depth.edgelength(toy_treeEx)))[-c(1:length(toy_treeEx$tip.label))]
toy_treeEx$node.label
which(age_arrival_toy[1, ] != 0)

divB_test <- DivB_metrics(W = W_toy, 
                          tree = toy_treeEx, 
                          ancestral.area = ancestral_area_toy, 
                          biogeo = biogeo_toy, 
                          diversification = c("jetz"),
                          dispersal.from = T,
                          age.arrival = T)
