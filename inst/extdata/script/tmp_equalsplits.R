internal.brlen_div * switch(type, equal.splits = sort(rep(0.5,
                                                          length(internal.brlen_div))^c(1:length(internal.brlen_div))),
                            fair.proportion = {
                              for (j in 1:length(nodes_div)) {
                                # j = 1
                                sons <- picante::.node.desc(tree, nodes_div[j])
                                n.descendents <- length(sons$tips)
                                if (j == 1) portion <- n.descendents else portion <- c(n.descendents,
                                                                                       portion)
                              }
                              1/portion
                            })


# equal splits metric

tree$tip.label

equal.splits <- function(tree) {
  tips <- tree$tip.label
  # tip = 11
  sapply(tips, function(tip) {
    path <- nodepath(tree, tip)[[11]]
    bls <- tree$edge.length[match(path[-1], tree$edge[,2])]
    sum(bls / 2 ^ (seq_along(bls) - 1))
  })
}

path <- nodepath(tree, tip)
bls <- tree$edge.length[match(path[-1], tree$edge[,2])]
sum(bls / 2 ^ (seq_along(bls) - 1))


internal.brlen_div * sort(rep(0.5,
         length(internal.brlen_div))^c(1:length(internal.brlen_div)))
