matrix("-",
       nrow = 10,
       ncol = 10,
       )


age_arrival<-  matrix(0,
                      nrow = nrow(W),
                      ncol = ncol(W),
                      dimnames = list(rownames(W), colnames(W)))

