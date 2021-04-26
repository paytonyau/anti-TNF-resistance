# https://stackoverflow.com/questions/37474672/looping-through-t-tests-for-data-frame-subsets-in-r/37479506#37479506

comps <- expand.grid(unique(A$Group)[-5],  names(A)[2:ncol(A)]) ## excluding Control group

comps$pval <- apply(comps, 1, function(x) {
  t.test(A[A$Group=="Ctrl", x[2]], A[A$Group==x[1], x[2]])$p.value 
} )

library(reshape2)
comps2 <- dcast(comps, Var2~Var1, value.var = "pval")