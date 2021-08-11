## The simple t-test loop was modified based on the scrpit from the URL below,
## https://stackoverflow.com/questions/37474672/looping-through-t-tests-for-data-frame-subsets-in-r/37479506#37479506

# input clinical data
sample.list = read.csv("....csv", header = TRUE, rowname = 1)
# input matrix, the column order should be the same as the inputted clinical data
comps = read.csv("matrix.csv", header = TRUE, rowname = 1)

# make groupwise comparisons for all groups
comps = expand.grid(unique(sample.list$Group)[-5], ## excluding Control group [-5], it can be adjusted based on the need
                    names(sample.list)[2:ncol(sample.list)]) 
head(comps)

comps$pval = apply(comps, 1, function(x) {
  t.test(sample.list[sample.list$Group=="Ctrl", x[2]], 
         sample.list[sample.list $Group==x[1], x[2]])$p.value}
                  )

# call reshape2 package to re-format the result
library(reshape2)
comps2 <- dcast(comps, Var2~Var1, value.var = "pval")

#output result
write.csv(comp2, "p_values.csv")
