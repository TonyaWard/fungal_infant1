#Visual the correlation network
#Collapse taxa that are correlating together

#Get Correlations for OTU  table and taxa table
Cor_out2 <- cor(t(taxa_table))
Cor_out_OTU <- cor(t(otutable_RA))
# Cor_out_OTU <- cor(t(otutable)) #option 1 only

#visualize as a heatmap
direct <- paste(main_fp, "correlation_network/", sep='/')
fp <- paste(direct, "original_correlation_taxa.pdf", sep="")
pdf(file = fp)
image(Cor_out2)
dev.off()

direct <- paste(main_fp, "correlation_network/", sep='/')
fp <- paste(direct, "original_correlation_otus.pdf", sep="")
pdf(file = fp)
image(Cor_out_OTU)
dev.off()

######################################################################
#Collapse by correlation (TAXA)
#Select rows where the correlation is > 0.95 
Cor_out3 <- abs(Cor_out2) > 0.95
new_taxatable <- taxa_table

for(i in 1:nrow(Cor_out3)){
  for(j in i:nrow(Cor_out3)){
    if(i==j){
      next
    }
    if(is.na(Cor_out3[j,i])){
      next
    }
    if(Cor_out3[j,i]){
      new_taxatable[j,] <- new_taxatable[j,] + new_taxatable[i,]
      new_taxatable[i,] <- rep(0, ncol(new_taxatable))
    }
  }
}

#Collapse
taxa_table.new <- new_taxatable[rowSums(new_taxatable) > 0, ]
#sanity check
#sum(taxa_table.new) == sum(taxa_table)

#New correlation
cor_out4 <- cor(t(taxa_table.new))

#print new heatmap showing anything correlating above 0.90
fp <- paste(direct, "remaining_correlations_taxa.pdf", sep="")
pdf(file = fp)
image(cor_out4 > .9)
dev.off()

######################################################################
#Collapse by correlation (OTU)
#Select rows where the correlation is > 0.95 
Cor_OTU2 <- abs(Cor_out_OTU) > 0.95
new_table <- otutable_RA
# new_table <- otutable # option 1 only

for(i in 1:nrow(Cor_OTU2)){
  for(j in i:nrow(Cor_OTU2)){
    if(i==j){
      next
    }
    if(is.na(Cor_OTU2[j,i])){
      next
    }
    if(Cor_OTU2[j,i]){
      new_table[j,] <- new_table[j,] + new_table[i,]
      new_table[i,] <- rep(0, ncol(new_table))
    }
  }
}

#Collapse
table.new <- new_table[rowSums(new_table) > 0, ]

#New correlation
cor_OTU3 <- cor(t(table.new))

#print new heatmap showing anything correlating above 0.90
fp <- paste(direct, "remaining_correlations_otus.pdf", sep="")
pdf(file = fp)
image(cor_OTU3 > 0.9)
dev.off()

######################################################################
#Check the correlation of each taxon and seq depth
#Add column that is map depth
depth_table <- otutable[,rownames(mapping)]
mapping$depth <- colSums(depth_table)

#re-order mapping file
taxa_map <- mapping[colnames(taxa_table.new),]
 
#Store p-vals for taxon-depth correlation analysis
pvals <- c()
for(i in 1:nrow(taxa_table.new)){
 pvals <- c(pvals,cor.test(as.numeric(taxa_table.new[i,]),taxa_map$depth, method="spearman")$p.val)
}
pvals <- p.adjust(pvals, method="BH")

#print correlations to stats file
fp <- paste(direct, "correlation_stats.txt", sep="")
sink(fp)
cat("BH corrected p-vals (top ten)\n\n")
print(pvals[order(pvals)][1:10])
print(rownames(taxa_table.new)[order(pvals)[1:5]])
sink()
 
#Plot the top 5 correlations
for(i in 1:5){
fp <- paste(direct, rownames(taxa_table.new)[order(pvals)[i]], ".pdf", sep="")
  pdf(file = fp)
  plot(as.numeric(taxa_table.new[order(pvals)[i],]) ~ taxa_map$depth)
  abline(lm(as.numeric(taxa_table.new[order(pvals)[i],]) ~ taxa_map$depth))
  dev.off()
}

#replace old taxa_table with the new taxa table
taxa_table <- taxa_table.new
otutable_RA <- table.new
# otutable <- table.new #option 1 only


