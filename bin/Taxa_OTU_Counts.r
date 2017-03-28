#Make Table of Taxa and OTU counts per body site

#Add taxonomy to the full OTU table (samples with 50 counts, otus in >1 sample)
otutable_here <- otutable2
#otutable_here <- st_otu
taxonomy_ids_keep <- intersect(rownames(taxonomy), rownames(otutable_here))
otutable_here <- otutable_here[taxonomy_ids_keep,] 
taxonomy_all <- taxonomy[taxonomy_ids_keep,]

rownames(otutable_here) <- taxonomy_all$taxonomy7

#Get total number of taxa at each body site, and the average number of taxa
body_sets <- list(Anal_B, Skin, Oral, Anal_M, Vagina)
taxa_totals <- c()
mean_taxa <- c()
site_names <- c("Anal", "Skin", "Oral", "Anal_M", "Vagina")
all_taxa_per_site <- c()
for(i in 1:length(body_sets)){
  body_set <- body_sets[[i]]
  subset_taxa_table <- taxa_table[,colnames(taxa_table) %in% body_set]
  subset_taxa_table <- subset_taxa_table[rowSums(subset_taxa_table) > 0,]
  all_taxa_per_site[[i]] <- rownames(subset_taxa_table)
  num_taxa <- nrow(subset_taxa_table)
  m_taxa <- round(mean(colSums(subset_taxa_table > 0)), digits=2)
  sd_taxa <- round(sd(colSums(subset_taxa_table > 0)),digits=2)
  m_taxa <- paste(m_taxa, " +/- ", sd_taxa)
  mean_taxa <- c(mean_taxa, m_taxa)
  taxa_totals <- c(taxa_totals, num_taxa)
  names(taxa_totals)[i] <- site_names[i]
  names(mean_taxa)[i] <- site_names[i]
}
taxa_totals <- c("Total", taxa_totals)
names(taxa_totals)[1] <- "Taxa"
mean_taxa <- c("Mean +/- SD", mean_taxa)
names(mean_taxa)[1] <- "Taxa"
  
#Get top 10 bugs and subset otu
top_taxa <- c()
body_sets <- list(Oral, Anal_B, Skin, Anal_M, Vagina)

for(i in 1:length(body_sets)){
  body_set <- body_sets[[i]]
  subset_taxa_table <- taxa_table[,colnames(taxa_table) %in% body_set]
  subset_taxa <- rownames(subset_taxa_table[order(rowSums(subset_taxa_table), decreasing=TRUE)[1:10],])
  top_taxa <- c(top_taxa, subset_taxa)
}
top_taxa <- unique(top_taxa)

#Subset otus to body site and mom/baby, remove otus with no counts
otus_baby_anal <- otutable_here[rownames(otutable_here) %in% top_taxa,Anal_B]
otus_baby_anal <- otus_baby_anal[rowSums(otus_baby_anal) > 0,]
otus_baby_skin <- otutable_here[rownames(otutable_here) %in% top_taxa,Skin]
otus_baby_skin <- otus_baby_skin[rowSums(otus_baby_skin) > 0,]
otus_baby_oral <- otutable_here[rownames(otutable_here) %in% top_taxa,Oral]
otus_baby_oral <- otus_baby_oral[rowSums(otus_baby_oral) > 0,]
otus_mom_anal <- otutable_here[rownames(otutable_here) %in% top_taxa,Anal_M]
otus_mom_anal <- otus_mom_anal[rowSums(otus_mom_anal) > 0,]
otus_mom_vagina <- otutable_here[rownames(otutable_here) %in% top_taxa,Vagina]
otus_mom_vagina <- otus_mom_vagina[rowSums(otus_mom_vagina) > 0,]

#Make output table
taxa_otu_counts <- data.frame(top_taxa)
taxa_otu_counts <- cbind(taxa_otu_counts, c(1:nrow(taxa_otu_counts)),c(1:nrow(taxa_otu_counts)),c(1:nrow(taxa_otu_counts)),c(1:nrow(taxa_otu_counts)),c(1:nrow(taxa_otu_counts)))
colnames(taxa_otu_counts) <- c("Taxa","Anal", "Skin", "Oral", "Anal_M", "Vagina")
tables <- list(otus_baby_anal, otus_baby_skin, otus_baby_oral, otus_mom_anal, otus_mom_vagina)

for(k in 1:length(tables)){
  current_table <- tables[[k]]
  for(i in 1:length(top_taxa)){
    current_t <- top_taxa[i]
    t_count <- sum(rownames(current_table) == current_t)
    taxa_otu_counts[i,(k+1)] <- t_count
  }
}
taxa_otu_counts$Taxa <- as.character(taxa_otu_counts$Taxa)
taxa_otu_counts <- rbind(taxa_otu_counts, taxa_totals)
taxa_otu_counts <- rbind(taxa_otu_counts, mean_taxa)
taxa_otu_counts[nrow(taxa_otu_counts),1] <- "Total Taxa"

#Write file to output
table_name <- paste(main_fp, "Taxa_otu_counts.txt", sep="/")
write.table(taxa_otu_counts, table_name, sep="\t", quote=FALSE, col.names=NA)

#which taxa are in all sites (mom and baby)
core_taxa <- Reduce(intersect, all_taxa_per_site)

core_baby <- Reduce(intersect, list(all_taxa_per_site[[1]], all_taxa_per_site[[2]], all_taxa_per_site[[3]]))  

core_mom <- Reduce(intersect, list(all_taxa_per_site[[4]], all_taxa_per_site[[5]]))  
