#Find taxa increasing or decreasing over time

#set output dir
taxa_dir <- paste(main_fp, "taxa_time/", sep='/')

ALPHA <- 0.5
kids <- c(Bodysites_B[[1]], Bodysites_B[[2]], Bodysites_B[[3]])
test_table <- taxa_table[, kids]
test_table <- test_table[rowSums(test_table) > 0,] #take out 3 taxa
#test_table <- test_table[rowSums(test_table > 0 )/ncol(test_table) > 0.1, ] #keep only taxa in 10% of all infants

for(b in 1:length(Bodysites_B)){
  body_samples <- Bodysites_B[[b]]
  test_table2 <- test_table[,body_samples]
  test_table2 <- test_table2[rowSums(test_table2) > 0.001,]
  test_table2 <- as.data.frame(t(test_table2))
  test_map <- mapping[body_samples,]
  pvals <- c()
  for(t in 1:ncol(test_table2)){
    taxon <- rownames(test_table2)[t]
    test_out <- cor.test(test_table2[,t], test_map[,"planned_sampling_day_0_1_3_7_14_30_60_90___365"], method="spearman")
    pvals <- c(pvals, test_out$p.value)
  }
  adjusted <- p.adjust(pvals, method="fdr")
  for(p in 1:length(adjusted)){
    if(adjusted[p] < ALPHA){
      plot <- ggplot(test_table2, aes(x=test_map[,"planned_sampling_day_0_1_3_7_14_30_60_90___365"], y= test_table2[,p], color=test_map[,"SuperbodysiteOralSkinNoseVaginaAnalsAureola"])) +
        geom_point(alpha=0.65, size=4)+
        geom_smooth(method=lm, se=FALSE) +
        labs(y="Relative Abundance", x="Day") +
        scale_color_manual(values=body_cols[b]) +
        guides(color=FALSE)
      taxon <- colnames(test_table2)[p]
      site <- names(Bodysites_B)[b]
      fp <- paste(taxa_dir, site, taxon, ".pdf", sep="")
      stats_file <- paste(taxa_dir, site, taxon, ".txt", sep="")
      save_plot(fp, plot)
      sink(stats_file)
      print("fdr adjusted pval:\n")
      cat(adjusted[p])
      sink()
    }
  }
}
