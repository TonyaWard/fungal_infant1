#### Test differences for within and between group distances for beta
######################################################################
#For PCoA by body sites sep, color by infant
pcoa_dir_baby <- paste(main_fp, "beta_div/PCOA_baby/", sep='/')
names(baby_colors) <- unique(mapping$subject_id)

for(b in 1:length(beta_tables)){
  beta_div <- beta_tables[[b]]
  beta_name <- beta_metrics[b]
  plot_list <- c()
  
  #make stats file
  file_name <- paste(pcoa_dir_baby, "Beta_Stats_", beta_name, ".txt", sep='')
  sink(file_name)
  sink()
  
  for(a in 1:length(Bodysites_B)){
    site <- Bodysites_B[[a]]
    beta_subset <- beta_div[site,site]
    
    #Set 1 will be one infant
    #Set 2 will be all other infants
    infant_ids <- unique(mapping[rownames(beta_subset),"subject_id"])
    within <- c()
    between <- c()
    for(n in 1:length(infant_ids)){
      set1 <- intersect(rownames(beta_subset), rownames(mapping[mapping$subject_id == infant_ids[n],]))
      set2 <- intersect(rownames(beta_subset), rownames(mapping[!mapping$subject_id == infant_ids[n],]))
      full_set <- c(set1, set2)
      set1_within <- c()
      between_sets <- c()
          
      if(length(set1) > 2 && length(set2) > 2){
        for(c in 1:length(colnames(beta_div))){
          for(r in (c+1):length(rownames(beta_div))){
            if(rownames(beta_div)[r] %in% set1 && colnames(beta_div)[c] %in% set1){
              set1_within <- c(set1_within, beta_div[r,c])
            } else {
              if(rownames(beta_div)[r] %in% set1 && colnames(beta_div)[c] %in% set2){
                between_sets <- c(between_sets, beta_div[r,c])
              } 
            }
          }
        }
      }
      within <- c(within, set1_within)
      between <- c(between_sets, between)
    }
     
    sets_test <- list(within, between)
    names(sets_test) <- c("Within", "Between")
    wtest <- wilcox.test(within, between)
    pval <- wtest$p.value
    distances <- melt(sets_test)
          
    #print stats to screen
    cat(sprintf('\n%s,%s:\n',names(Bodysites_B)[a], beta_name))
    print(wtest)
          
    #write stats to file
    sink(file_name, append =TRUE)
    cat(sprintf('\n%s,%s:\n',names(Bodysites_B)[a], beta_name))
    print(wtest)
    sink()
          
    #assign pdf name for plot
    name1 <- paste(pcoa_dir_baby, "Beta_dists_", beta_name, "_", names(Bodysites_B)[a], ".pdf", sep='')
    pdf(name1, height=4,width=6)
          
    #make beta div box plots
    plot1 <- ggplot(distances, aes(x=L1, y=value, fill=L1)) +
      geom_boxplot(outlier.shape = NA) +
      #geom_jitter(position=position_jitter(0.1), shape=1, size=1) +
      theme_cowplot(font_size=7) +
      annotate("text", x="Between", y=1.0, label= paste("(wilcox)P=", round(pval, digits=3)), size=2) +
      guides(fill=FALSE) +
      labs(x="", y = paste(beta_name, "distance")) +
      scale_fill_manual(values= c("#E88D3F", "#5F8CA2"))
    plot(plot1)
    dev.off()
  }
}
    


