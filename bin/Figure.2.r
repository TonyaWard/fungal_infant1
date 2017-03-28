##Making Figure 2
#Alpha diversity over time
#PCoAs of each bodysite colored by infant
#Distance to self over time

######################################################################
#Alpha diversity over time
#Use shannon

alpha_metric <- "shannon"
alpha_plots <- c()
min_a <- min(alpha_table_baby[,alpha_metric])
max_a <- max(alpha_table_baby[,alpha_metric])
for(t in 1:length(alpha_tissues)){
  working_alpha <- melt(alpha_tissues[t], id.vars = c('SampleID', 'Delivery_Vvaginal_Ccs_IcsInoc', 'planned_sampling_day_0_1_3_7_14_30_60_90___365', 'subject_id'), measure.vars = c(alpha_metric))
  working_alpha$Collection_Day <- as.numeric(working_alpha$planned_sampling_day_0_1_3_7_14_30_60_90___365)
  working_alpha$subject_id <- factor(working_alpha$subject_id)
  colnames(working_alpha)[6] <- alpha_metric
  fm2 <- lme(shannon ~ Collection_Day, data=working_alpha, random= ~ 1 | subject_id)
  alpha <- working_alpha$shannon
  subject <- working_alpha$subject_id
  time <- working_alpha$Collection_Day
  obs <- mean(sapply(split(1:nrow(working_alpha), subject), 
         function(ixx) if(length(ixx) < 3) 0 else cor.test(alpha[ixx], time[ixx], method='spear')$statistic))
  
  mc.stats <- replicate(999,mean(sapply(split(1:nrow(working_alpha), subject), 
                     function(ixx) if(length(ixx) < 3) 0 else cor.test(alpha[ixx], sample(time[ixx]), method='spear')$statistic)))

  print(mean(c(obs,mc.stats) >= obs))
  # Using a permutation test and a groupwise average of the test statistic of the spearman correlation, we find clinical relevance but not significance.
  #fm1 <- lme(shannon ~ Collection_Day, data=working_alpha, random= ~ Collection_Day | subject_id)
  #anova(fm2, fm1) #use fm2 based on anova
  r_sq <- fm2$coefficients["fixed"][[1]][[2]]
  pval <- data.frame(coef(summary(fm2)))[2,"p.value"]
  
  figure <- ggplot(working_alpha, aes_string(x="Collection_Day", y=alpha_metric)) +
    geom_jitter(alpha=0.65, size=2, color=body_cols[t], width = 0.25) +
    theme_cowplot(font_size = 7) +
    geom_smooth(method=lm, se=FALSE, color="#99897E", linetype = "dashed", size=0.5)+
    guides(fill=FALSE)+
    annotate("text", x=25, y=3.25, label= paste("P=", round(pval, digits=3)), size=2) +
    annotate("text", x=25, y=3, label= paste("R2=", round(r_sq, digits=3)), size=2) +
    labs(x="Day", y=alpha_metric )+
    expand_limits(y=c(min_a,max_a))
  name <- names(alpha_tissues[t])
  alpha_plots[[name]] <- figure
}
alpha_time <- plot_grid(alpha_plots[[1]], alpha_plots[[2]], alpha_plots[[3]],ncol=3)





######################################################################
#Dan's way? 
#mean  slope
#permute 999 times within subject
#each point keeps it's x-value, but scramble the y (diversity)
#what fraction of the time is the m.s. greater than permuted one


######################################################################
#PCoA for each bodysite, colored by infant
baby_colors <- sample(cols2(length(unique(mapping$subject_id))))
names(baby_colors) <- unique(mapping$subject_id)

#Use weighted unifrac
beta_div <- wunifrac
plot_list <- c()
for(a in 1:length(Bodysites_B)){
  site <- Bodysites_B[[a]]
  beta_subset <- beta_div[site,site]
  PCOA <- pcoa(beta_subset)$vectors
  for(i in 1:ncol(PCOA)){
    colnames(PCOA)[i] <- paste("PC",i, sep="")
  }
    
  PCOA <- cbind(PCOA, rownames(PCOA))
  colnames(PCOA)[ncol(PCOA)] <- "SampleID"
  mapping2 <- mapping
  mapping2 <- data.frame(lapply(mapping2, as.character), stringsAsFactors=FALSE)
  PCOA <- merge(PCOA, mapping2, by="SampleID")
  PCOA$PC1 <- as.numeric(levels(PCOA$PC1))[PCOA$PC1]
  PCOA$PC2 <- as.numeric(levels(PCOA$PC2))[PCOA$PC2]
  PCOA$PC3 <- as.numeric(levels(PCOA$PC3))[PCOA$PC3]
  PCOA$PC4 <- as.numeric(levels(PCOA$PC4))[PCOA$PC4]
  
  #Run stats for diff. centroids
  beta_dist = as.dist(beta_subset)
  map2 <- mapping[site,]
  ad = adonis(beta_dist ~ map2[,"subject_id"], data=map2, permutations=999)
  p_val <- ad$aov.tab[1,6]
  r_sq <- ad$aov.tab[1,5]
  #Run Stats for diff. dispersion
  beta_out <- betadisper(beta_dist, map2$subject_id)
  p_val_disp <- permutest(beta_out)$tab[1, 6]
    
  plot1 <- ggplot(PCOA) +
    geom_point(size = 2, alpha=0.65, aes_string(x = "PC1", y = "PC2", color = "subject_id")) + 
    scale_color_manual(values=baby_colors) +
    guides(color=F) +
    #guides(color=guide_legend(nrow = 2)) +
    annotate("text", x=-0.6, y=-0.2, label= paste("P=", p_val), size=2) +
    annotate("text", x=-0.6, y=-0.25, label= paste("R2=", round(r_sq, digits=3)), size=2) +
    theme_cowplot(font_size = 7) +
    theme(legend.position = "bottom")
  name <- names(Bodysites_B)[a]
  plot_list[[name]] <- plot1
}
plots <- plot_grid(plot_list[[1]], plot_list[[2]], plot_list[[3]],ncol=3)

#make family legend separate, add to final figure after
legend <- ggplot(PCOA) + 
  geom_point(size = 2, alpha=0.65, aes_string(x = "PC1", y = "PC2", color = "subject_id")) + 
  scale_color_manual(values=baby_colors) +
  theme_cowplot(font_size = 7) +
  guides(color=guide_legend(nrow=2))  +
  theme(legend.title = element_blank())

g <- ggplotGrob(legend)$grobs
legend2 <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]

pdf(paste(main_fp, "/Figure2_legend.pdf", sep=""), height=2, width=5)
print(grid.arrange(legend2))
dev.off()

######################################################################
#Distance to self by body site
#Use weighted unifrac

timeline_subjects <- as.character(unique(mapping[mapping$planned_sampling_day_0_1_3_7_14_30_60_90___365 == '30',"familyN"]))
mapping$familyN <- as.character(mapping$familyN)


beta_div <- wunifrac
distance_table <- t(data.frame(c("bodysite", "subject", "day", "distance", "sample")))
body_names <- c("Forehead","Oralmucosa","Anal")
for(b in 1:length(Bodysites_B)){
  body_now <- Bodysites_B[[b]]
  for(s in 1:length(timeline_subjects)){
    baby_now <- rownames(mapping[mapping$familyN == timeline_subjects[s],])
    thirty <- rownames(mapping[mapping$planned_sampling_day_0_1_3_7_14_30_60_90___365 == 30 & mapping$familyN == timeline_subjects[s] & mapping$bodysite_OralMucosa_Forehead_VolarRight_PalmRight_FootRight_Vaginal_Anal_Feces == body_names[b],])
    for(d in 1:(length(Days)-1)){
      day_now <- Days[[d]]
      working <- intersect(day_now, baby_now)
      working <- intersect(working, body_now)
      if(length(working) == 0){
        distances <- c(names(Bodysites_B)[b], timeline_subjects[s], names(Days)[d], NA, NA)
        distance_table <- rbind(distance_table, distances)
      } else {
        distances <- c(names(Bodysites_B)[b], timeline_subjects[s], names(Days)[d], beta_div[working[[1]], thirty], working[[1]])
        distance_table <- rbind(distance_table, distances)
      }
    }
  }
}
colnames(distance_table) <- c("Bodysite", "Subject", "Day", "Distance", "SampleID")
rownames(distance_table) <- c(1:nrow(distance_table))
distance_table <- as.data.frame(distance_table[2:nrow(distance_table),])
distance_table <- data.frame(matrix(unlist(distance_table), nrow=nrow(distance_table)),stringsAsFactors=FALSE)
colnames(distance_table) <- c("Bodysite", "Subject", "Day", "Distance", "SampleID")
distance_table$Day <- as.numeric(as.character(distance_table$Day))
distance_table$Distance <- as.numeric(as.character(distance_table$Distance))
distance_table$Subject <- as.character(distance_table$Subject)
distance_table$Bodysite <- as.character(distance_table$Bodysite)
max_distance <- max(distance_table$Distance, na.rm=T)
min_distance <- max(distance_table$Distance, na.rm=T)

fig_list <- c()
for(x in 1:length(unique(distance_table$Bodysite))){
  site <- unique(distance_table$Bodysite)[x]
  working_dists <- distance_table[distance_table$Bodysite == site,]
  working_dists <- working_dists[complete.cases(working_dists),]
  fm2 <- lme(Distance ~ Day, data=working_dists, random= ~ 1 | Subject)
  r_sq <- fm2$coefficients["fixed"][[1]][[2]]
  pval <- data.frame(coef(summary(fm2)))[2,"p.value"]
  
  figure <- ggplot(distance_table[distance_table$Bodysite == site,], aes(x=Day, y=Distance)) +
    geom_jitter(alpha=0.65, size=2, color=body_cols[x], width = 0.25) +
    geom_smooth(method=lm, se=FALSE, color="#99897E", linetype = "dashed", size=0.5) +
    theme_cowplot(font_size = 7) +
    expand_limits(y=c(min_distance,max_distance)) +
    annotate("text", x=15, y=1.18, label= paste("P=", round(pval, digits=3)), size=2) +
    annotate("text", x=15, y=1.1, label= paste("R2=", round(r_sq, digits=3)), size=2) +
    guides(fill=FALSE) +
    labs(y="Distance to Self")
  fig_list[[site]] <- figure
}
distance_self <- plot_grid(fig_list[["Skin"]], fig_list[["Oral"]], fig_list[["Anal_B"]], ncol=3)

######################################################################
#Compile Figure 2

together <- plot_grid(alpha_time, plots, distance_self,nrow=3, labels="auto", rel_heights = c(1,1,1), label_size = 12)
#save_plot(paste(main_fp, "/Figure2.pdf", sep=""), together)
pdf(paste(main_fp, "/Figure2.pdf", sep=""), width=6.69, height =6)
plot(together)
dev.off()

