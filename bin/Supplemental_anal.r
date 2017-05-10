#Distance to A (not paired) over time
# for each infant sample, plot the average distance of that sample to all adult A samples

beta_div <- wunifrac
plot_list_r_p <- c()
distance_table <- t(data.frame(c("bodysite", "subject", "day", "distance", "SD", "r_or_p")))

site <- Anal_B
for(i in 1:length(Families)){
  babies_working <- intersect(Families[[i]], babies)
  babies_working <- intersect(babies_working, site)
  for(d in 1:(length(Days))){
    day_now <- Days[[d]]
    working <- intersect(day_now, babies_working)
    if(length(working) == 0){
      distances <- c("anal", names(Families)[i], names(Days)[d], NA, NA, "paired")
      distance_table <- rbind(distance_table, distances)
    } else {
      av_dist <- mean(beta_div[Anal_M,working])
      dist_dev <- sd(beta_div[Anal_M,working])
      distances <- c("anal", names(Families)[i], names(Days)[d], av_dist, dist_dev, "paired")
      distance_table <- rbind(distance_table, distances)
    }
  }
}
colnames(distance_table) <- c("Bodysite", "Family", "Day", "Distance", "StDev", "R_P")
rownames(distance_table) <- c(1:nrow(distance_table))
distance_table <- as.data.frame(distance_table[2:nrow(distance_table),])
distance_table$Day <- as.numeric(as.character(distance_table$Day))
distance_table$Distance <- as.numeric(as.character(distance_table$Distance))
distance_table$R_P <- as.character(distance_table$R_P)
distance_table <- distance_table[complete.cases(distance_table),]
alpha <- distance_table$Distance
subject <- distance_table$Family
time <- distance_table$Day
obs <- -mean(sapply(split(1:nrow(distance_table), subject), 
                    function(ixx) if(length(ixx) < 3) 0 else cor.test(alpha[ixx], time[ixx], method='spear')$statistic))
mc.stats <- -replicate(999,mean(sapply(split(1:nrow(distance_table), subject), 
                                       function(ixx) if(length(ixx) < 3) 0 else cor.test(alpha[ixx], sample(time[ixx]), method='spear')$statistic)))
pval1 <- mean(c(obs,mc.stats) >= obs)

this_plot <- ggplot(distance_table, aes(x=Day, y=Distance)) +
  geom_jitter(alpha=0.65, size=2, width = 0.25, color=body_cols2[1]) +
  geom_smooth(method=lm, se=FALSE, linetype = "dashed", size=0.5, color=body_cols2[1]) +
  theme_cowplot(font_size = 7) +
  annotate("text", x=10, y=0.92, label= paste("P=", round(pval1, digits=3)), size=2, color= body_cols2[1]) +
  guides(color=FALSE) +
  labs(y="Distance to Adult") 

######################################################################
#PCoA of just A samples, mom and baby
#Set beta table, subset to just babies
beta_table <- wunifrac
keeps <- union(Anal_B, Anal_M)
beta_subset <- beta_table[keeps,keeps]
PCOA <- pcoa(beta_subset)$vectors
map2 <- mapping[keeps,]
#Run stats for diff. centroids
beta_dist = as.dist(beta_subset)
ad = adonis(beta_dist ~ map2[,"motherorbaby_M_B"], data=map2, permutations=999)
p_val <- ad$aov.tab[1,6]
r_sq <- ad$aov.tab[1,5]
#Run Stats for diff. dispersion
beta_out <- betadisper(beta_dist, map2$motherorbaby_M_B)
p_val_disp <- permutest(beta_out)$tab[1, 6]

mom_baby_cols <- c("#CB1B4A", "#99897e")
#Make PCOA
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

#Make final plot
body_PCOA <- ggplot(PCOA) +
  geom_point(size = 2, alpha=0.65, aes_string(x = "PC1", y = "PC2", color = "motherorbaby_M_B")) + 
  scale_color_manual(values=mom_baby_cols) +
  theme_cowplot(font_size = 7) +
  guides(color=F) +
  annotate("text", x=-0.3, y=-0.2, label= paste("P=", p_val), size=2) +
  annotate("text", x=-0.3, y=-0.25, label= paste("R2=", round(r_sq, digits=3)), size=2) +
  #guides(color=guide_legend(nrow=3)) +
  labs(x="", y="") +
  theme(axis.text.x = element_text(color=NA))

#Make boxplot of PCs
PC1_boxes <- ggplot(PCOA) +
  geom_boxplot(aes_string(x = "motherorbaby_M_B", y = "PC1", fill = "motherorbaby_M_B")) + 
  scale_fill_manual(values=mom_baby_cols) +
  theme_cowplot(font_size = 7) +
  guides(fill=F)+
  coord_flip() +
  labs(x="")

PC2_boxes <- ggplot(PCOA) +
  geom_boxplot(aes_string(x ="motherorbaby_M_B", y = "PC2", fill = "motherorbaby_M_B")) + 
  scale_fill_manual(values=mom_baby_cols) +
  theme_cowplot(font_size = 7) +
  guides(fill=F) +
  labs(x ="") +
  theme(axis.text.x = element_text(color=NA))

top2 <- plot_grid(PC2_boxes, body_PCOA, ncol=2, rel_widths=c(0.3, 1))
bottom2 <- plot_grid(NULL, PC1_boxes, ncol=2, rel_widths=c(0.31, 1))
together2 <- plot_grid(top2, bottom2, nrow=2, rel_heights=c(1, 0.3))

total_plot <- plot_grid(together2, this_plot, nrow=2, labels="auto", label_size=12, rel_heights=c(0.65,0.35))
pdf(paste(main_fp, "/supplemental_anal.pdf", sep=""), width=3.34, height =4)
plot(total_plot)
dev.off()
