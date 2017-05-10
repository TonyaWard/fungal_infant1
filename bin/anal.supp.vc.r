##Birth mode using anal
#separation by birth mode
#distance to mom by c-section or vaginal

######################################################################
#plot body site samples according to birth mode
#Use only early samples
###
beta_table <- wunifrac
orals <- Bodysites_B[[2]]
anals <- Bodysites_B[[3]]
body_samples <- anals
beta_subset <- beta_table[body_samples,body_samples]
PCOA <- pcoa(beta_subset)$vectors

beta_dist = as.dist(beta_subset)
map2 <- mapping[body_samples,]
ad = adonis(beta_dist ~ map2[,"Delivery_Vvaginal_Ccs_IcsInoc"], data=map2, permutations=999)
p_val <- ad$aov.tab[1,6]
r_sq <- ad$aov.tab[1,5]
#Run Stats for diff. dispersion


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
plot(PCOA$PC2, PCOA$PC3)
plot1 <- ggplot(PCOA) +
  geom_point(size = 2, alpha=0.65, aes_string(x = "PC1", y = "PC2", color = "Delivery_Vvaginal_Ccs_IcsInoc")) + 
  scale_color_manual(values=c("#c3c823", "#053d58")) +
  theme_cowplot(font_size = 7) +
  #guides(color=guide_legend(nrow=3)) +
  guides(color=F)  +
  annotate("text", x=-0.5, y=-0.2, label= paste("P=", p_val), size=2) +
  annotate("text", x=-0.5, y=-0.25, label= paste("R2=", round(r_sq, digits=3)), size=2) +
  theme(legend.title=element_blank()) +
  labs(x=" ", y=" ") +
  theme(axis.text.x = element_text(color=NA), axis.text.y = element_text(color=NA))

#Make box plots for sides
#Make boxplot of PCs
PC1_boxes <- ggplot(PCOA) +
  geom_boxplot(aes_string(x = "Delivery_Vvaginal_Ccs_IcsInoc", y = "PC1", fill = "Delivery_Vvaginal_Ccs_IcsInoc")) + 
  scale_fill_manual(values=c("#c3c823", "#053d58")) +
  theme_cowplot(font_size = 7) +
  guides(fill=F)+
  coord_flip() +
  labs(x=" ")

PC2_boxes <- ggplot(PCOA) +
  geom_boxplot(aes_string(x ="Delivery_Vvaginal_Ccs_IcsInoc", y = "PC2", fill = "Delivery_Vvaginal_Ccs_IcsInoc")) + 
  scale_fill_manual(values=c("#c3c823", "#053d58")) +
  theme_cowplot(font_size = 7) +
  guides(fill=F) +
  labs(x =" ") +
  theme(axis.text.x = element_text(color=NA))

top2 <- plot_grid(PC2_boxes, plot1, ncol=2, rel_widths=c(0.3, 1))
bottom2 <- plot_grid(NULL, PC1_boxes, ncol=2, rel_widths=c(0.3, 1))
together2 <- plot_grid(top2, bottom2, nrow=2, rel_heights=c(1, 0.3))


######################################################################
#Average distance to mom's Vaginal Samples to body site, by birth mode
#Use weighted unifrac

beta_div <- wunifrac

firsts <- c(Days[['1']], Days[['3']])
early_babies <- intersect(firsts, babies)

firsts <- babies #not early only - use all

mom <- Vagina
site <- Anal_B
dist_plots <- c()
v_distances <- t(as.data.frame(c("SampleID", "Distance_MomV")))
v_distances <- ddf <- data.frame(matrix(ncol = 2, nrow = 1))
for(i in 1:length(Families)){
  mom_v <- intersect(Families[[i]], mom)
  if(length(mom_v) > 0){
    babies_working <- intersect(Families[[i]], site)
    babies_working <- intersect(babies_working, firsts)
    if(length(babies_working) > 0){
      for(a in 1:length(babies_working)){
        working_sample <- babies_working[[a]]
        new_row <- c(working_sample, beta_div[working_sample,mom_v])
        v_distances <- rbind(v_distances, new_row)
      }
    }
  }
}


colnames(v_distances) <- c("SampleID", "Distance_MomV")
v_distances <- merge(v_distances, mapping, by="SampleID")
v_distances$Distance_MomV <- as.numeric(as.character(v_distances$Distance_MomV))
wilcox.pval <- wilcox.test(v_distances$Distance_MomV ~ v_distances$Delivery_Vvaginal_Ccs_IcsInoc)$p.value

dist_mom <- ggplot(v_distances, aes(x=Delivery_Vvaginal_Ccs_IcsInoc, y=Distance_MomV, fill=Delivery_Vvaginal_Ccs_IcsInoc)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(position=position_jitter(0.1), shape=1, size=1) +
  labs(y="Distance to Mom", x="Birth Mode") +
  theme_cowplot(font_size=7) +
  annotate("text", x="V", y=1.2, label= paste("(wilcox)P=", round(wilcox.pval, digits=3)), size=2) +
  scale_fill_manual(values=c("#c3c823", "#053d58")) +
  guides(fill=FALSE)

######################################################################

###Plot final
side <- plot_grid(dist_mom, NULL, nrow=2)
final_plot <- plot_grid(together2, side, ncol=2, rel_widths = c(0.55, 0.45), labels="auto", label_size = 12)
pdf(paste(main_fp, "/Supplement_Anal_VC.pdf", sep=""), width=6.64, height=3)
plot(final_plot)
dev.off()


