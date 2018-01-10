##Making Figure 1
#PCoA of Bodysites, and colored by Candida spp. and boxplots of RA

######################################################################
#Make PCOA of all infants colored by bodysite (from the PCOA.r script)
#we will use weighted Unifrac

#Set beta table, subset to just babies
beta_table <- wunifrac
beta_subset <- beta_table[babies,babies]
PCOA <- pcoa(beta_subset)$vectors
var_exp <- pcoa(beta_subset)$values
map2 <- mapping[babies,]
#Run stats for diff. centroids
beta_dist = as.dist(beta_subset)
ad = adonis(beta_dist ~ map2[,"SuperbodysiteOralSkinNoseVaginaAnalsAureola"], data=map2, permutations=999)
p_val <- ad$aov.tab[1,6]
r_sq <- ad$aov.tab[1,5]
#Run Stats for diff. dispersion
beta_out <- betadisper(beta_dist, map2$SuperbodysiteOralSkinNoseVaginaAnalsAureola)
p_val_disp <- permutest(beta_out)$tab[1, 6]


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
  geom_point(size = 2, alpha=0.65, aes_string(x = "PC1", y = "PC2", color = "SuperbodysiteOralSkinNoseVaginaAnalsAureola")) + 
  scale_color_manual(values=body_cols2) +
  theme_cowplot(font_size = 7) +
  guides(color=F) +
  annotate("text", x=-0.6, y=-0.2, label= paste("P=", p_val), size=2) +
  annotate("text", x=-0.6, y=-0.25, label= paste("R2=", round(r_sq, digits=3)), size=2) +
  #guides(color=guide_legend(nrow=3)) +
  labs(x="", y="") +
  theme(axis.text.x = element_text(color=NA), axis.text.y = element_text(color=NA))

#Make boxplot of PCs
PC1_boxes <- ggplot(PCOA) +
  geom_boxplot(aes_string(x = factor(PCOA$SuperbodysiteOralSkinNoseVaginaAnalsAureola, levels=c("Anal", "Oral", "Skin")), y = "PC1", fill = "SuperbodysiteOralSkinNoseVaginaAnalsAureola")) + 
  scale_fill_manual(values=body_cols2) +
  theme_cowplot(font_size = 7) +
  guides(fill=F)+
  coord_flip() +
  labs(x="", y= paste("PC1 (", round(var_exp$Rel_corr_eig[1], digits=3)*100, "%)", sep=""))


PC2_boxes <- ggplot(PCOA) +
  geom_boxplot(aes_string(x =factor(PCOA$SuperbodysiteOralSkinNoseVaginaAnalsAureola, levels=c("Anal", "Oral", "Skin")), y = "PC2", fill = "SuperbodysiteOralSkinNoseVaginaAnalsAureola")) + 
  scale_fill_manual(values=body_cols2) +
  theme_cowplot(font_size = 7) +
  guides(fill=F) +
  labs(x ="",  y= paste("PC2 (", round(var_exp$Rel_corr_eig[2], digits=3)*100, "%)", sep="")) +
  theme(axis.text.x = element_text(color=NA))

#Compile the PCoA and boxes

top2 <- plot_grid(PC2_boxes, body_PCOA, ncol=2, rel_widths=c(0.3, 1))
bottom2 <- plot_grid(NULL, PC1_boxes, ncol=2, rel_widths=c(0.3, 1))
together2 <- plot_grid(top2, bottom2, nrow=2, rel_heights=c(1, 0.3))

######################################################################
#For the sample PCoA, color by relative abundance of Candida spp.
#These are significant gradients pulled out of the earlier 
#exhaustive analysis

PCOA$Candida_albicans <- as.numeric(PCOA$Candida_albicans)
albicans_PCOA <- ggplot(PCOA) +
  geom_point(size = 1.5, alpha=0.85, aes(x = PC1, y = PC2, color = Candida_albicans)) + 
  scale_color_gradient(low="#f49f3f", high= "#3f1f6b", guide="colorbar") +
  geom_point(data= PCOA[PCOA$SuperbodysiteOralSkinNoseVaginaAnalsAureola == "Anal",], alpha=0.55, size = 1.5, pch=21, aes(x = PC1, y = PC2)) + 
  theme_cowplot(font_size = 7) +
  guides(color=F)

PCOA$Candida_parapsilosis <- as.numeric(PCOA$Candida_parapsilosis)
parap_PCOA <- ggplot(PCOA) +
  geom_point(size = 1.5, alpha=0.85, aes(x = PC1, y = PC2, color = Candida_parapsilosis)) +
  scale_color_gradient(low="#f49f3f", high= "#3f1f6b", guide="colorbar") +
  geom_point(data= PCOA[PCOA$SuperbodysiteOralSkinNoseVaginaAnalsAureola == "Anal",], alpha=0.55, size = 1.5, pch=21, aes(x = PC1, y = PC2)) + 
  theme_cowplot(font_size = 7) +
  guides(color=F)

#Print with legend and add this later
forlegend <- ggplot(PCOA) +
  geom_point(size = 1.5, alpha=0.85, aes_string(x = "PC1", y = "PC2", color = "Candida_parapsilosis")) + 
  scale_color_gradient(low="#f49f3f", high= "#3f1f6b", guide="colorbar") +
  theme_cowplot(font_size = 7)
pdf(paste(main_fp, "/Fig1_grad_legend.pdf", sep=""), height=6, width=6.69)
print(forlegend)
dev.off()

taxa_pcoas <- plot_grid(albicans_PCOA, parap_PCOA, ncol=2)

######################################################################
#Add body site diff taxa
kids <- c(Bodysites_B[[1]], Bodysites_B[[2]], Bodysites_B[[3]])
test_table <- taxa_table[, kids]
test_table <- test_table[rowSums(test_table) > 0,]
test_table <- test_table[rowSums(test_table > 0 )/ncol(test_table) > 0.1, ]
map_test <- mapping[colnames(test_table),]
map_test$SuperbodysiteOralSkinNoseVaginaAnalsAureola <- factor(map_test$SuperbodysiteOralSkinNoseVaginaAnalsAureola, levels=c("Skin", "Oral", "Anal"))
working_table <- data.frame(t(test_table))

taxon1 <- ggplot(working_table, aes(y=working_table[,"Candida.albicans"], x=map_test$SuperbodysiteOralSkinNoseVaginaAnalsAureola, fill=map_test$SuperbodysiteOralSkinNoseVaginaAnalsAureola)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(position=position_jitter(0.1), shape=1, size=1) +
  theme(legend.position = 'bottom') + 
  theme_cowplot(font_size = 7) +
  labs(x="", y = "C. albicans Relative Abundance") +
  guides(fill=F) +
  scale_fill_manual(values=body_cols)

taxon2 <- ggplot(working_table, aes(y=working_table[,"Candida.parapsilosis"], x=map_test$SuperbodysiteOralSkinNoseVaginaAnalsAureola, fill=map_test$SuperbodysiteOralSkinNoseVaginaAnalsAureola)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(position=position_jitter(0.1), shape=1, size=1) +
  theme(legend.position = 'bottom') + 
  theme_cowplot(font_size = 7) +
  guides(fill=F) +
  labs(x="", y = "C. parapsilosis Relative Abundance") +
  scale_fill_manual(values=body_cols)

boxes <- plot_grid(taxon1,taxon2, ncol=2)

taxa_together <- plot_grid(taxa_pcoas, boxes, nrow=2)
######################################################################
#Compile Figure

right_side <- plot_grid(albicans_PCOA, taxon1, taxon2, body_plot, ncol=3)
right_side2 <- plot_grid(right_side, NULL, nrow=2, rel_heights = c(0.85, 0.15))

fig1 <- ggdraw() +
  draw_plot(PC2_boxes, 0, 0.2, 0.12, 0.75) +
  draw_plot(PC1_boxes, 0.08, 0.0, 0.4, 0.2)+
  draw_plot(body_PCOA, 0.08, 0.2, 0.4, 0.75) +
  draw_plot(taxa_together, 0.5, 0, 0.5, 0.95) +
  draw_plot_label(c("a", "b", "c", "d", "e"), c(0,0.5,0.75,0.5,0.75), c(1,1,1,0.6,0.6), size = 12)

pdf(paste(main_fp, "/Figure1.pdf", sep=""), height=3.3, width=6.6)
print(fig1)
dev.off()
