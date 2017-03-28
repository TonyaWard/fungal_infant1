##Making Figure 1
#PCoA of Bodysites
#Taxa Summaries and Lengend

######################################################################
#Make PCOA of all infants colored by bodysite (from the PCOA.r script)
#we will use weighted Unifrac

#Set beta table, subset to just babies
beta_table <- wunifrac
beta_subset <- beta_table[babies,babies]
PCOA <- pcoa(beta_subset)$vectors
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
  geom_boxplot(aes_string(x = "SuperbodysiteOralSkinNoseVaginaAnalsAureola", y = "PC1", fill = "SuperbodysiteOralSkinNoseVaginaAnalsAureola")) + 
  scale_fill_manual(values=body_cols2) +
  theme_cowplot(font_size = 7) +
  guides(fill=F)+
  coord_flip() +
  labs(x="")


PC2_boxes <- ggplot(PCOA) +
  geom_boxplot(aes_string(x =factor(PCOA$SuperbodysiteOralSkinNoseVaginaAnalsAureola, levels=c("Skin", "Oral", "Anal")), y = "PC2", fill = "SuperbodysiteOralSkinNoseVaginaAnalsAureola")) + 
  scale_fill_manual(values=body_cols2) +
  theme_cowplot(font_size = 7) +
  guides(fill=F) +
  labs(x ="") +
  theme(axis.text.x = element_text(color=NA))

#########WORKING HERE

top2 <- plot_grid(PC2_boxes, body_PCOA, ncol=2, rel_widths=c(0.3, 1))
bottom2 <- plot_grid(NULL, PC1_boxes, ncol=2, rel_widths=c(0.3, 1))
together2 <- plot_grid(top2, bottom2, nrow=2, rel_heights=c(1, 0.3))

######################################################################
#Make Taxa Summaries by bodysite

b_bodysite <- c(Skin, Oral, Anal_B)
otu <- make_taxa_sums(taxa_table, b_bodysite)
otu$SuperbodysiteOralSkinNoseVaginaAnalsAureola <- factor(otu$SuperbodysiteOralSkinNoseVaginaAnalsAureola, levels=c("Skin", "Oral", "Anal"))

taxa_plot <- ggplot(otu, aes_string(x = "SampleID", y = "Relative_Abundance", fill="Taxa")) + 
  geom_bar(stat="identity", position="fill") +
  facet_wrap(facets=~SuperbodysiteOralSkinNoseVaginaAnalsAureola, scales = "free_x") +
  guides(fill=FALSE) +
  theme_cowplot(font_size = 7) +
  scale_fill_manual(name= names(taxa_cols), values= taxa_cols)+
  labs(y="Relative Abundance", x='')+
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

#Make separate legend
legend <- ggplot(otu, aes_string(x = "SampleID", y = "Relative_Abundance", fill="Taxa")) + 
  geom_bar(stat="identity", position="fill") +
  facet_wrap(facets=~SuperbodysiteOralSkinNoseVaginaAnalsAureola, scales = "free_x") +
  theme_bw() +
  guides(fill=guide_legend(ncol=3))  +
  theme_cowplot(font_size = 7) +
  scale_fill_manual(name= names(taxa_cols), values= taxa_cols)+
  theme(legend.title = element_blank())
g <- ggplotGrob(legend)$grobs
legend2 <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]

######################################################################
#Compile Figure 1
#Legend must be added after becuase it behaves weird with plot_grid

fig1 <- ggdraw() +
  #draw_plot(together2, 0, 0.5, 0.5, 0.5) +
  draw_plot(PC2_boxes, 0, 0.7, 0.12, 0.3) +
  draw_plot(PC1_boxes, 0.08, 0.6, 0.3, 0.12)+
  draw_plot(body_PCOA, 0.08, 0.7, 0.3,0.3) +
  draw_plot(taxa_plot, 0, 0, 1, 0.6) +
  draw_plot_label(c("a", "b"), c(0, 0), c(1, 0.6), size = 15)

top <- plot_grid(together2, NULL, ncol=2, rel_widths = c(1,1))
together <- plot_grid(top, taxa_plot, nrow=2, rel_heights = c(1,1), labels="auto", label_size = 12)

pdf(paste(main_fp, "/Figure1.pdf", sep=""), width=6.69, height=6)
plot(fig1)
dev.off()

pdf(paste(main_fp, "/Figure1_Legend.pdf", sep=""), height=6, width=6.69)
print(grid.arrange(legend2))
dev.off()
