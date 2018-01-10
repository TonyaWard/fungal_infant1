##Making Figure 2
#Taxa Summary of infant samples & alpha diversity

######################################################################
#Make Taxa Summaries by bodysite

b_bodysite <- c(Skin, Oral, Anal_B)
otu <- make_taxa_sums(taxa_table, b_bodysite)
otu$SuperbodysiteOralSkinNoseVaginaAnalsAureola <- factor(otu$SuperbodysiteOralSkinNoseVaginaAnalsAureola, levels=c("Skin", "Oral", "Anal"))

otu2 <- otu[order(otu$subject_id),] 
otu$SampleID <- factor(otu$SampleID, levels=otu2$SampleID)
sample_labels <- table(otu$subject_id)
sample_levels <- names(sample_labels)[order(sample_labels)]
otu$subject_id <- factor(otu$subject_id, levels = sample_levels)


new_cols <- c( "#AC4E6A", "#DC7B3F","#00ffff", "#959D9C", "#6d7e56", "#b5cc90", "#d698ab", "#F0F1E9",  "#AAA3A0", "#300d2e", "#95cec8",
               "#1e3596", "#BC9415", "#4e6aac", "#51707F", "#7F4464", "#E88D3F", "#8f8984",
  "#1B9E77", "#4f2bdb", "#771b9e", "#E7298A", "#edc44d", "#800000", "#c978a3", "#A6761D", "#ff0000", "#7570B3",  
               "#fad4e7", "#D95F02", "#8eb2d2","#66722a", 
               "#66A61E","#a5a2ce" , "#a3c978","#e3e2ef",
               "#b37570", "#d3d3d3")

taxa_plot <- ggplot(otu, aes_string(x = "SampleID", y = "Relative_Abundance", fill="Taxa")) + 
  geom_bar(stat="identity", position="fill") +
  facet_wrap(~ SuperbodysiteOralSkinNoseVaginaAnalsAureola, scales = "free_x") +
  guides(fill=FALSE) +
  theme_cowplot(font_size = 7) +
  #scale_fill_manual(name= names(taxa_cols), values= taxa_cols)+
  scale_fill_manual(values= new_cols) +
  labs(y="Relative Abundance", x='')+
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

#Make separate legend
legend <- ggplot(otu, aes_string(x = "SampleID", y = "Relative_Abundance", fill="Taxa")) + 
  geom_bar(stat="identity", position="fill") +
  facet_wrap(facets=~SuperbodysiteOralSkinNoseVaginaAnalsAureola, scales = "free_x") +
  theme_bw() +
  guides(fill=guide_legend(ncol=3))  +
  theme_cowplot(font_size = 7) +
  #scale_fill_manual(name= names(taxa_cols), values= taxa_cols)+
  scale_fill_manual(values= new_cols) +
  theme(legend.title = element_blank())
g <- ggplotGrob(legend)$grobs
legend2 <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]

######################################################################
#Alpha Div by body site 

working_table <- rbind(alpha_skin, alpha_anal)
working_table <- rbind(alpha_oral, working_table)
working_table$SuperbodysiteOralSkinNoseVaginaAnalsAureola <- factor(working_table$SuperbodysiteOralSkinNoseVaginaAnalsAureola, levels=c("Skin", "Oral", "Anal"))

body_plot <- ggplot(working_table, aes_string(x="SuperbodysiteOralSkinNoseVaginaAnalsAureola", y="observed_species", fill="SuperbodysiteOralSkinNoseVaginaAnalsAureola")) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(position=position_jitter(0.1), shape=1, size=1) +
  theme(legend.position = 'bottom') + 
  theme_cowplot(font_size = 7) +
  labs(x="") +
  guides(fill=F) +
  scale_fill_manual(values=body_cols)


######################################################################
#Compile Figure 2
#Legend must be added after because it behaves weird with plot_grid

fig2 <- ggdraw() +
  draw_plot(body_plot, 0, 0.75, 0.25, 0.25) +
  draw_plot(taxa_plot, 0, 0, 1, 0.75) +
  draw_plot_label(c("a", "b"), x=c(0, .3), y=c(1, 1), size = 12)

pdf(paste(main_fp, "/Figure2.pdf", sep=""), width=6.69, height=5, useDingbats = F)
plot(fig2)
dev.off()

pdf(paste(main_fp, "/Figure2_Legend.pdf", sep=""), height=6, width=6.69, useDingbats=F)
print(grid.arrange(legend2))
dev.off()
