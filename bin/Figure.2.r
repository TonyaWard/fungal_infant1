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

taxa_plot <- ggplot(otu, aes_string(x = "SampleID", y = "Relative_Abundance", fill="Taxa")) + 
  geom_bar(stat="identity", position="fill") +
  facet_wrap(~ SuperbodysiteOralSkinNoseVaginaAnalsAureola, scales = "free_x") +
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

pdf(paste(main_fp, "/Figure2.pdf", sep=""), width=6.69, height=5)
plot(fig2)
dev.off()

pdf(paste(main_fp, "/Figure2_Legend.pdf", sep=""), height=6, width=6.69)
print(grid.arrange(legend2))
dev.off()
