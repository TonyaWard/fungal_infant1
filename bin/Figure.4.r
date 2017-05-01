#separation by birth mode
#distance to mom by c-section or vaginal

######################################################################
#plot skin samples according to birth mode
###
beta_table <- wunifrac
body_samples <- Bodysites_B[[1]]
body_samples <- intersect(body_samples, earlies)
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
    geom_point(size = 4, alpha=0.65, aes_string(x = "PC1", y = "PC2", color = "Delivery_Vvaginal_Ccs_IcsInoc")) + 
    scale_color_manual(values=c("#c3c823", "#053d58")) +
    theme_cowplot(font_size = 7) +
    #guides(color=guide_legend(nrow=3)) +
    guides(color=F)  +
    theme(legend.title=element_blank()) +
    labs(x=" ", y=" ") +
    theme(axis.text.x = element_text(color=NA), axis.text.y = element_text(color=NA))

#Make box plots for sides
#Make boxplot of PCs
PC1_boxes <- ggplot(PCOA) +
  geom_boxplot(aes_string(x = "Delivery_Vvaginal_Ccs_IcsInoc", y = "PC2", fill = "Delivery_Vvaginal_Ccs_IcsInoc")) + 
  scale_fill_manual(values=c("#c3c823", "#053d58")) +
  theme_cowplot(font_size = 7) +
  guides(fill=F)+
  coord_flip() +
  labs(x=" ")

PC2_boxes <- ggplot(PCOA) +
  geom_boxplot(aes_string(x ="Delivery_Vvaginal_Ccs_IcsInoc", y = "PC1", fill = "Delivery_Vvaginal_Ccs_IcsInoc")) + 
  scale_fill_manual(values=c("#c3c823", "#053d58")) +
  theme_cowplot(font_size = 7) +
  guides(fill=F) +
  labs(x =" ") +
  theme(axis.text.x = element_text(color=NA))

#########WORKING HERE

top2 <- plot_grid(PC2_boxes, plot1, ncol=2, rel_widths=c(0.3, 1))
bottom2 <- plot_grid(NULL, PC1_boxes, ncol=2, rel_widths=c(0.3, 1))
together2 <- plot_grid(top2, bottom2, nrow=2, rel_heights=c(1, 0.3))



for(i in 1:nrow(taxa_table)){
  header <- rownames(taxa_table)[[i]]
  color_by <- gsub(" ", "_", header)
  print(header)
  if(mean(as.numeric(PCOA[,color_by])) > 0.001){
    test_cor1 <- cor.test(PCOA[,"PC3"], as.numeric(PCOA[,color_by]), method="spearman")
    if(! is.na(test_cor1$p.value) & test_cor1$p.value < 0.05){
      plot1 <- ggplot(PCOA) +
        geom_point(size = 4, alpha=0.65, aes_string(x = "PC2", y = "PC3", color = color_by)) +
        theme_bw() +
        theme(panel.background = element_blank(),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              plot.title = element_text(size = 10),
              legend.title = element_blank(),
              legend.key.size = unit(0.2, "in"),
              legend.text = element_text(size=5),
              legend.position = 'bottom',
              axis.text = element_text(size=5),
              axis.title = element_text(size=8)) +
        scale_color_manual(values=cols_grad2(length(unique(PCOA[,color_by])))) +
        guides(color=guide_legend(nrow=3))
      name1 <- paste(color_by, "1", ".pdf", sep="")
      #fp <- paste(pcoa_dir_baby, name1, sep="")
      pdf(name1, height=4,width=6)
      print(plot1)
      dev.off()
    }
  }
}








######################################################################
#PCoA of mom's Vaginal Samples & Baby Skin, by birth mode

beta_div <- wunifrac
mom <- Vagina
baby_site <- Skin
keeps <- c()
families_plotted <- c()

for(f in 1:length(Families)){
  fam_working <- Families[[f]]
  babies_working <- intersect(fam_working, babies_site)
  momv <- intersect(fam_working, mom)
  if(length(babies_working) > 0 & length(momv) > 0){
    keeps <- c(keeps,babies_working,momv)
    families_plotted <- c(families_plotted,names(Families)[f])
  }
}
beta_table <-  beta_div[keeps, keeps]
PCOA <- pcoa(beta_table)$vectors
for(c in 1:ncol(PCOA)){
  colnames(PCOA)[c] <- paste("PC",c, sep="")
}
PCOA <- cbind(PCOA, rownames(PCOA))
colnames(PCOA)[ncol(PCOA)] <- "SampleID"
PCOA <- merge(PCOA, mapping, by="SampleID")
PCOA$PC1 <- as.numeric(levels(PCOA$PC1))[PCOA$PC1]
PCOA$PC2 <- as.numeric(levels(PCOA$PC2))[PCOA$PC2]

connections <- data.frame(x1 = NA, y1 = NA, x2 = NA, y2 = NA, familyN=NA)
for(t in 1:length(families_plotted)){
  fam_working <- Families[[families_plotted[t]]]
  momv <- intersect(fam_working, mom)
  babies_working <- intersect(fam_working, babies_site)
  for(a in 1:length(babies_working)){
    baby_now <- babies_working[a]
    mom_x <- PCOA[PCOA$SampleID==momv,"PC1"]
    mom_y <- PCOA[PCOA$SampleID==momv,"PC2"]
    baby_x <- PCOA[PCOA$SampleID==baby_now,"PC1"]
    baby_y <- PCOA[PCOA$SampleID==baby_now,"PC2"]
    new_row <- c(mom_x, mom_y, baby_x, baby_y, families_plotted[t])
    connections <- rbind(connections, new_row)
  }
}
connections <- connections[2:nrow(connections),]
connections$x1 <- as.numeric(connections$x1)
connections$x2 <- as.numeric(connections$x2)
connections$y1 <- as.numeric(connections$y1)
connections$y2 <- as.numeric(connections$y2)
v_list <- c()
for(v in 1:length(connections$familyN)){
  working_row <- which(mapping$familyN == connections$familyN[v])[[1]]
  v_list <- c(v_list, mapping[working_row,"Delivery_Vvaginal_Ccs_IcsInoc"])
}
connections$DelMode <- v_list

current_plot <- ggplot() + #plot with lines
  geom_point(data=PCOA,size = 2, alpha=0.65, aes_string(x = "PC1", y = "PC2", color = "Delivery_Vvaginal_Ccs_IcsInoc", shape = "motherorbaby_M_B")) + 
  scale_color_manual(values=c("#c3c823", "#053d58")) +
  theme_cowplot(font_size=7) +
  geom_segment(data=connections,aes(x=x1, y=y1, xend=x2, yend=y2, color=DelMode), alpha=0.65)+
  guides(color=FALSE) +
  guides(shape=FALSE) +
  theme(legend.title = element_blank())


######################################################################
#Average distance to mom's Vaginal Samples to Skin, by birth mode
#Use weighted unifrac

beta_div <- wunifrac

firsts <- c(Days[['1']], Days[['3']])
early_babies <- intersect(firsts, babies)

mom <- Vagina
site <- Skin
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

together <- plot_grid(current_plot, dist_mom, ncol=2, labels="auto", label_size = 12)

pdf(paste(main_fp, "/Figure4.pdf", sep=""), width=3.34, height=2)
plot(together)
dev.off()
