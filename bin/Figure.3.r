##Making Figure 3
#separation by birth mode
#distance to mom by c-section or vaginal

######################################################################
#plot skin samples according to birth mode
#Use only early samples
###
beta_table <- wunifrac
body_samples <- Bodysites_B[[1]]
#body_samples <- intersect(body_samples, earlies)
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


#This is Candida tropicalis on PC2
for(i in 1:nrow(taxa_table)){
  header <- rownames(taxa_table)[[i]]
  color_by <- gsub(" ", "_", header)
  print(header)
  if(mean(as.numeric(PCOA[,color_by])) > 0.001){
    test_cor1 <- cor.test(PCOA[,"PC1"], as.numeric(PCOA[,color_by]), method="spearman")
    if(! is.na(test_cor1$p.value) & test_cor1$p.value < 0.05){
      plot1 <- ggplot(PCOA) +
        geom_point(size = 4, alpha=0.65, aes_string(x = "PC1", y = "PC2", color = color_by)) +
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
      name1 <- paste("birthmode_skin", name1, sep="")
      #fp <- paste(pcoa_dir_baby, name1, sep="")
      pdf(paste(main_fp, "/", name1, sep=""), height=4,width=6)
      print(plot1)
      dev.off()
    }
  }
}


######################################################################
#PCoA of mom's Vaginal Samples & Baby Skin, by birth mode
#NOT USING THIS

beta_div <- wunifrac
mom <- Vagina
baby_site <- Skin
#babies_site <- intersect(baby_site, earlies)
babies_site <- baby_site
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
PCOA$PC3 <- as.numeric(levels(PCOA$PC3))[PCOA$PC3]
PCOA$PC4 <- as.numeric(levels(PCOA$PC4))[PCOA$PC4]

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

PCOA$birth_sampletype <- paste(PCOA$motherorbaby_M_B, PCOA$Delivery_Vvaginal_Ccs_IcsInoc, sep="_")
PCOA[PCOA$birth_sampletype == "M_C", "birth_sampletype"] <- "M_V"

current_plot <- ggplot() + #plot with lines
  geom_point(data=PCOA,size = 2, alpha=0.65, aes_string(x = "PC1", y = "PC2", color = "birth_sampletype", shape = "motherorbaby_M_B")) + 
  scale_color_manual(values=c("#c3c823", "#053d58", "#b75f6d")) +
  theme_cowplot(font_size=7) +
  #geom_segment(data=connections,aes(x=x1, y=y1, xend=x2, yend=y2, color=DelMode), alpha=0.65)+
  guides(color=FALSE) +
  #guides(shape=FALSE) +
  theme(legend.title = element_blank())


######################################################################
#Average distance to mom's Vaginal Samples to Skin, by birth mode
#Use weighted unifrac

beta_div <- wunifrac

firsts <- c(Days[['1']], Days[['3']])
early_babies <- intersect(firsts, babies)

firsts <- babies #not early only - use all

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


######################################################################
#### Distance to mom over time for Vs vs Cs
#NOT USING THIS
beta_div <- wunifrac
plot_list_r_p <- c()
distance_table <- t(data.frame(c("bodysite", "subject", "day", "distance", "sample", "r_or_p")))
baby_sites <- list(Skin)
baby_site_names <- c("skin")
mom_sites <- list(Vagina)

for(s in 1:length(mom_sites)){
  site <- baby_sites[[s]]
  mom <- mom_sites[[s]]
  for(i in 1:length(Families)){
    mom_d <- intersect(Families[[i]], mom)
    babies_working <- intersect(Families[[i]], babies)
    babies_working <- intersect(babies_working, site)
    if(length(mom_d) > 0){
      for(d in 1:length(Days)){
        day_now <- Days[[d]]
        working <- intersect(day_now, babies_working)
        if(length(working) == 0){
          if(babies_working[i] %in% Vaginal){
            distances <- c(baby_site_names[s], names(Families)[i], names(Days)[d], NA, NA, "Vaginal")
            distance_table <- rbind(distance_table, distances)
          } else {
            distances <- c(baby_site_names[s], names(Families)[i], names(Days)[d], NA, NA, "Csection")
            distance_table <- rbind(distance_table, distances)
          } 
        }
        if(length(working) > 0){
          if(babies_working[1] %in% Vaginal){
            distances <- c(baby_site_names[s], names(Families)[i], names(Days)[d], beta_div[working[[1]], mom_d], working[[1]], "Vaginal")
            distance_table <- rbind(distance_table, distances)
          } else {
            distances <- c(baby_site_names[s], names(Families)[i], names(Days)[d], beta_div[working[[1]], mom_d], working[[1]], "Csection")
            distance_table <- rbind(distance_table, distances)
          }
        }
      }
    }
  }
  
  colnames(distance_table) <- c("Bodysite", "Family", "Day", "Distance", "SampleID", "R_P")
  rownames(distance_table) <- c(1:nrow(distance_table))
  distance_table <- as.data.frame(distance_table[2:nrow(distance_table),])
  distance_table$Day <- as.numeric(as.character(distance_table$Day))
  distance_table$Distance <- as.numeric(as.character(distance_table$Distance))
  distance_table$R_P <- as.character(distance_table$R_P)
  working_dists1 <- distance_table[distance_table$R_P == "Csection",]
  working_dists1 <- working_dists1[complete.cases(working_dists1),]
  working_dists2 <- distance_table[distance_table$R_P == "Vaginal",]
  working_dists2 <- working_dists2[complete.cases(working_dists2),]
  #fm1 <- lme(Distance ~ Day, data=working_dists1, random= ~ 1 | Family)
  #fm2 <- lme(Distance ~ Day, data=working_dists2, random= ~ 1 | Family)
  #r_sq1 <- fm1$coefficients["fixed"][[1]][[2]]
  #pval1 <- data.frame(coef(summary(fm1)))[2,"p.value"]
  #r_sq2 <- fm2$coefficients["fixed"][[1]][[2]]
  #pval2 <- data.frame(coef(summary(fm2)))[2,"p.value"]
  
  for(i in 1:2){
    if(i == 1){
      working_dists <- working_dists1
      alpha <- working_dists$Distance
      subject <- working_dists$Family
      time <- working_dists$Day
      obs <- -mean(sapply(split(1:nrow(working_dists), subject), 
                          function(ixx) if(length(ixx) < 3) 0 else cor.test(alpha[ixx], time[ixx], method='spear')$statistic))
      mc.stats <- -replicate(999,mean(sapply(split(1:nrow(working_dists), subject), 
                                             function(ixx) if(length(ixx) < 3) 0 else cor.test(alpha[ixx], sample(time[ixx]), method='spear')$statistic)))
      pval1 <- mean(c(obs,mc.stats) >= obs)
    } else {
      working_dists <- working_dists2
      alpha <- working_dists$Distance
      subject <- working_dists$Family
      time <- working_dists$Day
      obs <- -mean(sapply(split(1:nrow(working_dists), subject), 
                          function(ixx) if(length(ixx) < 3) 0 else cor.test(alpha[ixx], time[ixx], method='spear')$statistic))
      mc.stats <- -replicate(999,mean(sapply(split(1:nrow(working_dists), subject), 
                                             function(ixx) if(length(ixx) < 3) 0 else cor.test(alpha[ixx], sample(time[ixx]), method='spear')$statistic)))
      pval2 <- mean(c(obs,mc.stats) >= obs)
    }
  }
  
  
  r_p_plot <- ggplot(distance_table, aes(x=Day, y=Distance, color=R_P)) +
    geom_jitter(alpha=0.65, size=2, width = 0.25) +
    geom_smooth(method=lm, se=FALSE, linetype = "dashed", size=0.5) +
    theme_cowplot(font_size = 7) +
    annotate("text", x=10, y=1.25, label= paste("P=", round(pval1, digits=3)), size=2, color= "#c3c823") +
    #annotate("text", x=10, y=1.2, label= paste("R2=", round(r_sq1, digits=3)), size=2, color= "#5F8CA2") +
    annotate("text", x=15, y=1.25, label= paste("P=", round(pval2, digits=3)), size=2, color="#053d58") +
    #annotate("text", x=15, y=1.2, label= paste("R2=", round(r_sq2, digits=3)), size=2, color="#e88d3f") +
    guides(color=FALSE) +
    labs(y="Distance to Mom") +
    scale_color_manual(values= c("#c3c823", "#053d58"))
  
  
  plot_list_r_p[[s]] <- r_p_plot
}

######################################################################
#Plot candida albicans on day 30 in skin
map2 <- mapping[mapping$planned_sampling_day_0_1_3_7_14_30_60_90___365 == 30 & mapping$SuperbodysiteOralSkinNoseVaginaAnalsAureola == "Skin",]
candida <- ggplot(map2, aes(y=Candida_albicans, x=Delivery_Vvaginal_Ccs_IcsInoc, fill=Delivery_Vvaginal_Ccs_IcsInoc)) +
  geom_boxplot(outlier.shape = NA) +
  theme_cowplot(font_size=7) +
  geom_jitter(position=position_jitter(0.1), shape=1, size=1) +
  guides(fill=FALSE) +
  labs(y = "Candida albicans (relative abundance)", x= "Birth Mode") +
  scale_fill_manual(values=c("#c3c823", "#053d58")) +
  annotate("text", x="C", y=0.60, label= paste("P=0.004"), size=2)

###Plot final
side <- plot_grid(dist_mom, candida, nrow=2, labels="auto", label_size=12)
final_plot <- plot_grid(together2, side, ncol=2, rel_widths = c(0.55, 0.45), labels="auto", label_size = 12)
pdf(paste(main_fp, "/Figure3.pdf", sep=""), width=6.64, height=3)
plot(final_plot)
dev.off()


