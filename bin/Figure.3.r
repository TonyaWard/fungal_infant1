##Making Figure 3
#Distances of mom to baby
#Skin to vagina, anal to anal
#distance to mom or random mom


######################################################################
#PCoA of mom and babies, use Skin and Anal only
#Use weighted unifrac
#Make PCOa joined to mom - one set of plots for V one for A

naming_list <- c("Skin to Vaginal", "Anal to Anal")
mom_sites <- list(Vagina, Anal_M)
baby_sites <- list(Skin, Anal_B)
beta_div <- wunifrac
plots <- c()
for(s in 1:length(mom_sites)){
  mom <- mom_sites[[s]]
  baby_site <- baby_sites[[s]]
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
  
  current_plot <- ggplot() + #plot with lines
    geom_point(data=PCOA,size = 2, alpha=0.65, aes_string(x = "PC1", y = "PC2", color = "familyN", shape = "motherorbaby_M_B")) + 
    scale_color_manual(values=fam_colors) +
    theme_cowplot(font_size=7)+
    geom_segment(data=connections,aes(x=x1, y=y1, xend=x2, yend=y2, color=familyN), alpha=0.65)+
    guides(color=FALSE)+
    guides(shape=FALSE) +
    theme(legend.title = element_blank()) 
  plots[[s]] <- current_plot
}
pcoa_to_mom <- plot_grid(plots[[1]], plots[[2]], ncol=2)

######################################################################
#Average distance to mom's sample's vs random mom
#Use w_unifrac table

beta_div <- wunifrac
plot_list <- c()
for(s in 1:length(mom_sites)){
  site <- baby_sites[[s]]
  mom <- mom_sites[[s]]
  paired_distances <- c()
  r_distances <- c()
  
  for(i in 1:length(Families)){
    mom_d <- intersect(Families[[i]], mom)
    babies_working <- intersect(Families[[i]], babies)
    babies_working <- intersect(babies_working, site)
    if(i == length(Families)){
      babies2 <- intersect(Families[[1]], babies)
      babies2 <- intersect(babies2, site)
    } else {
      babies2 <- intersect(Families[[i+1]], babies)
      babies2 <- intersect(babies2, site)
    }
    if(length(mom_d) > 0){
      dists1 <- beta_div[babies_working, mom_d]
      paired_distances <- c(paired_distances, dists1)
      rdists1 <- beta_div[babies2, mom_d]
      r_distances <- c(r_distances, rdists1)
    }
  }
  distance_table <- data.frame(rep(NA, 170))
  distance_table$r_distances <- c(r_distances, rep(NA, nrow(distance_table)-length(r_distances)))
  distance_table$paired_distances <- c(paired_distances, rep(NA, nrow(distance_table)-length(paired_distances)))
  distance_table <- distance_table[,2:ncol(distance_table)]
  
  distance_table <- melt(distance_table, na.rm = TRUE)
  colnames(distance_table)[2] <- "Distance"
  wilcox.pval <- wilcox.test(distance_table$Distance ~ distance_table$variable)$p.value
  
  paired_plot <- ggplot(distance_table, aes_string(x="variable", y="Distance", fill="variable")) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(position=position_jitter(0.1), shape=1, size=1) +
    theme_cowplot(font_size=7) +
    annotate("text", x="r_distances", y=1.2, label= paste("(wilcox)P=", round(wilcox.pval, digits=3)), size=2) +
    guides(fill=FALSE) +
    labs(x="", y = "Distance to Mom") +
    scale_x_discrete(labels=c("random", "paired")) +
    scale_fill_manual(values= c("#5F8CA2", "#E88D3F"))
  
  plot_list[[s]] <- paired_plot
}
random_paired_plot <- plot_grid(plot_list[[1]], plot_list[[2]], ncol=2)


######################################################################
#Distance to mom over time
#Use weighted unifrac

beta_div <- wunifrac
plot_list_r_p <- c()
distance_table <- t(data.frame(c("bodysite", "subject", "day", "distance", "sample", "r_or_p")))
baby_site_names <- c("skin", "anal")
for(s in 1:length(mom_sites)){
  site <- baby_sites[[s]]
  mom <- mom_sites[[s]]
  for(i in 1:length(Families)){
    mom_d <- intersect(Families[[i]], mom)
    babies_working <- intersect(Families[[i]], babies)
    babies_working <- intersect(babies_working, site)
    if(i == length(Families)){
      babies2 <- intersect(Families[[1]], babies)
      babies2 <- intersect(babies2, site)
    } else {
      babies2 <- intersect(Families[[i+1]], babies)
      babies2 <- intersect(babies2, site)
    }
    if(length(mom_d) > 0){
      for(d in 1:(length(Days)-1)){
        day_now <- Days[[d]]
        working <- intersect(day_now, babies_working)
        working_r <- intersect(day_now, babies2)
        if(length(working) == 0){
          distances <- c(baby_site_names[s], names(Families)[i], names(Days)[d], NA, NA, "paired")
          distance_table <- rbind(distance_table, distances)
        } else {
          distances <- c(baby_site_names[s], names(Families)[i], names(Days)[d], beta_div[working[[1]], mom_d], working[[1]], "paired")
          distance_table <- rbind(distance_table, distances)
        }
        if(length(working_r) == 0){
          distances <- c(baby_site_names[s], names(Families)[i], names(Days)[d], NA, NA, "random")
          distance_table <- rbind(distance_table, distances)
        } else {
          distances <- c(baby_site_names[s], names(Families)[i], names(Days)[d], beta_div[working_r[[1]], mom_d], working_r[[1]], "random")
          distance_table <- rbind(distance_table, distances)
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
  working_dists1 <- distance_table[distance_table$R_P == "random",]
  working_dists1 <- working_dists1[complete.cases(working_dists1),]
  working_dists2 <- distance_table[distance_table$R_P == "paired",]
  working_dists2 <- working_dists2[complete.cases(working_dists2),]
  fm1 <- lme(Distance ~ Day, data=working_dists1, random= ~ 1 | Family)
  fm2 <- lme(Distance ~ Day, data=working_dists2, random= ~ 1 | Family)
  r_sq1 <- fm1$coefficients["fixed"][[1]][[2]]
  pval1 <- data.frame(coef(summary(fm1)))[2,"p.value"]
  r_sq2 <- fm2$coefficients["fixed"][[1]][[2]]
  pval2 <- data.frame(coef(summary(fm2)))[2,"p.value"]
  
  r_p_plot <- ggplot(distance_table, aes(x=Day, y=Distance, color=R_P)) +
    geom_jitter(alpha=0.65, size=2, width = 0.25) +
    geom_smooth(method=lm, se=FALSE, linetype = "dashed", size=0.5) +
    theme_cowplot(font_size = 7) +
    annotate("text", x=10, y=1.25, label= paste("P=", round(pval1, digits=3)), size=2, color= "#5F8CA2") +
    annotate("text", x=10, y=1.2, label= paste("R2=", round(r_sq1, digits=3)), size=2, color= "#5F8CA2") +
    annotate("text", x=15, y=1.25, label= paste("P=", round(pval2, digits=3)), size=2, color="#e88d3f") +
    annotate("text", x=15, y=1.2, label= paste("R2=", round(r_sq2, digits=3)), size=2, color="#e88d3f") +
    guides(color=FALSE) +
    labs(y="Distance to Mom") +
    scale_color_manual(values= c("#E88D3F", "#5F8CA2"))
  
  plot_list_r_p[[s]] <- r_p_plot
}


######################################################################
#Compile Figure 3

together <- plot_grid(plots[[1]], plot_list[[1]], plot_list_r_p[[1]], plots[[2]], plot_list[[2]], plot_list_r_p[[2]], ncol=3, nrow=2, rel_widths = c(0.50, 0.45, 0.75), labels="auto", label_size = 12)

pdf(paste(main_fp, "/Figure3.pdf", sep=""), width=6.69, height =4)
plot(together)
dev.off()
