#plot PCOAs
beta_tables <- list(bray,unifrac,wunifrac)
beta_metrics <- c("bray_curtis", "unifrac", "weight_unifrac")

######################################################################
#Make PCoAs of body site
pcoa_dir_baby <- paste(main_fp, "beta_div/PCOA_bodysite/", sep='/')
early_babies <- intersect(babies, earlies)
late_babies <- intersect(babies, lates)
early_late <- list(early_babies, late_babies)
names(early_late) <- c("early", "late")

for(b in 1:length(beta_tables)){
  beta_table <- beta_tables[[b]]
  beta_name <- beta_metrics[b]
  plot_list <- c()
  for(t in 1:length(early_late)){
    beta_subset <- beta_table[early_late[[t]],early_late[[t]]]
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
    plot1 <- ggplot(PCOA) +
      geom_point(size = 4, alpha=0.65, aes_string(x = "PC1", y = "PC2", color = "SuperbodysiteOralSkinNoseVaginaAnalsAureola")) + 
      scale_color_manual(values=body_cols2) +
      guides(color=guide_legend(nrow=3)) +
      theme(legend.title=element_blank())
    name <- names(early_late)[t]
    plot_list[[name]] <- plot1
  }
  plot_together <- plot_grid(plot_list[[1]], plot_list[[2]], labels=names(early_late), ncol=2)
  plot_name <- paste(beta_name, "_PCOA_bodysiteEL.pdf", sep="")
  save_plot(paste(pcoa_dir_baby, plot_name), plot_together, ncol=2)
}

for(b in 1:length(beta_tables)){
  beta_table <- beta_tables[[b]]
  beta_name <- beta_metrics[b]
  plot_list <- c()
  beta_subset <- beta_table[babies,babies]
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
  plot1 <- ggplot(PCOA) +
      geom_point(size = 4, alpha=0.65, aes_string(x = "PC1", y = "PC2", color = "SuperbodysiteOralSkinNoseVaginaAnalsAureola")) + 
      scale_color_manual(values=body_cols2) +
      guides(color=guide_legend(nrow=3)) +
      theme(legend.title=element_blank())
  plot_name <- paste(beta_name, "_PCOA_bodysite.pdf", sep="")
  save_plot(paste(pcoa_dir_baby, plot_name), plot1)
}


######################################################################
#Make PCOA for each bodysite, divide by early-late, color by birth mode
pcoa_dir_baby <- paste(main_fp, "beta_div/PCOA_birthmode/", sep='/')

for(b in 1:length(beta_tables)){
  beta_table <- beta_tables[[b]]
  beta_name <- beta_metrics[b]
  plot_list2 <- c()
  for(s in 1:length(Bodysites_B)){
    plot_list <- c()
    body_samples <- Bodysites_B[[s]]
    for(t in 1:length(early_late)){
      time_body <- intersect(body_samples, early_late[[t]])
      beta_subset <- beta_table[time_body,time_body]
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
      plot1 <- ggplot(PCOA) +
        geom_point(size = 4, alpha=0.65, aes_string(x = "PC1", y = "PC2", color = "Delivery_Vvaginal_Ccs_IcsInoc")) + 
        scale_color_manual(values=c("#c3c823", "#053d58")) +
        guides(color=guide_legend(nrow=3)) +
        theme(legend.title=element_blank())
      name <- names(early_late)[t]
      plot_list[[name]] <- plot1
    }
    plot_together <- plot_grid(plot_list[[1]], plot_list[[2]], nrow=2)
    name2 <- names(Bodysites_B)[s]
    plot_list2[[name2]] <- plot_together 
  }
  plot_this <- plot_grid(plot_list2[[1]], plot_list2[[2]], plot_list2[[3]], ncol=3)
  plot_name <- paste(beta_name, "_PCOA_bmodeEL.pdf", sep="")
  save_plot(paste(pcoa_dir_baby, plot_name), plot_this, ncol=3, nrow=2, base_aspect_ratio = 1.2)
}

### Do the same but don't separate by early late
for(b in 1:length(beta_tables)){
  beta_table <- beta_tables[[b]]
  beta_name <- beta_metrics[b]
  plot_list <- c()
  for(s in 1:length(Bodysites_B)){
    body_samples <- Bodysites_B[[s]]
    beta_subset <- beta_table[body_samples,body_samples]
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
    plot1 <- ggplot(PCOA) +
      geom_point(size = 4, alpha=0.65, aes_string(x = "PC1", y = "PC2", color = "Delivery_Vvaginal_Ccs_IcsInoc")) + 
      scale_color_manual(values=c("#c3c823", "#053d58")) +
      guides(color=guide_legend(nrow=3)) +
      theme(legend.title=element_blank())
    name <- names(Bodysites_B)[s]
    plot_list[[name]] <- plot1
  }
  plot_this <- plot_grid(plot_list[[1]], plot_list[[2]], plot_list[[3]], ncol=3)
  plot_name <- paste(beta_name, "_PCOA_bmode.pdf", sep="")
  save_plot(paste(pcoa_dir_baby, plot_name), plot_this, ncol=3, base_aspect_ratio = 1)
}



######################################################################
#Plot bodysites sep, color by family
pcoa_dir_baby <- paste(main_fp, "beta_div/PCOA_baby/", sep='/')

baby_colors <- sample(cols2(length(unique(mapping$subject_id))))
#mapping$subject_id <- factor(mapping$subject_id)
names(baby_colors) <- unique(mapping$subject_id)

for(b in 1:length(beta_tables)){
  beta_div <- beta_tables[[b]]
  beta_name <- beta_metrics[b]
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
    
    plot1 <- ggplot(PCOA) +
      geom_point(size = 4, alpha=0.65, aes_string(x = "PC1", y = "PC2", color = "subject_id")) + 
      scale_color_manual(values=baby_colors) +
      #guides(color=guide_legend(nrow=3)) +
      guides(color=F)
      theme(legend.title=element_blank())
    name <- names(Bodysites_B)[a]
    plot_list[[name]] <- plot1
  }
  plots <- plot_grid(plot_list[[1]], plot_list[[2]], plot_list[[3]], labels=names(plot_list),ncol=3)
  plot_name <- paste(beta_name, "_Subject_Bodysites.pdf", sep="")
  fp <- paste(pcoa_dir_baby, plot_name, sep="/")
  save_plot(fp, plots, ncol=3,base_aspect_ratio = 1.2)
}


#Now color by time
mapping$planned_sampling_day_0_1_3_7_14_30_60_90___365 <- as.numeric(mapping$planned_sampling_day_0_1_3_7_14_30_60_90___365)
time_colors <- cols_grad(length(unique(mapping$planned_sampling_day_0_1_3_7_14_30_60_90___365)))
names(time_colors) <- sort(unique(mapping$planned_sampling_day_0_1_3_7_14_30_60_90___365))

for(b in 1:length(beta_tables)){
  beta_div <- beta_tables[[b]]
  beta_name <- beta_metrics[b]
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
    PCOA$planned_sampling_day_0_1_3_7_14_30_60_90___365 <- factor(PCOA$planned_sampling_day_0_1_3_7_14_30_60_90___365, levels= c("1", "3", "7", "14", "21", "30"))
    PCOA$PC1 <- as.numeric(levels(PCOA$PC1))[PCOA$PC1]
    PCOA$PC2 <- as.numeric(levels(PCOA$PC2))[PCOA$PC2]
    PCOA$PC3 <- as.numeric(levels(PCOA$PC3))[PCOA$PC3]
    PCOA$PC4 <- as.numeric(levels(PCOA$PC4))[PCOA$PC4]
    
    plot1 <- ggplot(PCOA) +
      geom_point(size = 4, alpha=0.65, aes_string(x = "PC1", y = "PC2", color = "planned_sampling_day_0_1_3_7_14_30_60_90___365")) + 
      scale_color_manual(values=time_colors) +
      #guides(color=F) +
      theme(legend.title=element_blank())
    name <- names(Bodysites_B)[a]
    plot_list[[name]] <- plot1
  }
  
  plots <- plot_grid(plot_list[[1]], plot_list[[2]], plot_list[[3]], labels=names(plot_list),ncol=3)
  plot_name <- paste(beta_name, "_Day.pdf", sep="")
  fp <- paste(pcoa_dir_baby, plot_name, sep="/")
  save_plot(fp, plots, ncol=3,base_aspect_ratio =1.2)
}

######################################################################
#Loop through the taxa, plot those correlated with PC1 or PC2
#just use wunifrac.. you can change that if you want
pcoa_dir_baby <- paste(main_fp, "beta_div/PCOA_baby_all_cats/", sep='/')

#for(b in 1:length(beta_tables)){
#  beta_div <- beta_tables[[b]]
#  beta_name <- beta_metrics[b]
  beta_div <- wunifrac
  beta_subset <- beta_div[babies,babies]
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
  for(i in 1:nrow(taxa_table)){
    header <- rownames(taxa_table)[[i]]
    color_by <- gsub(" ", "_", header)
    print(header)
    test_cor1 <- cor.test(PCOA[,"PC1"], as.numeric(PCOA[,color_by]), method="spearman")
    test_cor2 <- cor.test(PCOA[,"PC2"], as.numeric(PCOA[,color_by]), method="spearman")
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
      fp <- paste(pcoa_dir_baby, name1, sep="")
      pdf(fp, height=4,width=6)
      print(plot1)
      dev.off()
    }
    if(! is.na(test_cor2$p.value) & test_cor2$p.value < 0.05){
      plot2 <- ggplot(PCOA) +
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
      name2 <- paste(color_by, "2", ".pdf", sep="")
      fp <- paste(pcoa_dir_baby, name2, sep="")
      pdf(fp, height=4,width=6)
      print(plot1)
      dev.off()
    }
  }
#}
#print color range:
cols_grad2
fp <- paste(pcoa_dir_baby, "color_scale.pdf", sep="")
pdf(fp, height=2, width=4, useDingbats = F)
plot(rep(0,10),col=cols_grad2(10),pch=19,cex=3)
dev.off()
######################################################################
##Plot each bodysite separate, with mom_v added in. Baby=circle, mom=triangle

pcoa_dir_baby <- paste(main_fp, "beta_div/PCOA_family/", sep='/')

#Make PCOa joined to mom - one set of plots for V one for A
for(m in 1:length(Bodysites_M)){
  mom <- Bodysites_M[[m]]
  mom_site <- names(Bodysites_M)[m]
  for(b in 1:length(beta_tables)){
    beta_div <- beta_tables[[b]]
    beta_name <- beta_metrics[b]
    plots<-c()
    plots2 <- c()
    for(i in 1:length(Bodysites_B)){
      keeps <- c()
      families_plotted <- c()
      babies_site <- Bodysites_B[[i]]
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
      
      plot1 <- ggplot() + #plot with lines
        geom_point(data=PCOA,size = 4, alpha=0.65, aes_string(x = "PC1", y = "PC2", color = "familyN", shape = "motherorbaby_M_B")) + 
        scale_color_manual(values=fam_colors) +
        guides(color=guide_legend(nrow=3))+
        geom_segment(data=connections,aes(x=x1, y=y1, xend=x2, yend=y2, color=familyN))+
        guides(color=FALSE)+
        guides(shape=FALSE) +
        theme(legend.title = element_blank())
      
      plot2 <- ggplot() + #plot with no lines
        geom_point(data=PCOA,size = 4, alpha=0.65, aes_string(x = "PC1", y = "PC2", color = "familyN", shape = "motherorbaby_M_B")) + 
        scale_color_manual(values=fam_colors) +
        guides(color=FALSE)+
        guides(shape=FALSE) +
        theme(legend.title = element_blank())
      
      name <- names(Bodysites_B[i])
      plots[[name]] <- plot1
      plots2[[name]] <- plot2
    }
    joined_pcoa <- plot_grid(plots[["Skin"]], plots[["Oral"]], plots[["Anal_B"]], ncol=3)
    file_name <- paste(beta_name, "_", mom_site, ".pdf", sep="")
    file_path <- (paste(pcoa_dir_baby, file_name, sep=""))
    save_plot(file_path, joined_pcoa, ncol=3)
    notjoined_pcoa <- plot_grid(plots2[["Skin"]], plots2[["Oral"]], plots2[["Anal_B"]], ncol=3)
    file_name <- paste(beta_name, "_NOTJOINED", mom_site, ".pdf", sep="")
    file_path <- (paste(pcoa_dir_baby, file_name, sep=""))
    save_plot(file_path, notjoined_pcoa, ncol=3)
  }
}


######################################################################
#PCOa mom samples
pcoa_dir_mom <- paste(main_fp, "beta_div/PCOA_mom/", sep='/')
moms_body <- c(Bodysites_M[[1]], Bodysites_M[[2]])

for(b in 1:length(beta_tables)){
  beta_div <- beta_tables[[b]]
  beta_name <- beta_metrics[b]
  beta_table <-  beta_div[moms_body, moms_body]
  PCOA <- pcoa(beta_table)$vectors
  for(c in 1:ncol(PCOA)){
    colnames(PCOA)[c] <- paste("PC",c, sep="")
  }
  PCOA <- cbind(PCOA, rownames(PCOA))
  colnames(PCOA)[ncol(PCOA)] <- "SampleID"
  PCOA <- merge(PCOA, mapping2, by="SampleID")
  PCOA$PC1 <- as.numeric(levels(PCOA$PC1))[PCOA$PC1]
  PCOA$PC2 <- as.numeric(levels(PCOA$PC2))[PCOA$PC2]
  
  momPCOa <- ggplot(PCOA) +
    geom_point(size = 4, alpha=0.65, aes(x = PC1, y = PC2, color = SuperbodysiteOralSkinNoseVaginaAnalsAureola)) + 
    scale_color_manual(values=c("#b75f6d","#99897e")) +
    theme(legend.title = element_blank())+
    guides(color=guide_legend(nrow=2))
    #guides(color=FALSE)
  file_name <- paste(beta_name, "_momVA.pdf", sep="")
  file_path <- paste(pcoa_dir_mom, file_name, sep="")
  save_plot(file_path, momPCOa, base_aspect_ratio = 1.2)
}

######################################################################

#Plot just mom's v or A samples, look at impact of birth mode (abx)

for(m in 1:length(Bodysites_M)){
  mom <- Bodysites_M[[m]]
  mom_site <- names(Bodysites_M)[m]
  for(b in 1:length(beta_tables)){
    beta_div <- beta_tables[[b]]
    beta_name <- beta_metrics[b]
    beta_table <- beta_div[mom,mom]
    PCOA <- pcoa(beta_table)$vectors
    for(c in 1:ncol(PCOA)){
      colnames(PCOA)[c] <- paste("PC",c, sep="")
    }
    PCOA <- cbind(PCOA, rownames(PCOA))
    colnames(PCOA)[ncol(PCOA)] <- "SampleID"
    PCOA <- merge(PCOA, mapping2, by="SampleID")
    PCOA$PC1 <- as.numeric(levels(PCOA$PC1))[PCOA$PC1]
    PCOA$PC2 <- as.numeric(levels(PCOA$PC2))[PCOA$PC2]
    
    momPCOa <- ggplot(PCOA) +
      geom_point(size = 4, alpha=0.65, aes(x = PC1, y = PC2, color = Delivery_Vvaginal_Ccs_IcsInoc)) + 
      scale_color_manual(values=c("#c3c823", "#053d58")) +
      theme(legend.title = element_blank())+
      guides(color=guide_legend(nrow=2))
    #guides(color=FALSE)
    file_name <- paste(mom_site,beta_name, "_momBirthMode.pdf", sep="")
    file_path <- paste(pcoa_dir_mom, file_name, sep="")
    save_plot(file_path, momPCOa, base_aspect_ratio = 1.2)
  }
}
