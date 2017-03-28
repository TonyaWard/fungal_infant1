#Add alpha to mapping file
alpha <- alpha[rownames(mapping),]
mapping$shannon <- alpha$shannon
mapping$observed_species <- alpha$observed_species
mapping$simpson <- alpha$simpson
mapping$PDwholetree <- alpha$PD_whole_tree
alpha_metrics <- c("shannon", "observed_species", "simpson", "PDwholetree")

alpha_oral <- mapping[Bodysites_B[["Oral"]],]
alpha_skin <- mapping[Bodysites_B[["Skin"]],]
alpha_anal <- mapping[Bodysites_B[["Anal_B"]],]
alpha_tissues <- list(alpha_skin, alpha_oral, alpha_anal)
names(alpha_tissues) <- c("Skin", "Oral", "Anal")


######################################################################
#Compare bodysite alpha diversity for each time point, use mom and babies

alpha_table_baby <- rbind(alpha_skin, alpha_oral, alpha_anal)

alpha_analM <- mapping[Bodysites_M[["Anal_M"]],]
alpha_vagina <- mapping[Bodysites_M[["Vagina"]],]
alpha_table_mom <- rbind(alpha_analM, alpha_vagina)
alpha_table_mom <- alpha_table_mom[alpha_table_mom$planned_sampling_day_0_1_3_7_14_30_60_90___365 == 1,]

alpha_table_baby$planned_sampling_day_0_1_3_7_14_30_60_90___365 <- as.numeric(alpha_table_baby$planned_sampling_day_0_1_3_7_14_30_60_90___365)
alpha_table_baby$SuperbodysiteOralSkinNoseVaginaAnalsAureola <- factor(alpha_table_baby$SuperbodysiteOralSkinNoseVaginaAnalsAureola, levels=c("Skin", "Oral", "Anal"))

for(i in 1:length(alpha_metrics)){
  alpha_use <- alpha_metrics[i]
  merged <- rbind(alpha_table_mom, alpha_table_baby)
  max_val <- max(merged[,alpha_use])
  min_val <- min(merged[,alpha_use])
  alpha_plot1 <- ggplot(alpha_table_baby, aes_string(x="SuperbodysiteOralSkinNoseVaginaAnalsAureola", y=alpha_use, fill="SuperbodysiteOralSkinNoseVaginaAnalsAureola"), color="black") +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(position=position_jitter(0.1), shape=1, size=3) +
    facet_wrap(facets=~planned_sampling_day_0_1_3_7_14_30_60_90___365, scales="free_x", nrow=1) +
    labs(x="Body Site") +
    expand_limits(y=c(min_val,max_val)) +
    guides(fill=F) +
    scale_fill_manual(values=body_cols)
  
  alpha_plot2 <- ggplot(alpha_table_mom, aes_string(x="SuperbodysiteOralSkinNoseVaginaAnalsAureola", y=alpha_use, fill="SuperbodysiteOralSkinNoseVaginaAnalsAureola")) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(position=position_jitter(0.1), shape=1, size=3) +
    facet_wrap(facets=~planned_sampling_day_0_1_3_7_14_30_60_90___365, scales="free_x", nrow=1) +
    theme(legend.position = 'bottom') + 
    labs(x="Body Site") +
    expand_limits(y=c(min_val,max_val)) +
    guides(fill=F) +
    scale_fill_manual(values=c("#b75f6d", "#99897E"))
  
  plot_together <- plot_grid(alpha_plot2, alpha_plot1, ncol=2,rel_widths=c(1,4))
  alpha_name <- paste("alpha_bodysites_w_mom_", alpha_use, ".pdf", sep='')
  plot_this <- paste(main_fp, "alpha_div/Baby_Sites", alpha_name, sep='/')
  save_plot(plot_this, plot_together,
            ncol = 2,
            nrow = 1,
            base_aspect_ratio = 1.8)
}

#Run stats to see if any are different
alpha_dir <- paste(main_fp, "alpha_div/Baby_Sites/", sep='/')
file_name <- paste(alpha_dir, "Bodysites_Alpha_Stats_FDRadjust.txt", sep='')
sink(file_name)
sink()

####FINISH THIS!!!
# for(a in 1:length(alpha_metrics)){
#   alpha_metric <- alpha_metrics[a]
#   pvals <- c()
#   working_alphas <- data.frame(rep(1:17, 6), c(rep(1, 17), rep(3,17), rep(7, 17), rep(14, 17), rep(21, 17), rep(30,17)), rep(NA, 102), rep(NA, 102), rep(NA, 102))
#   colnames(working_alphas) <- c("subject_id", "day", "skin", "oral", "anal")
#   
#   
#   
#   for(d in 1:length(Days)){
#     working_day <- Days[[d]]
#     name_day <- names(Days[d])
#     for(s in 1:length(unique(alpha_table_baby$subject_id))){
#       working_babe <- unique(alpha_table_baby$subject_id)[s]
#       working_table <- alpha_table_baby[alpha_table_baby$subject_id == working_babe,]
#       working_table <- working_table[intersect(rownames(working_table), working_day),]
#       alphas <- c(working_table[working_table$SuperbodysiteOralSkinNoseVaginaAnalsAureola == "Skin", alpha_metric], working_table[working_table$SuperbodysiteOralSkinNoseVaginaAnalsAureola == "Oral", alpha_metric], working_table[working_table$SuperbodysiteOralSkinNoseVaginaAnalsAureola == "Anal", alpha_metric], working_babe)
#       working_alphas <- cbind(working_alphas,alphas)
#     }
#     #Run Tests for each day. 
#   }
#     test_p <- t.test(working_table[,alpha_metric]~working_table[,"SuperbodysiteOralSkinNoseVaginaAnalsAureola"])$p.value
#     pvals <- c(pvals, test_p)
#     names(pvals)[length(pvals)] <- paste(name_day, name1, name2, sep="_")
# 
#   adjusted_pvals <- p.adjust(pvals, method="fdr")
#   sink(file_name, append=TRUE)
#   print(alpha_metric)
#   print(adjusted_pvals)
#   sink()
# }


for(a in 1:length(alpha_metrics)){
  alpha_metric <- alpha_metrics[a]
  pvals <- c()
  for(d in 1:length(Days)){
    working_day <- Days[[d]]
    name_day <- names(Days[d])
    for(t in 1:(length(Bodysites_B)-1)){
      for(x in (t+1):length(Bodysites_B)){
        test1 <- intersect(Bodysites_B[[t]],working_day)
        name1 <- names(Bodysites_B)[t]
        test2 <- intersect(Bodysites_B[[x]],working_day)
        name2 <- names(Bodysites_B)[x]
        working_table <- mapping[union(test1,test2),]
        test_p <- t.test(working_table[,alpha_metric]~working_table[,"SuperbodysiteOralSkinNoseVaginaAnalsAureola"])$p.value
        pvals <- c(pvals, test_p)
        names(pvals)[length(pvals)] <- paste(name_day, name1, name2, sep="_")
      }
    }
  }
  adjusted_pvals <- p.adjust(pvals, method="fdr")
  sink(file_name, append=TRUE)
  print(alpha_metric)
  print(adjusted_pvals)
  sink()
}

#Test mom's alpha
file_name <- paste(alpha_dir, "Mom_Alpha.txt", sep='')
for(a in 1:length(alpha_metrics)){
  alpha_metric <- alpha_metrics[a]
  working_table <- mapping[union(Anal_M,Vagina),]
  test_p <- t.test(working_table[,alpha_metric]~working_table[,"SuperbodysiteOralSkinNoseVaginaAnalsAureola"])$p.value
  sink(file_name, append=TRUE)
  print(alpha_metric)
  print(test_p)
  sink()
}


######################################################################
#Plot alpha diversity over time, one plot per bodysite, color by birth mode

for(a in 1:length(alpha_metrics)){
  alpha_plots  <- c()
  alpha_use <- alpha_metrics[a]
  for(t in 1:length(alpha_tissues)){
    working_alpha <- melt(alpha_tissues[t], id.vars = c('SampleID', 'Delivery_Vvaginal_Ccs_IcsInoc', 'planned_sampling_day_0_1_3_7_14_30_60_90___365'), measure.vars = c(alpha_use))
    working_alpha$Collection_Day <- as.numeric(as.character(working_alpha$planned_sampling_day_0_1_3_7_14_30_60_90___365))
    min_a <- min(alpha_table_baby[,alpha_use])
    max_a <- max(alpha_table_baby[,alpha_use])
    figure <- ggplot(working_alpha, aes(x=Collection_Day, y=value, color=Delivery_Vvaginal_Ccs_IcsInoc, group=Delivery_Vvaginal_Ccs_IcsInoc)) +
      geom_jitter(width = 0.25, size=4, alpha=0.65) +
      geom_smooth(se=FALSE) +
      scale_color_manual(values=c("#c3c823", "#053d58")) +
      theme(legend.title=element_blank()) +
      labs(x="Day", y=alpha_metric ) +
      expand_limits(y=c(min_a,max_a))
    name <- names(alpha_tissues[t])
    alpha_plots[[name]] <- figure
  }
  plot3 <- plot_grid(alpha_plots[[1]], alpha_plots[[2]], alpha_plots[[3]],
                     labels=c(names(alpha_plots)),hjust=-2, ncol=3)
  alpha_name <- paste("alpha_time_", alpha_use, ".pdf", sep='')
  plot_this <- paste(main_fp, "alpha_div/vc_time", alpha_name, sep='/')
  save_plot(plot_this, plot3,
            ncol = 3,
            nrow = 1,
            base_aspect_ratio = 1.8)
}

######################################################################
#Plot alpha diversity over time, one plot per bodysite
##alpha over time
for(a in 1:length(alpha_metrics)){
  alpha_metric <- alpha_metrics[a]
  alpha_plots <- c()
  min_a <- min(alpha_table_baby[,alpha_metric])
  max_a <- max(alpha_table_baby[,alpha_metric])
  for(t in 1:length(alpha_tissues)){
    working_alpha <- melt(alpha_tissues[t], id.vars = c('SampleID', 'Delivery_Vvaginal_Ccs_IcsInoc', 'planned_sampling_day_0_1_3_7_14_30_60_90___365', 'subject_id'), measure.vars = c(alpha_metric))
    working_alpha$Collection_Day <- as.numeric(working_alpha$planned_sampling_day_0_1_3_7_14_30_60_90___365)
    working_alpha$subject_id <- factor(working_alpha$subject_id)
    colnames(working_alpha)[6] <- alpha_metric
    figure <- ggplot(working_alpha, aes_string(x="Collection_Day", y=alpha_metric)) +
      geom_jitter(alpha=0.65, size=4, color=body_cols[t], width = 0.25) +
      geom_smooth(se=FALSE, color="#99897E")+
      guides(fill=FALSE)+
      labs(x="Day", y=alpha_metric )+
      expand_limits(y=c(min_a,max_a))
    name <- names(alpha_tissues[t])
    alpha_plots[[name]] <- figure
  }
  alpha_time <- plot_grid(alpha_plots[[1]], alpha_plots[[2]], alpha_plots[[3]],ncol=3)
  alpha_name <- paste("alpha_time_", alpha_metric, ".pdf", sep='')
  plot_this <- paste(main_fp, "alpha_div/time", alpha_name, sep='/')
  save_plot(plot_this, alpha_time,
            ncol = 3,
            nrow = 1,
            base_aspect_ratio = 1.8)
}


######################################################################
#Test for change in alpha over time by body site
#Print only those significant
for(i in 1:length(Bodysites_B)){
  for(j in 1:ncol(alpha)){
    color <- body_cols[i]
    site <- Bodysites_B[[i]]
    site_name <- names(Bodysites_B[i])
    metric_name <- colnames(alpha)[j]
    lmfit <- lm(alpha[site, metric_name] ~ as.numeric(mapping[site,"planned_sampling_day_0_1_3_7_14_30_60_90___365"]))
    cor.tests <- cor.test(alpha[site, metric_name], as.numeric(mapping[site,"planned_sampling_day_0_1_3_7_14_30_60_90___365"]), method="spearman")
    stat_name <- paste(paste(main_fp, "alpha_div/linear/", sep="/"), site_name, metric_name, ".txt")
    if(cor.tests$p.value < 0.05){
      sink(stat_name)
      cat("lm values for:")
      print(site_name)
      print(metric_name)
      print(lmfit)
      print(cor.tests)
      sink()
      pdf_name <- paste(paste(main_fp, "alpha_div/linear/", sep="/"), site_name, metric_name, ".pdf")
      pdf(pdf_name, height=4,width=6)
      plot(alpha[site, metric_name] ~ as.numeric(mapping[site,"planned_sampling_day_0_1_3_7_14_30_60_90___365"]),
           main= site_name,
           xlab="Day",
           ylab= metric_name,
           col=color, pch=16, cex=1.5)
      abline(lmfit, lty=2)
      dev.off()
    }
  }
}

#Test for change in alpha over time by body site, according to birth mode
#Print only those significant
for(i in 1:length(Bodysites_B)){
  for(j in 1:ncol(alpha)){
    color <- body_cols[i]
    site <- Bodysites_B[[i]]
    site <- intersect(site, Vaginal)
    site_name <- paste(names(Bodysites_B[i]), "vaginal", sep="_")
    metric_name <- colnames(alpha)[j]
    lmfit <- lm(alpha[site, metric_name] ~ as.numeric(mapping[site,"planned_sampling_day_0_1_3_7_14_30_60_90___365"]))
    cor.tests <- cor.test(alpha[site, metric_name], as.numeric(mapping[site,"planned_sampling_day_0_1_3_7_14_30_60_90___365"]), method="spearman")
    stat_name <- paste(paste(main_fp, "alpha_div/linear/", sep="/"), site_name, metric_name, ".txt")
    if(cor.tests$p.value < 0.05){
      sink(stat_name)
      cat("lm values for:")
      print(site_name)
      print(metric_name)
      print(lmfit)
      print(cor.tests)
      sink()
      pdf_name <- paste(paste(main_fp, "alpha_div/linear/", sep="/"), site_name, metric_name, ".pdf")
      pdf(pdf_name, height=4,width=6)
      plot(alpha[site, metric_name] ~ as.numeric(mapping[site,"planned_sampling_day_0_1_3_7_14_30_60_90___365"]),
           main= site_name,
           xlab="Day",
           ylab= metric_name,
           col=color, pch=16, cex=1.5)
      abline(lmfit, lty=2)
      dev.off()
    }
  }
}

for(i in 1:length(Bodysites_B)){
  for(j in 1:ncol(alpha)){
    color <- body_cols[i]
    site <- Bodysites_B[[i]]
    site <- intersect(site, Csection)
    site_name <- paste(names(Bodysites_B[i]), "Csection", sep="_")
    metric_name <- colnames(alpha)[j]
    lmfit <- lm(alpha[site, metric_name] ~ as.numeric(mapping[site,"planned_sampling_day_0_1_3_7_14_30_60_90___365"]))
    cor.tests <- cor.test(alpha[site, metric_name], as.numeric(mapping[site,"planned_sampling_day_0_1_3_7_14_30_60_90___365"]), method="spearman")
    stat_name <- paste(paste(main_fp, "alpha_div/linear/", sep="/"), site_name, metric_name, ".txt")
    if(cor.tests$p.value < 0.05){
      sink(stat_name)
      cat("lm values for:")
      print(site_name)
      print(metric_name)
      print(lmfit)
      print(cor.tests)
      sink()
      pdf_name <- paste(paste(main_fp, "alpha_div/linear/", sep="/"), site_name, metric_name, ".pdf")
      pdf(pdf_name, height=4,width=6)
      plot(alpha[site, metric_name] ~ as.numeric(mapping[site,"planned_sampling_day_0_1_3_7_14_30_60_90___365"]),
           main= site_name,
           xlab="Day",
           ylab= metric_name,
           col=color, pch=16, cex=1.5)
      abline(lmfit, lty=2)
      dev.off()
    }
  }
}



######################################################################
#### Testing Babies by Birth Mode, Controling for Bodysite####
#set output directory
alpha_dir <- paste(main_fp, "alpha_div/Baby_VC/", sep='/')

for(i in 1:length(alpha_metrics)){
  alpha_use <- alpha_metrics[i]
  plot_list <- c()
  for(b in 1:length(Bodysites_B)){
    alpha_working <- alpha_table_baby[Bodysites_B[[b]],]
    alpha_plot1 <- ggplot(alpha_working, aes_string(x="Delivery_Vvaginal_Ccs_IcsInoc", y=alpha_use, fill="Delivery_Vvaginal_Ccs_IcsInoc")) +
      geom_boxplot(outlier.shape = NA) +
      geom_jitter(position=position_jitter(0.1), shape=1, size=3) +
      facet_wrap(facets=~planned_sampling_day_0_1_3_7_14_30_60_90___365, scales="free_x", nrow=1) +
      labs(x="Birth Mode", title=names(Bodysites_B[b])) +
      expand_limits(y=c(min_val,max_val)) +
      guides(fill=F) +
      scale_fill_manual(values=c("#c3c823", "#053d58"))
    plot_list[[b]] <- alpha_plot1
  }
  plot_together <- plot_grid(plot_list[[1]], plot_list[[2]], plot_list[[3]], ncol=1, nrow=3)
  pdf(paste(alpha_dir, alpha_use, ".pdf", sep=""), width=7, height=7, useDingbats=FALSE)
  print(plot_together)
  dev.off()
}

#Run stats to see if any are different
file_name <- paste(alpha_dir, "Alpha_Stats_VC_FDRadjust.txt", sep='')
sink(file_name)
sink()

for(a in 1:length(alpha_metrics)){
  alpha_metric <- alpha_metrics[a]
  pvals <- c()
  for(d in 1:length(Days)){
    working_day <- Days[[d]]
    name_day <- names(Days[d])
    for(b in 1:length(Bodysites_B)){
      test1 <- intersect(Bodysites_B[[b]],working_day)
      test1 <- intersect(test1, Vaginal)
      name1 <- paste(names(Bodysites_B)[b], "vaginal", sep="_")
      test2 <- intersect(Bodysites_B[[b]],working_day)
      test2 <- intersect(test2, Csection)
      name2 <- paste(names(Bodysites_B)[b], "Csection", sep="_")
      working_table <- mapping[union(test1,test2),]
      if(length(test1) > 2 && length(test2) >2){
        test_p <- t.test(working_table[,alpha_metric]~working_table[,"Delivery_Vvaginal_Ccs_IcsInoc"])$p.value
      }
      pvals <- c(pvals, test_p)
      names(pvals)[length(pvals)] <- paste(name_day, name1, name2, sep="_")
    }
  }
  adjusted_pvals <- p.adjust(pvals, method="fdr")
  sink(file_name, append=TRUE)
  print(alpha_metric)
  print(adjusted_pvals)
  sink()
}
