#Distance plots
######################################################################
#distance to self over time

output_dir <- paste(main_fp, "beta_div/distance_baby/", sep='/')

timeline_subjects <- as.character(unique(mapping[mapping$planned_sampling_day_0_1_3_7_14_30_60_90___365 == '30',"familyN"]))
mapping$familyN <- as.character(mapping$familyN)

for(t in 1:length(beta_tables)){
  beta_div <- beta_tables[[t]]
  beta_name <- beta_metrics[t]
  distance_table <- t(data.frame(c("bodysite", "subject", "day", "distance", "sample")))
  body_names <- c("Forehead","Oralmucosa","Anal")
  for(b in 1:length(Bodysites_B)){
    body_now <- Bodysites_B[[b]]
    for(s in 1:length(timeline_subjects)){
      baby_now <- rownames(mapping[mapping$familyN == timeline_subjects[s],])
      thirty <- rownames(mapping[mapping$planned_sampling_day_0_1_3_7_14_30_60_90___365 == 30 & mapping$familyN == timeline_subjects[s] & mapping$bodysite_OralMucosa_Forehead_VolarRight_PalmRight_FootRight_Vaginal_Anal_Feces == body_names[b],])
      for(d in 1:(length(Days)-1)){
        day_now <- Days[[d]]
        working <- intersect(day_now, baby_now)
        working <- intersect(working, body_now)
        if(length(working) == 0){
          distances <- c(names(Bodysites_B)[b], timeline_subjects[s], names(Days)[d], NA, NA)
          distance_table <- rbind(distance_table, distances)
        } else {
          distances <- c(names(Bodysites_B)[b], timeline_subjects[s], names(Days)[d], beta_div[working[[1]], thirty], working[[1]])
          distance_table <- rbind(distance_table, distances)
        }
      }
    }
  }
  colnames(distance_table) <- c("Bodysite", "Subject", "Day", "Distance", "SampleID")
  rownames(distance_table) <- c(1:nrow(distance_table))
  distance_table <- as.data.frame(distance_table[2:nrow(distance_table),])
  distance_table$Day <- as.numeric(as.character(distance_table$Day))
  distance_table$Distance <- as.numeric(as.character(distance_table$Distance))
  distance_table$Subject <- as.character(distance_table$Subject)
  distance_table$Bodysite <- as.character(distance_table$Bodysite)
  max_distance <- max(distance_table$Distance, na.rm=T)
  min_distance <- max(distance_table$Distancem, na.rm=T)
  
  fig_list <- c()
  for(x in 1:length(unique(distance_table$Bodysite))){
    site <- unique(distance_table$Bodysite)[x]
    figure <- ggplot(distance_table[distance_table$Bodysite == site,], aes(x=Day, y=Distance)) +
      geom_jitter(alpha=0.65, size=4, color=body_cols[x], width = 0.25) +
      geom_smooth(se=FALSE, color="#99897E") +
      expand_limits(y=c(min_distance,max_distance)) +
      guides(fill=FALSE) +
      labs(y="Distance to Self")
    fig_list[[site]] <- figure
  }
  distance_self <- plot_grid(fig_list[["Skin"]], fig_list[["Oral"]], fig_list[["Anal_B"]], ncol=3)
  plot_name <- paste(beta_name, "_distanceSelf.pdf", sep="")
  save_plot(paste(output_dir, plot_name), distance_self, ncol=3, base_aspect_ratio = 1.8)
}


######################################################################
#Average distance to mom's samples, by birth mode
output_dir <- paste(main_fp, "beta_div/distance_toMom/", sep='/')

for(m in 1:length(Bodysites_M)){
  mom <- Bodysites_M[[m]]
  mom_site <- names(Bodysites_M)[m]
  for(b in 1:length(beta_tables)){
    beta_div <- beta_tables[[b]]
    beta_name <- beta_metrics[b]
    
    dist_plots <- c()
    for(s in 1:length(Bodysites_B)){
      site <- Bodysites_B[[s]]
      v_distances <- t(data.frame(c("SampleID", "Distance_MomV")))
      for(i in 1:length(Families)){
        mom_v <- intersect(Families[[i]], mom)
        if(length(mom_v) > 0){
          babies_working <- intersect(Families[[i]], site)
          if(length(babies_working) > 0){
            for(a in 1:length(babies_working)){
              working_sample <- babies_working[[a]]
              new_row <- c(working_sample, beta_div[working_sample,mom_v])
              v_distances <- rbind(v_distances, new_row)
            }
          }
        }
      }
      colnames(v_distances) <- v_distances[1,]
      v_distances <- merge(v_distances, mapping, by="SampleID")
      v_distances$Distance_MomV <- as.numeric(as.character(v_distances$Distance_MomV))
      print(names(Bodysites_B)[s])
      print(t.test(v_distances$Distance_MomV~v_distances$Delivery_Vvaginal_Ccs_IcsInoc))
      
      #plot distances
      dist_mom <- ggplot(v_distances, aes(x=Delivery_Vvaginal_Ccs_IcsInoc, y=Distance_MomV, fill=Delivery_Vvaginal_Ccs_IcsInoc)) +
        geom_boxplot(outlier.shape = NA) +
        geom_jitter(position=position_jitter(0.1), shape=1, size=3) +
        labs(x="Distance to Mom") +
        scale_fill_manual(values=c("#c3c823", "#053d58")) +
        guides(fill=FALSE) +
        labs(x='', y='Distance to Mom')
      t.test(v_distances[,"Distance_MomV"] ~ v_distances[,"Delivery_Vvaginal_Ccs_IcsInoc"])
      name <- names(Bodysites_B)[s]
      dist_plots[[name]] <- dist_mom
    }
    distances_mom <- plot_grid(dist_plots[["Skin"]], dist_plots[["Oral"]], dist_plots[["Anal_B"]], ncol=3)
    
    file_name <- paste(beta_name, "_", mom_site, "_BirthMode", ".pdf", sep="")
    file_path <- (paste(output_dir, file_name, sep=""))
    save_plot(file_path, distances_mom, ncol=3)
  }
}


######################################################################
#Average distance to mom's sample's vs random mom
output_dir <- paste(main_fp, "beta_div/distance_toMom/", sep='/')

#Use w_unifrac table:
beta_div <- wunifrac

for(s in 1:length(Bodysites_B)){
  site <- Bodysites_B[[s]]
  site_name <- names(Bodysites_B)[s]
  
  paired_a_distances <- c()
  paired_v_distances <- c()
  r_a_distances <- c()
  r_v_distances <- c()
  
  for(i in 1:length(Families)){
    mom_a <- intersect(Families[[i]], Anal_M)
    mom_v <- intersect(Families[[i]], Vagina)
    babies_working <- intersect(Families[[i]], babies)
    babies_working <- intersect(babies_working, site)
    if(i == length(Families)){
      babies2 <- intersect(Families[[1]], babies)
      babies2 <- intersect(babies2, site)
    } else {
      babies2 <- intersect(Families[[i+1]], babies)
      babies2 <- intersect(babies2, site)
    }
    if(length(mom_a) > 0){
      dists1 <- beta_div[babies_working, mom_a]
      paired_a_distances <- c(paired_a_distances, dists1)
      rdists1 <- beta_div[babies2, mom_a]
      r_a_distances <- c(r_a_distances, rdists1)
    }
    if(length(mom_v) > 0){
      dists2 <- beta_div[babies_working, mom_v]
      paired_v_distances <- c(paired_v_distances, dists2)
      rdists2 <- beta_div[babies2, mom_v]
      r_v_distances <- c(r_v_distances, rdists2)
    }
  }
  distance_table <- data.frame(rep(NA, 170))
  distance_table$r_v_distances <- c(r_v_distances, rep(NA, nrow(distance_table)-length(r_v_distances)))
  distance_table$paired_v_distances <- c(paired_v_distances, rep(NA, nrow(distance_table)-length(paired_v_distances)))
  distance_table$r_a_distances <- c(r_a_distances, rep(NA, nrow(distance_table)-length(r_a_distances)))
  distance_table$paired_a_distances <- c(paired_a_distances, rep(NA, nrow(distance_table)-length(paired_a_distances)))
  distance_table <- distance_table[,2:ncol(distance_table)]
  
  distance_table <- melt(distance_table, na.rm = TRUE)
  which_site <- strsplit(as.character(distance_table$variable), "_", fixed=TRUE)
  sites <- c()
  for(i in 1:length(which_site)){
    site <- which_site[[i]][2]
    sites <- c(sites, site)
  }
  distance_table$site <- sites
  colnames(distance_table)[2] <- "Distance"
  
  #run stats and print
  wtest_v <- wilcox.test(r_v_distances, paired_v_distances)
  wtest_a <- wilcox.test(r_a_distances, paired_a_distances)
  sink(paste(output_dir,"paired_tests.txt", sep=''), append = TRUE)
  print(site_name)
  print(wtest_v)
  print(wtest_a)
  sink()
  
  paired_plot <- ggplot(distance_table, aes_string(x="variable", y="Distance", fill="variable")) +
    geom_boxplot(outlier.shape = NA) +
    facet_wrap(facets=~site, scales="free_x", nrow=1) +
    geom_jitter(position=position_jitter(0.1), shape=1, size=3) +
    theme(legend.position = 'bottom') + 
    guides(fill=FALSE) +
    labs(x="", y = "Distance to Mom") +
    scale_x_discrete(labels=c("random", "paired", "random", "paired")) +
    scale_fill_manual(values= c("#5F8CA2", "#E88D3F", "#5F8CA2", "#E88D3F"))
  
  plot_name <- paste(site_name, "_Paired_MomtoBaby.pdf", sep="")
  file_path <- paste(output_dir, plot_name, sep='')
  pdf(file_path, height=4,width=6)
  print(paired_plot)
  dev.off()
}

