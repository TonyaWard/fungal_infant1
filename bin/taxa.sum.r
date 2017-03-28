#Make Taxa Summaries
source("bin/make_taxa_sums.r")

######################################################################
#Make Early - Late by bodysite plot
taxa_dir <- paste(main_fp, "taxa_sum/Baby_EL/", sep='/')

for(f in 1:length(Families)){
  fam_working <- Families[[f]]
  for(b in 1:length(Bodysites_B)){
    body <- Bodysites_B[[b]]
    for(t in 1:length(times)){
      tp <- times[[t]]
      working_samples <- intersect(fam_working, body)
      working_samples <- intersect(working_samples, tp)
      if(length(working_samples) == 0){
        print("no samples")
      } else {
        subset_otu <- taxa_table[,working_samples,drop=F]
        name <- paste(names(Families)[f],names(Bodysites_B)[b], names(times)[t], sep=".")
        EL_taxa[,name] <- rowSums(subset_otu)
      }
      
    }
  }
}
EL_taxa <- EL_taxa[,2:ncol(EL_taxa)]

#convert to RA
subset_otu <- sweep(EL_taxa,2,colSums(EL_taxa),`/`)
subset_otu <- as.data.frame(t(subset_otu))

#if less than 5%, set to 0 and change name to other
for(i in 1:ncol(subset_otu)){
  for(k in 1:nrow(subset_otu)){
    if(subset_otu[k,i] < 0.05){
      subset_otu[k,i] <- 0
    }
  }
}
#Keep only taxa present (remove 0 counts)
subset_otu <- subset_otu[,colSums(subset_otu) > 0]

subset_otu$Other <- 0
#Set other to be the counts subtracted earlier
for(v in 1:nrow(subset_otu)){
  subset_otu[v,"Other"] <- 1 - rowSums(subset_otu[v,])
}

#Find taxa present at highest numbers
subset_otu_na <- subset_otu[,1:ncol(subset_otu)-1]
subset_otu_na[subset_otu_na == 0] <- NA
max_abund <- colMeans(subset_otu_na, na.rm=TRUE)
names(max_abund) <- colnames(subset_otu_na)
max_abund <- c(max_abund, 0)
names(max_abund)[length(max_abund)] <- "Other"
max_abund <- sort(max_abund, decreasing=TRUE)

#add sample IDs to otu table and melt table
subset_otu <- subset_otu[,names(max_abund)]
subset_otu$SampleID <- rownames(subset_otu)

otu <- melt(subset_otu, id.vars = "SampleID", variable.name = "Taxa", value.name = "Relative_Abundance")

#Splitsample ID name to be three columns
otu$Baby <- NA
otu$BodySite <- NA
otu$Time <- NA
for(r in 1:nrow(otu)){
  otu[r,"Baby"] <- strsplit(otu[r,"SampleID"], "[.]")[[1]][1]
  otu[r,"BodySite"] <- strsplit(otu[r,"SampleID"], "[.]")[[1]][2]
  otu[r,"Time"] <- strsplit(otu[r,"SampleID"], "[.]")[[1]][3]
}

#Assign taxa levels (order in bar plot)
otu$Taxa <- factor(otu$Taxa, levels = c(names(max_abund)), ordered=TRUE)
otu <- otu[order(otu$Taxa),]


##Make plots
skin_plot <- ggplot(otu[otu$BodySite == "Skin",], aes_string(x = "SampleID", y = "Relative_Abundance", fill="Taxa")) + 
  geom_bar(stat="identity", position="fill") +
  facet_wrap(facets=~Time, scales = "free_x") +
  guides(fill=FALSE) +
  scale_fill_manual(name= names(taxa_cols), values= taxa_cols) +
  labs(y="Relative Abundance", x='')+
  theme(axis.text = element_blank(), axis.ticks.x = element_blank())

oral_plot <- ggplot(otu[otu$BodySite == "Oral",], aes_string(x = "SampleID", y = "Relative_Abundance", fill="Taxa")) + 
  geom_bar(stat="identity", position="fill") +
  facet_wrap(facets=~Time, scales = "free_x") +
  guides(fill=FALSE) +
  scale_fill_manual(name= names(taxa_cols), values= taxa_cols)+
  labs(y="Relative Abundance", x='')+
  theme(axis.text = element_blank(), axis.ticks.x = element_blank())

anal_plot <- ggplot(otu[otu$BodySite == "Anal_B",], aes_string(x = "SampleID", y = "Relative_Abundance", fill="Taxa")) + 
  geom_bar(stat="identity", position="fill") +
  facet_wrap(facets=~Time, scales = "free_x") +
  guides(fill=FALSE) +
  scale_fill_manual(name= names(taxa_cols), values= taxa_cols)+
  labs(y="Relative Abundance", x='')+
  theme(axis.text = element_blank(), axis.ticks.x = element_blank())

legend <- ggplot(otu[otu$BodySite == "Skin",], aes_string(x = "SampleID", y = "Relative_Abundance", fill="Taxa")) + 
  geom_bar(stat="identity", position="fill") +
  facet_wrap(facets=~Time, scales = "free_x") +
  theme_bw() +
  guides(fill=guide_legend(ncol=3))  +
  scale_fill_manual(name= names(taxa_cols), values= taxa_cols)+
  theme(legend.title = element_blank())
#Make separate legend
g <- ggplotGrob(legend)$grobs
legend2 <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
#Make plot
taxa_sums <- plot_grid(skin_plot, oral_plot, anal_plot, ncol=3)
#save plots
save_plot(paste(taxa_dir, "babies_EL.pdf", sep=""), taxa_sums, ncol=3)
save_plot(paste(taxa_dir, "babies_EL_guide.pdf", sep=""), legend2, ncol=3)


######################################################################
#Taxa summary of baby samples by bodysite, no time sep
taxa_dir <- paste(main_fp, "taxa_sum/baby_site/", sep='/')

b_bodysite <- c(Skin, Oral, Anal_B)
otu <- make_taxa_sums(taxa_table, b_bodysite)

otu$SuperbodysiteOralSkinNoseVaginaAnalsAureola <- factor(otu$SuperbodysiteOralSkinNoseVaginaAnalsAureola, levels=c("Skin", "Oral", "Anal"))
taxa_plot <- ggplot(otu, aes_string(x = "SampleID", y = "Relative_Abundance", fill="Taxa")) + 
  geom_bar(stat="identity", position="fill") +
  facet_wrap(facets=~SuperbodysiteOralSkinNoseVaginaAnalsAureola, scales = "free_x") +
  guides(fill=FALSE) +
  scale_fill_manual(name= names(taxa_cols), values= taxa_cols)+
  labs(y="Relative Abundance", x='')+
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

legend <- ggplot(otu, aes_string(x = "SampleID", y = "Relative_Abundance", fill="Taxa")) + 
  geom_bar(stat="identity", position="fill") +
  facet_wrap(facets=~SuperbodysiteOralSkinNoseVaginaAnalsAureola, scales = "free_x") +
  theme_bw() +
  guides(fill=guide_legend(ncol=3))  +
  scale_fill_manual(name= names(taxa_cols), values= taxa_cols)+
  theme(legend.title = element_blank())

#Make separate legend
g <- ggplotGrob(legend)$grobs
legend2 <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]

save_plot(paste(taxa_dir, "baby_sites.pdf", sep=""), taxa_plot, base_aspect_ratio = 2.0)
save_plot(paste(taxa_dir, "baby_sites_guide.pdf", sep=""), legend2, ncol=3)

######################################################################
#Taxa summary of mom samples

taxa_dir <- paste(main_fp, "taxa_sum/Mom_VA/", sep='/')
m_bodysite <- c(Vagina, Anal_M)
otu <- make_taxa_sums(taxa_table, m_bodysite)

taxa_plot <- ggplot(otu, aes_string(x = "SampleID", y = "Relative_Abundance", fill="Taxa")) + 
  geom_bar(stat="identity", position="fill") +
  facet_wrap(facets=~SuperbodysiteOralSkinNoseVaginaAnalsAureola, scales = "free_x") +
  guides(fill=FALSE) +
  scale_fill_manual(name= names(taxa_cols), values= taxa_cols)+
  labs(y="Relative Abundance", x='')+
  theme(axis.text = element_blank(), axis.ticks.x = element_blank())

legend <- ggplot(otu, aes_string(x = "SampleID", y = "Relative_Abundance", fill="Taxa")) + 
  geom_bar(stat="identity", position="fill") +
  facet_wrap(facets=~SuperbodysiteOralSkinNoseVaginaAnalsAureola, scales = "free_x") +
  theme_bw() +
  guides(fill=guide_legend(ncol=3))  +
  scale_fill_manual(name= names(taxa_cols), values= taxa_cols)+
  theme(legend.title = element_blank())

#Make separate legend
g <- ggplotGrob(legend)$grobs
legend2 <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]

save_plot(paste(taxa_dir, "mom_VA.pdf", sep=""), taxa_plot)
save_plot(paste(taxa_dir, "mom_VA_guide.pdf", sep=""), legend2, ncol=3)


######################################################################
#facet by day, rows = bodysite
taxa_dir <- paste(main_fp, "taxa_sum/Baby_Day_Site/", sep='/')

otu <- make_taxa_sums(taxa_table, babies)
otu$planned_sampling_day_0_1_3_7_14_30_60_90___365 <- factor(otu$planned_sampling_day_0_1_3_7_14_30_60_90___365, levels= c("1", "3","7", "14", "21", "30"))

plot_title <- "Babies_Bodysites"

taxa_plot_a <- ggplot(otu[otu$SampleID %in% Anal_B,], aes_string(x = "SampleID", y = "Relative_Abundance", fill="Taxa")) + 
  geom_bar(stat="identity", position="fill") +
  facet_wrap(facets=~planned_sampling_day_0_1_3_7_14_30_60_90___365, nrow=1,
             scales = "free_x") +
  guides(fill=FALSE) +
  scale_fill_manual(name= names(taxa_cols), values= taxa_cols)+
  labs(y="Relative Abundance", x='')+
  theme(axis.text = element_blank(), axis.ticks.x = element_blank())
  
taxa_plot_o <- ggplot(otu[otu$SampleID %in% Oral,], aes_string(x = "SampleID", y = "Relative_Abundance", fill="Taxa")) + 
  geom_bar(stat="identity", position="fill") +
  facet_wrap(facets=~planned_sampling_day_0_1_3_7_14_30_60_90___365, nrow=1,
             scales = "free_x") +
  guides(fill=FALSE) +
  scale_fill_manual(name= names(taxa_cols), values= taxa_cols)+
  labs(y="Relative Abundance", x='')+
  theme(axis.text = element_blank(), axis.ticks.x = element_blank())

taxa_plot_s <- ggplot(otu[otu$SampleID %in% Skin,], aes_string(x = "SampleID", y = "Relative_Abundance", fill="Taxa")) + 
  geom_bar(stat="identity", position="fill") +
  facet_wrap(facets=~planned_sampling_day_0_1_3_7_14_30_60_90___365, nrow=1,
             scales = "free_x") +
  guides(fill=FALSE) +
  scale_fill_manual(name= names(taxa_cols), values= taxa_cols)+
  labs(y="Relative Abundance", x='')+
  theme(axis.text = element_blank(), axis.ticks.x = element_blank())

name <- paste("Baby_Bodysites", ".pdf", sep='')
file_path <- paste(taxa_dir, name, sep='')
plot_this <- plot_grid(taxa_plot_s, taxa_plot_o, taxa_plot_a, nrow=3)
save_plot(file_path, plot_this, nrow=2)


######################################################################
#Compare Mom to HER baby

family_dir <- paste(main_fp, "taxa_sum/family/", sep='/')

for(f in 1:length(Families)){
  otu <- make_taxa_sums(taxa_table, Families[[f]])
  mom_v <- otu[otu$SampleID %in% Vagina,]
  mom_a <- otu[otu$SampleID %in% Anal_M,]
  baby_table <- rbind(otu[otu$SampleID %in% Anal_B,], otu[otu$SampleID %in% Skin,], otu[otu$SampleID %in% Oral,])
  baby_table$planned_sampling_day_0_1_3_7_14_30_60_90___365 <- factor(baby_table$planned_sampling_day_0_1_3_7_14_30_60_90___365, levels= c("1", "3","7", "14", "21", "30"))
  
  to_plot <- list(mom_v, mom_a, baby_table)
  names(to_plot) <- c("mom_v", "mom_a", "baby_table")
  
  v_plot <- ggplot(mom_v, aes_string(x = "SampleID", y = "Relative_Abundance", fill="Taxa")) + 
    geom_bar(stat="identity", position="fill") +
    guides(fill=FALSE) +
    scale_fill_manual(name= names(taxa_cols), values= taxa_cols)+
    labs(y="Relative Abundance", x='')+
    theme(axis.text = element_blank(), axis.ticks.x = element_blank())
  
  a_plot <- ggplot(mom_a, aes_string(x = "SampleID", y = "Relative_Abundance", fill="Taxa")) + 
    geom_bar(stat="identity", position="fill") +
    guides(fill=FALSE) +
    scale_fill_manual(name= names(taxa_cols), values= taxa_cols)+
    labs(y="Relative Abundance", x='')+
    theme(axis.text = element_blank(), axis.ticks.x = element_blank())
  
  baby_plot <- ggplot(baby_table, aes_string(x = "planned_sampling_day_0_1_3_7_14_30_60_90___365", y = "Relative_Abundance", fill="Taxa")) + 
    geom_bar(stat="identity", position="fill") +
    facet_wrap(facets=~SuperbodysiteOralSkinNoseVaginaAnalsAureola, nrow=3) +
    facet_wrap(facets=~SuperbodysiteOralSkinNoseVaginaAnalsAureola, scales = "free_x") +
    guides(fill=FALSE) +
    scale_fill_manual(name= names(taxa_cols), values= taxa_cols)+
    labs(y="Relative Abundance", x='')+
    theme(axis.text = element_blank(), axis.ticks.x = element_blank())
  
  name <- paste(names(Families[f]), ".pdf", sep='')
  file_path <- paste(family_dir, name, sep='')
  plot_this <- plot_grid(v_plot, baby_plot, a_plot, ncol=3, nrow=1, rel_widths = c(1,3,1))
  save_plot(file_path, plot_this, ncol=3)
}

####Control taxa summary###

taxa_dir <- paste(main_fp, "taxa_sum/", sep='/')
neg_control <- t(control_taxa)
neg_control <- neg_control[1,, drop=FALSE]
neg_control <- melt(neg_control)
names(neg_control) <- c("SampleID", "Taxa", "Count")
neg_control <- neg_control[neg_control$Count > 0,]
neg_control <- neg_control[order(neg_control$Count, decreasing=TRUE),]
neg_control$Taxa <- factor(neg_control$Taxa, levels=neg_control$Taxa)

con_plot <- ggplot(neg_control, aes_string(x = "Taxa", y = "Count", fill="Taxa")) + 
  geom_bar(stat="identity") +
  guides(fill=guide_legend(ncol=3)) +
  scale_fill_manual(name= names(taxa_cols), values= taxa_cols)+
  labs(y="Normalized Counts", x='')+
  theme(axis.text.x = element_blank(), legend.position = 'bottom', legend.title=element_blank())
  
#assign pdf name for plot
name <- "negative_control.pdf"
file_path <- paste(taxa_dir, name, sep='')
pdf(file_path, height=4,width=6)
print(con_plot)
dev.off()

