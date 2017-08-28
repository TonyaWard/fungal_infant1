## Making Figure 4
######################################################################
#Alpha diversity over time
#Use shannon

# alpha <- alpha[rownames(mapping),]
# mapping$shannon <- alpha$shannon
# mapping$observed_species <- alpha$observed_species
# mapping$simpson <- alpha$simpson
# mapping$PDwholetree <- alpha$PD_whole_tree
alpha_metrics <- c("shannon", "observed_species", "simpson", "PDwholetree")

alpha_oral <- mapping[Bodysites_B[["Oral"]],]
alpha_skin <- mapping[Bodysites_B[["Skin"]],]
alpha_anal <- mapping[Bodysites_B[["Anal_B"]],]
alpha_tissues <- list(alpha_skin, alpha_oral, alpha_anal)
names(alpha_tissues) <- c("Skin", "Oral", "Anal")


set.seed(68)

scaleFUN <- function(x) sprintf("%.1f", x)

alpha_metric <- "shannon"
alpha_plots <- c()
min_a <- min(mapping[,alpha_metric])
max_a <- max(mapping[,alpha_metric])
for(t in 1:length(alpha_tissues)){
  working_alpha <- melt(alpha_tissues[t], id.vars = c('SampleID', 'Delivery_Vvaginal_Ccs_IcsInoc', 'planned_sampling_day_0_1_3_7_14_30_60_90___365', 'subject_id'), measure.vars = c(alpha_metric))
  working_alpha$Collection_Day <- as.numeric(working_alpha$planned_sampling_day_0_1_3_7_14_30_60_90___365)
  working_alpha$subject_id <- factor(working_alpha$subject_id)
  colnames(working_alpha)[6] <- alpha_metric
  
  ##linear regression accounting for subject
  #fm2 <- lme(shannon ~ Collection_Day, data=working_alpha, random= ~ 1 | subject_id)
  #fm1 <- lme(shannon ~ Collection_Day, data=working_alpha, random= ~ Collection_Day | subject_id)
  #anova(fm2, fm1) #use fm2 based on anova
  #r_sq <- fm2$coefficients["fixed"][[1]][[2]]
  #pval <- data.frame(coef(summary(fm2)))[2,"p.value"]
  
  ##Permuation based test to see if different from random
  # Using a permutation test and a groupwise average of the test statistic of the spearman correlation, we find clinical relevance but not significance.
  alpha <- working_alpha$shannon
  subject <- working_alpha$subject_id
  time <- working_alpha$Collection_Day
  obs <- -mean(sapply(split(1:nrow(working_alpha), subject), 
                      function(ixx) if(length(ixx) < 3) 0 else cor.test(alpha[ixx], time[ixx], method='spear')$statistic))
  mc.stats <- -replicate(999,mean(sapply(split(1:nrow(working_alpha), subject), 
                                         function(ixx) if(length(ixx) < 3) 0 else cor.test(alpha[ixx], sample(time[ixx]), method='spear')$statistic)))
  pval <- mean(c(obs,mc.stats) >= obs)
  
  figure <- ggplot(working_alpha, aes_string(x="Collection_Day", y=alpha_metric)) +
    geom_jitter(alpha=0.65, size=1, color=body_cols[t], width = 0.25) +
    theme_cowplot(font_size = 7) +
    geom_smooth(method=lm, se=FALSE, color="#99897E", linetype = "dashed", size=0.5)+
    guides(fill=FALSE)+
    annotate("text", x=25, y=3.25, label= paste("P=", round(pval, digits=3)), size=1) +
    #annotate("text", x=25, y=3, label= paste("R2=", round(r_sq, digits=3)), size=2) +
    labs(x="Day", y=alpha_metric )+
    expand_limits(y=c(min_a,max_a)) +
    scale_y_continuous(labels=scaleFUN)
  name <- names(alpha_tissues[t])
  alpha_plots[[name]] <- figure
}
alpha_time <- plot_grid(alpha_plots[[1]], alpha_plots[[2]], alpha_plots[[3]],ncol=3, labels="auto", label_size=12)

pdf(paste(main_fp, "/Figure4.pdf", sep="/"), height=1, width=3.5)
plot(alpha_time)
dev.off()
