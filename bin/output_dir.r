####Make output directories####

now<- format(Sys.time(), "%H-%M-%S")
main_fp <- paste("analysis",now, sep="-")
dir.create(main_fp)

dir.create(paste(main_fp, "correlation_network", sep='/'))

dir.create(paste(main_fp, "alpha_div", sep='/'))
dir.create(paste(main_fp, "alpha_div/Baby_VC", sep="/"))
dir.create(paste(main_fp, "alpha_div/Baby_Sites", sep="/"))
dir.create(paste(main_fp, "alpha_div/linear", sep="/"))
dir.create(paste(main_fp, "alpha_div/time", sep="/"))
dir.create(paste(main_fp, "alpha_div/vc_time", sep="/"))

dir.create(paste(main_fp, "beta_div", sep='/'))
dir.create(paste(main_fp, "beta_div/PCOA_bodysite", sep='/'))
dir.create(paste(main_fp, "beta_div/PCOA_birthmode", sep='/'))
dir.create(paste(main_fp, "beta_div/PCOA_baby", sep='/'))
dir.create(paste(main_fp, "beta_div/PCOA_baby_all_cats", sep='/'))
dir.create(paste(main_fp, "beta_div/PCOA_family", sep='/'))
dir.create(paste(main_fp, "beta_div/PCOA_mom", sep='/'))

dir.create(paste(main_fp, "beta_div/distance_baby", sep='/'))
dir.create(paste(main_fp, "beta_div/distance_toMom", sep='/'))

dir.create(paste(main_fp, "taxa_sum", sep='/'))
dir.create(paste(main_fp, "taxa_sum/Baby_EL", sep='/'))
dir.create(paste(main_fp, "taxa_sum/Mom_VA", sep='/'))
dir.create(paste(main_fp, "taxa_sum/Baby_Day_Site", sep='/'))
dir.create(paste(main_fp, "taxa_sum/family", sep='/'))
dir.create(paste(main_fp, "taxa_sum/baby_site", sep='/'))


dir.create(paste(main_fp, "diff_taxa", sep='/'))
dir.create(paste(main_fp, "diff_taxa/Baby_VC", sep='/'))
dir.create(paste(main_fp, "diff_taxa/Body_Sites", sep='/'))
dir.create(paste(main_fp, "diff_taxa/Body_Sites_gen", sep='/'))
dir.create(paste(main_fp, "diff_taxa/Body_Sites_gen/pairwise", sep='/'))
dir.create(paste(main_fp, "diff_taxa/Mom_Sites", sep='/'))
dir.create(paste(main_fp, "diff_taxa/Baby_Mom_Parap", sep='/'))

dir.create(paste(main_fp, "taxa_time", sep='/'))




