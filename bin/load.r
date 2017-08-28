library("vegan")
library("biom")
library("RColorBrewer")
library("vegan")
library("ggplot2")
library("reshape2")
library("plyr")
library('beeswarm')
library('ape')
library("grid")
library("gridExtra")
library("cowplot")
library("stringr")
library("nlme")

source('bin/output_dir.r')

#Find the data
data_dir <- "data/"

##Option 1 - DeSeq2 --> Not using this.
##OTU table keeps samples with min 100 counts and CSS normalized with
##QIIME 1.9.1's normalize_table.py (css) 
#otu_fp <- paste(data_dir, "ninTable100_CSS.json.biom", sep='')

##Option 2 - Full OTU table - Must normalize with relative abundance here
##OTU table is not normalized or filtered at all
#dimensions = 324 samples (one is taxonomy) and  232 otus
otu_fp <- paste(data_dir, "ninja_otutable_newCut_meta2.biom", sep='')

######################################################################
##Load OTU table and metadata
#otu table is OTU ID (rows) x sample ID (columns)
#dimensions = 260 samples and  232 otus
otutable <- as.matrix(biom_data(read_biom(otu_fp)))
colnames(otutable) <- gsub(".R1", "", colnames(otutable))
colnames(otutable) <- paste(colnames(otutable), "s", sep="_")

##only for option 2
#metadata is samples (rows) x metadata category (columns) (324 samples)
metadata <- sample_metadata(read_biom(otu_fp))
#write.table(metadata, paste(data_dir, "map.txt", sep=""), quote=FALSE, sep="\t",col.names=NA)
rownames(metadata) <- paste(rownames(metadata), "s", sep="_")
metadata$SampleID <- as.character(rownames(metadata))
#taxonomy is OTU_ID by taxonomy level
taxonomy <- observation_metadata(read_biom(otu_fp))
taxonomy$taxonomy7 <- sapply(strsplit(taxonomy$taxonomy7, split='__', fixed=TRUE), function(x) (x[2]))

#run source tracker (this is done on it's own)
#Rscript sourcetracker_for_qiime.r -i ../data/otu-table-source-tracker.txt -m ../data/map-source-tracker.txt -o st
#load source tracker otu table
#st_otu <- t(read.table("data/Unknown_contributions.txt", sep="\t", header=TRUE, comment='', row.names=1))

# metadata <- read.table(map_fp, sep='\t', header=T, check.names = F, comment='')
# colnames(metadata)[1] <- "sample_id"
# rownames(metadata) <- paste(metadata$sample_id, "s", sep="_")
# metadata$sample_id <- rownames(metadata)

######################################################################
##Filter and normalize - (for option 2 only)
#Filter to keep samples that have a minimum of 50 counts
#keep singletons that occur only in more than one sample #option2:(232 to 142)
otutable1 <- otutable[rowSums(otutable > 0) > 1,]
##keep samples with at least 50 sequence counts (option 2: 324 to 279)
otutable2 <- otutable1[,colSums(otutable1) > 50]

##convert to relative abundance
otutable_RA <- sweep(otutable2,2,colSums(otutable2),`/`) #142 OTUs, 279 samples
otutable_counts <- otutable2 # 142 OTUs, 279 samples

######################################################################
##Make an OTU table output that is compatible with the ghost tree 
#Use this for alpha and beta
# max OTU counts: 618,044    min: 54
otutable_output <- round(otutable_RA * 100000) #unifrac doesn't like decimals!
#Must replace under scores with spaces in OTU IDS, because newick doesn't read them properly
new_names <- rownames(otutable_output)
rownames(otutable_output) <- str_replace_all(new_names,"_"," ")

#Write file to output
sink("data/fungal_RA_wholenum.txt") #option2
# sink("data/fungal_norm_wholenum.txt") #option 1 only
cat("#OTUID")
write.table(otutable_output, 
            sep="\t", #tell R to make is tab-delimited
            quote=F, #tell R not to put quotes
            col.names=NA) #formats the column headers properly
sink()

######################################################################
##Remove positive and negative control
#Store the positive and negative values
controls <- c("Negative_s", "Positive_s")
control_table <- cbind(otutable_RA[,"Negative_s"], otutable_RA[,"Positive_s"])
#control_table <- cbind(otutable[,"Negative_s"], otutable[,"Positive_s"]) #For option 1
colnames(control_table) <- controls
otutable_RA <- otutable_RA[, ! colnames(otutable_RA) %in% controls] #Now 277 samples
otutable_counts <- otutable_counts[, ! colnames(otutable_counts) %in% controls]
# otutable <- otutable[, ! colnames(otutable) %in% controls] #For option 1 

######################################################################
#Filter the metadata to keep only samples in the new otutable
ids_keep <- intersect(rownames(metadata), colnames(otutable_RA))
# ids_keep <- intersect(rownames(metadata), colnames(otutable)) ##For option 2
mapping <- metadata[ids_keep,]

##OPTION 2 ONLY
####Add taxonomy at species level to OTU table####
t_table <- otutable_RA
# # t_table <- otutable #opton 1 only
taxonomy2 <- taxonomy[intersect(rownames(taxonomy), rownames(t_table)),]
rownames(t_table) <- taxonomy2$taxonomy7
# 
#Collapse same taxonomies (from 142 to 123)
#If option 1  - it goes from 232 to 189
taxa_table <- aggregate(t_table[,1:ncol(t_table)], by=list(rownames(t_table)), FUN = sum)
rownames(taxa_table) <- taxa_table[,1]
taxa_table <- taxa_table[,2:ncol(taxa_table)] #Taxa table is 123 otus x 277 samples
mapping$nTaxa <- as.factor(colSums(taxa_table > 0))
mapping$nOTUs <- as.factor(colSums(otutable_counts > 0))
#mapping$nOTUs <- as.factor(colSums(otutable > 0)) #option 1 only
# 
#Collapse taxa by correlation
#move forward with this collapsed taxa table (it replaces taxa_table and otutable_counts with a new versions)
# 
#took out those with a .95 correlation
source('bin/corr.network.r')
#If option 1: new otutable = 190 OTUs, taxa_table = 163 taxa
#If option 2: new otutable: start = 142, end = 110 OTUs, taxa_table: start = 123 taxa, end = 96 taxa
# 
####Add taxonomy for control table####
taxonomy3 <- taxonomy[intersect(rownames(taxonomy), rownames(control_table)),]
rownames(control_table) <- taxonomy3$taxonomy7
control_taxa <- aggregate(control_table[,1:ncol(control_table)], by=list(rownames(control_table)), FUN = sum)
rownames(control_taxa) <- control_taxa[,1]
control_taxa <- control_taxa[,2:ncol(control_taxa)]
control_taxa <- control_taxa[rowSums(control_taxa) > 0,] #There are 18 taxa in these samples, all of which overlap with taxa table

####Remove potential contaminants
#potential_contams <- rownames(control_taxa[control_taxa$Negative_s >2,])
#taxa_table <- taxa_table[! rownames(taxa_table) %in% potential_contams,]

######################################################################
##Generate alpha and beta diveristy and load them
#Change the output otu table to biom from command line
#with qiime/1.9.0
#option 3 (RA)
# biom convert  -i fungal_denovo_RA_wholenum.txt -o fungal_denovo_RA_wholenum.biom --to-json --table-type "OTU table"
#option 2 (RA)
# biom convert -i fungal_RA_wholenum.txt -o fungal_RA_wholenum.biom --to-json --table-type "OTU table"
#option 1 (DeSeq2)
# biom convert -i fungal_norm_wholenum.txt -o fungal_norm_wholenum.biom --to-json --table-type "OTU table"

#option closed ref
#use the get_otus_from_ghost_tree.py script (mod the path to the tree you care about)
#the tree we want to use is the 0116_s/dynamic_100_ghost_tree tree
# python get_otus_from_ghost_tree.py

#### !!! !!!!replace underscores with spaces in file

#transfer the tree, RA OTU table and trim list to MSI
#qiime 1.8.0
#Option 2 (RA)
# filter_otus_from_otu_table.py -i fungal_RA_wholenum.biom -o fungal_trimmed_RA_wholenum.biom -e ghost_tree_tips_underscore_fix.txt --negate_ids_to_exclude
#Option 1 (DeSeq2)
# filter_otus_from_otu_table.py -i fungal_norm_wholenum.biom -o fungal_trimmed_norm_wholenum.biom -e ghost_tree_tips_underscore_fix.txt --negate_ids_to_exclude

#run Beta with trimmed output
# beta_diversity.py -i fungal_denovo_RA_wholenum.biom -o MG_Fungal_Beta_denovo -m bray_curtis
# beta_diversity.py -i fungal_trimmed_RA_wholenum.biom -o MG_Fungal_Beta -m bray_curtis,unifrac,weighted_unifrac -t ghost_tree.nwk
# beta_diversity.py -i fungal_trimmed_norm_wholenum.biom -o MG_Fungal_Beta_norm -m bray_curtis,unifrac,weighted_unifrac -t ghost_tree.nwk

#run alpha with trimmed output
# alpha_diversity.py -i fungal_denovo_RA_wholenum.biom -o MG_Fungal_Alpha_denovo.txt -m shannon,simpson,observed_species
# alpha_diversity.py -i fungal_trimmed_RA_wholenum.biom -o MG_Fungal_Alpha.txt -m shannon,simpson,observed_species,PD_whole_tree -t ghost_tree.nwk
# alpha_diversity.py -i fungal_trimmed_norm_wholenum.biom -o MG_Fungal_Alpha_norm.txt -m shannon,simpson,observed_species,PD_whole_tree -t ghost_tree.nwk

#Transfer these to the data/ directory

# #Option 2 (RA)
alpha <- read.table("data/MG_Fungal_Alpha.txt",
                    sep='\t',
                    header=T,
                    row=1)
bray <- read.table("data/MG_Fungal_Beta/bray_curtis_fungal_trimmed_RA_wholenum.txt",
                   sep="\t",
                   header=T,
                   row=1,
                   check.names=F)
unifrac <- read.table("data/MG_Fungal_Beta/unifrac_fungal_trimmed_RA_wholenum.txt",
                      sep="\t",
                      header=T,
                      row=1,
                      check.names=F)
wunifrac <- read.table("data/MG_Fungal_Beta/weighted_unifrac_fungal_trimmed_RA_wholenum.txt",
                       sep="\t",
                       header=T,
                       row=1,
                       check.names=F)

#Option 1 (DeSeq2)
# alpha <- read.table("data/MG_Fungal_Alpha_norm.txt",
#                     sep='\t',
#                     header=T,
#                     row=1)
# bray <- read.table("data/MG_Fungal_Beta_norm/bray_curtis_fungal_trimmed_norm_wholenum.txt",
#                    sep="\t",
#                    header=T,
#                    row=1,
#                    check.names=F)
# unifrac <- read.table("data/MG_Fungal_Beta_norm/unifrac_fungal_trimmed_norm_wholenum.txt",
#                       sep="\t",
#                       header=T,
#                       row=1,
#                       check.names=F)
# wunifrac <- read.table("data/MG_Fungal_Beta_norm/weighted_unifrac_fungal_trimmed_norm_wholenum.txt",
#                        sep="\t",
#                        header=T,
#                        row=1,
#                        check.names=F)
######################################################################

####Get sample IDs for testing####
#Take out the innoculated babies
Innoc <- rownames(mapping[mapping$Delivery_Vvaginal_Ccs_IcsInoc == "I",])
mapping <- mapping[!rownames(mapping) %in% Innoc,] #258 to 209 samples
otutable_RA<- otutable_RA[,rownames(mapping)]
#otutable <- otutable[,rownames(mapping)]
#taxa_table <- taxa_table[,rownames(mapping)]
alpha <- alpha[rownames(mapping),]
bray <- bray[rownames(mapping), rownames(mapping)]
#unifrac <- unifrac[rownames(mapping), rownames(mapping)]
#wunifrac <- wunifrac[rownames(mapping), rownames(mapping)]

#Store bodysites
Skin <- rownames(mapping[mapping$SuperbodysiteOralSkinNoseVaginaAnalsAureola =="Skin",])
Oral_All <- rownames(mapping[mapping$SuperbodysiteOralSkinNoseVaginaAnalsAureola =="Oral",])
Vagina <- rownames(mapping[mapping$SuperbodysiteOralSkinNoseVaginaAnalsAureola =="Vaginal",])
Anal_B <- intersect(rownames(mapping[mapping$SuperbodysiteOralSkinNoseVaginaAnalsAureola =="Anal",]), rownames(mapping[mapping$motherorbaby_M_B == "B",]))
Anal_M <- intersect(rownames(mapping[mapping$SuperbodysiteOralSkinNoseVaginaAnalsAureola =="Anal",]),rownames(mapping[mapping$motherorbaby_M_B == "M",]))

#Store mom and babies
Mom <- rownames(mapping[mapping$motherorbaby_M_B == "M",])
Infant <- rownames(mapping[mapping$motherorbaby_M_B == "B",])

#Store birth mode
Vaginal <- rownames(mapping[mapping$Delivery_Vvaginal_Ccs_IcsInoc =="V",])
Csection <- rownames(mapping[mapping$Delivery_Vvaginal_Ccs_IcsInoc =="C",])
#Innoc <- rownames(mapping[mapping$Delivery_Vvaginal_Ccs_IcsInoc == "I",])

#Store Day Collected
Days_avail <- unique(mapping$planned_sampling_day_0_1_3_7_14_30_60_90___365)
Days <- list()
for(i in 1:length(Days_avail)){
  working_day <- Days_avail[i]
  Day_samples <- list(rownames(mapping[mapping$planned_sampling_day_0_1_3_7_14_30_60_90___365 == working_day,]))
  Days <- c(Days, Day_samples)
  names(Days)[i] <- working_day
}
day_1 <- Days[[1]]
day_3 <- Days[[2]]

#Keep the "before feed" to be oral samples UNLESS there is no before feed sample
BeforeFeed <- rownames(mapping[mapping$Description == "Mouthswabbeforefeeding",])
AfterFeed <- rownames(mapping[mapping$Description == "Mouthswabafterfeeding",])
Oral <- BeforeFeed
for(i in 1:length(Days)){
  day_samples <- Days[[i]]
  mapsubset <- mapping[rownames(mapping) %in% day_samples,]
  babies_before <- mapsubset[mapsubset$Description == "Mouthswabbeforefeeding",]$subject_id
  babies_after <- mapsubset[mapsubset$Description == "Mouthswabafterfeeding",]$subject_id
  duplicates <- intersect(babies_before, babies_after)
  ndup_samples <- rownames(mapsubset[!mapsubset$subject_id %in% duplicates,])
  Keep_after <- intersect(ndup_samples, AfterFeed)
  Oral <- c(Oral, Keep_after)
}

#Store families
fams <- unique(mapping$familyN)
Families <- list()
for(f in 1:length(fams)){
  Families <- c(Families,  list(rownames(mapping[mapping$familyN == fams[f],])))
  names(Families)[f] <- fams[f]
}

#List Baby and mom body sites and delivery modes
Bodysites_B <- list(Skin, Oral, Anal_B)
names(Bodysites_B) <- c("Skin", "Oral", "Anal_B")
Bodysites_M <- list(Anal_M, Vagina)
names(Bodysites_M) <- c("Anal_M", "Vagina")
#D_Mode <- list(Vaginal, Csection, Innoc)
#names(D_Mode) <- c("Vaginal", "Csection", "Innoc")
D_Mode <- list(Vaginal, Csection)
names(D_Mode) <- c("Vaginal", "Csection")

# #Make the Day 1 timepoint full with Day 3 and Day 1 with no overlap
# day1 <- Days[[1]]
# for(i in 1:length(Bodysites_B)){
#   working_bs <- Bodysites_B[[i]]
#   mapsubset <- mapping[rownames(mapping) %in% working_bs,]
#   babies_day1 <- mapsubset[mapsubset$planned_sampling_day_0_1_3_7_14_30_60_90___365 == "1",]$subject_id
#   babies_day3 <- mapsubset[mapsubset$planned_sampling_day_0_1_3_7_14_30_60_90___365 == "3",]$subject_id
#   duplicates <- intersect(babies_day1, babies_day3)
#   ndup_samples <- rownames(mapsubset[!mapsubset$subject_id %in% duplicates,])
#   Keep_after <- intersect(ndup_samples, Days[[2]])
#   day1 <- c(day1, Keep_after)
# }
# 
# Days <- Days[2:6]
# Days[[1]] <- day1
# names(Days)[1] <- "1"

####Add Taxa quantiles to mapping###
ranges <- c(0, 0.001, 0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)
taxa_table <- taxa_table[,rownames(mapping)]
rownames(taxa_table) <- gsub("-", "_", rownames(taxa_table))
mapping <- cbind(mapping, t(taxa_table))
taxa <- colnames(t(taxa_table))
for(i in 1:nrow(mapping)){
  sample_id <- rownames(mapping)[i]
  for(k in 1:length(taxa)){
    taxon <- taxa[k]
    for(m in 1:13){
      if((mapping[sample_id, taxon] >= ranges[m]) && (mapping[sample_id, taxon] < ranges[m])){
        mapping[sample_id, taxon] <- ranges[m]
        next
      } 
    }
  }
}
colnames(mapping) <- gsub(" ", "_", colnames(mapping))

EL_taxa <- taxa_table[,1,drop=F]
earlies <- c(Days[['1']], Days[['3']], Days[['7']])
lates <- c(Days[['14']], Days[['21']], Days[['30']])
times <- list(earlies, lates)
names(times) <- c("early", "late")

babies <- c(Bodysites_B[[1]], Bodysites_B[[2]], Bodysites_B[[3]])

######################################################################
#Set Colors####
#cols <- brewer.pal(8,'Paired')
cols <- c("#cb1b4a", "#42aeb8", "#FDB316", "#c3c823", "#00797e", "#053058", "#aaada6", "#ae1848", "#368b90", "#ca9012", "#9ba11e", "#f1f2ec", "#d9d3df", "#348fbe", "#ff8340", "#ffAf40", "#bf503f", "#951b72", "#b75f6d")
cols2 <- colorRampPalette(cols)
cols_grad <- colorRampPalette(c("#cb1b4a", "#42aeb8"))
cols_grad2 <- colorRampPalette(c("#ca9012", "#cb1b4a"))

# If you want to color by family
fam_colors <- sample(cols2(length(unique(mapping$familyN))))
names(fam_colors) <- unique(mapping$familyN)

#set body site colors
#Set for yellow, teal, pink (skin, oral, anal)
body_cols <- c(cols2(5)[3], cols2(5)[2], cols2(5)[1])
body_cols2 <- c(cols2(5)[1], cols2(5)[2], cols2(5)[3])

######################################################################
#make full sample table of expected 
#There should be 17 families, (6 timepoints and 3 bodysites for infants, 2 bodysites and 1 timepoint for mom)

family_ids <- c()
fams <- sort(as.numeric(fams))
for(i in 1:length(fams)){
  family_working <- fams[i]
  mom_row <- rownames(mapping[mapping$familyN==family_working & mapping$motherorbaby_M_B == "M",])
  mom_id <- mapping[mom_row[1], "subject_id"]
  baby_row <- rownames(mapping[mapping$familyN==family_working & mapping$motherorbaby_M_B == "B",])
  baby_id <- mapping[baby_row[1], "subject_id"]
  family_ids[[i]] <- c(baby_id, mom_id)
  names(family_ids)[i] <- family_working
}

full_samples <- data.frame(rep(NA, 340), rep(NA, 340), rep(NA, 340), rep(NA, 340))
colnames(full_samples) <- c("Family", "Subject", "Time", "Bodysite")
full_samples$Family <- rep(as.numeric(fams), 20)
full_samples$Family <- sort(full_samples$Family)
full_samples$Time <- c(rep(c(rep(1, 3), rep(3, 3), rep(7, 3), rep(14, 3), rep(21, 3), rep(30, 3), rep(1, 2)), 17))
full_samples$Bodysite <- c(rep(c(rep(c("Skin", "Oral", "Anal"), 6), "Vaginal", "Anal"), 17))

subject_list <- c()
babes <- c()
moms <- c()
for(i in 1:length(family_ids)){
  subject_list <- c(subject_list, rep(family_ids[[i]][1], 18), rep(family_ids[[i]][2], 2))
  babes <- c(babes, family_ids[[i]][1])
  moms <- c(moms, family_ids[[i]][2])
}
full_samples$Subject <- subject_list
full_samples$Sample_ID <- c(rep(0, 340))
full_samples[is.na(full_samples)] <- 0 

subjects_try <- list(babes, moms)
sites_try <- list(c("Skin", "Oral", "Anal"),  c("Vaginal", "Anal")) 
for(x in 1:2){
  sample_types <- sites_try[[x]]
  subjects <- subjects_try[[x]]
  times_now <- list(c(1,3,7,14,21,30), c(1))[[x]]
  for(i in 1:length(unique(subjects))){
    if(! is.na(subjects[i])){
      now <- subjects[i]
      for(t in 1:length(times_now)){
        time <- times_now[t]
        for(s in 1:length(sample_types)){
          site <- sample_types[s]
          mapping_row <- rownames(mapping[mapping$subject_id ==now & mapping$planned_sampling_day_0_1_3_7_14_30_60_90___365 == time & mapping$SuperbodysiteOralSkinNoseVaginaAnalsAureola == site,])
          if(length(mapping_row) > 1){
            mapping_row <- mapping_row[1]
          }
          if(length(mapping_row) < 1){
            full_samples[full_samples$Subject==now & full_samples$Time == time & full_samples$Bodysite==site, "Sample_ID"] <- 0
          } else {
            full_samples[full_samples$Subject==now & full_samples$Time == time & full_samples$Bodysite==site, "Sample_ID"] <- mapping_row
          }
        }
      }
    }
  }
}
print(paste("number of samples missing: ", length(which(full_samples$Sample_ID == 0))))

