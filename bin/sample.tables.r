# Make Supplemental Tables
######################################################################
#Predominant taxa by sample type
sites <- list(Skin, Oral, Anal_B, Anal_M, Vagina)
names(sites) <- c("Skin", "Oral", "Anal_B", "Anal_M", "Vagina")
for(s in 1:length(sites)){
  site <- sites[[s]]
  table_now <- taxa_table[,site]
  table_now$mean <- rowMeans(table_now)
  table_now2 <- table_now[with(table_now, order(-mean)),]
  file_name <- paste(names(sites)[s], ".txt", sep="")
  write.table(table_now2[1:5, "mean", drop=F], file_name, sep="\t", quote=F, row.names = T)
}

######################################################################
#Sample info - each body site: min, mean, max sequence info, samples passing filter, samples dropped, #OTUs, #Taxa

meta_orginal <- sample_metadata(read_biom(otu_fp))
meta_orginal$MapDepth <- as.numeric(meta_orginal$MapDepth)
meta_orginal <- meta_orginal[!meta_orginal$Delivery_Vvaginal_Ccs_IcsInoc == "I" & ! meta_orginal$bodysite_OralMucosa_Forehead_VolarRight_PalmRight_FootRight_Vaginal_Anal_Feces =="Feces",]


# original table doesn't make sense. Let's do:
# Mean Depth for sample types (Raw) +/- SE
# Mean sequences aligned for sample types +/- SE <- from metadata  (take out innocs)
# Total Number of Samples
# Number of samples dropped (<50 counts aligned) <- from metadata
# Mean OTUs +/- SE
# Mean Taxa +/- SE

## Get number of raw read and number of reads passing QC
raws <- read.table("data/raw_linecounts_split.tsv", sep="\t", header=T) #divide these by 4
raws$ReadCnt <- raws$ReadCnt/4
raws$SampID <- as.character(raws$SampID)
rownames(raws) <- raws$SampID
qcs <- read.table("data/qc_linecounts_fasta_nin.tsv", sep="\t", header=T) #divide these by 2
qcs$ReadCnt <- qcs$ReadCnt/2
qcs$SampID <- as.character(qcs$SampID)
rownames(qcs) <- qcs$SampID

#Subset to keep only samples in the meta data original --> no innoc/fec
raws <- raws[rownames(meta_orginal),]
qcs <- qcs[rownames(meta_orginal),]
meta_orginal$raws <- raws$ReadCnt
meta_orginal$qcs <- qcs$ReadCnt

#Raw reads
raw_skin <- mean(meta_orginal[meta_orginal$SuperbodysiteOralSkinNoseVaginaAnalsAureola == "Skin","raws"])
raw_oral <- mean(meta_orginal[meta_orginal$SuperbodysiteOralSkinNoseVaginaAnalsAureola == "Oral","raws"])
raw_analB <- mean(meta_orginal[meta_orginal$SuperbodysiteOralSkinNoseVaginaAnalsAureola == "Anal" & meta_orginal$motherorbaby_M_B == "B","raws"])
raw_analM <- mean(meta_orginal[meta_orginal$SuperbodysiteOralSkinNoseVaginaAnalsAureola == "Anal" & meta_orginal$motherorbaby_M_B == "M","raws"])
raw_vaginal <- mean(meta_orginal[meta_orginal$SuperbodysiteOralSkinNoseVaginaAnalsAureola == "Vaginal" & meta_orginal$motherorbaby_M_B == "M","raws"])

sd_rawskin <- sd(meta_orginal[meta_orginal$SuperbodysiteOralSkinNoseVaginaAnalsAureola == "Skin","raws"])
sd_raworal <- sd(meta_orginal[meta_orginal$SuperbodysiteOralSkinNoseVaginaAnalsAureola == "Oral","raws"])
sd_rawanalB <- sd(meta_orginal[meta_orginal$SuperbodysiteOralSkinNoseVaginaAnalsAureola == "Anal" & meta_orginal$motherorbaby_M_B == "B","raws"])
sd_rawanalM <- sd(meta_orginal[meta_orginal$SuperbodysiteOralSkinNoseVaginaAnalsAureola == "Anal" & meta_orginal$motherorbaby_M_B == "M","raws"])
sd_rawvaginal <- sd(meta_orginal[meta_orginal$SuperbodysiteOralSkinNoseVaginaAnalsAureola == "Vaginal" & meta_orginal$motherorbaby_M_B == "M","raws"])

#Reads QC
qc_skin <- mean(meta_orginal[meta_orginal$SuperbodysiteOralSkinNoseVaginaAnalsAureola == "Skin","qcs"])
qc_oral <- mean(meta_orginal[meta_orginal$SuperbodysiteOralSkinNoseVaginaAnalsAureola == "Oral","qcs"])
qc_analB <- mean(meta_orginal[meta_orginal$SuperbodysiteOralSkinNoseVaginaAnalsAureola == "Anal" & meta_orginal$motherorbaby_M_B == "B","qcs"])
qc_analM <- mean(meta_orginal[meta_orginal$SuperbodysiteOralSkinNoseVaginaAnalsAureola == "Anal" & meta_orginal$motherorbaby_M_B == "M","qcs"])
qc_vaginal <- mean(meta_orginal[meta_orginal$SuperbodysiteOralSkinNoseVaginaAnalsAureola == "Vaginal" & meta_orginal$motherorbaby_M_B == "M","qcs"])

sd_qcskin <- sd(meta_orginal[meta_orginal$SuperbodysiteOralSkinNoseVaginaAnalsAureola == "Skin","qcs"])
sd_qcoral <- sd(meta_orginal[meta_orginal$SuperbodysiteOralSkinNoseVaginaAnalsAureola == "Oral","qcs"])
sd_qcanalB <- sd(meta_orginal[meta_orginal$SuperbodysiteOralSkinNoseVaginaAnalsAureola == "Anal" & meta_orginal$motherorbaby_M_B == "B","qcs"])
sd_qcanalM <- sd(meta_orginal[meta_orginal$SuperbodysiteOralSkinNoseVaginaAnalsAureola == "Anal" & meta_orginal$motherorbaby_M_B == "M","qcs"])
sd_qcvaginal <- sd(meta_orginal[meta_orginal$SuperbodysiteOralSkinNoseVaginaAnalsAureola == "Vaginal" & meta_orginal$motherorbaby_M_B == "M","qcs"])


#Reads Aligning
mean_skin <- mean(meta_orginal[meta_orginal$SuperbodysiteOralSkinNoseVaginaAnalsAureola == "Skin","MapDepth"])
mean_oral <- mean(meta_orginal[meta_orginal$SuperbodysiteOralSkinNoseVaginaAnalsAureola == "Oral","MapDepth"])
mean_analB <- mean(meta_orginal[meta_orginal$SuperbodysiteOralSkinNoseVaginaAnalsAureola == "Anal" & meta_orginal$motherorbaby_M_B == "B","MapDepth"])
mean_analM <- mean(meta_orginal[meta_orginal$SuperbodysiteOralSkinNoseVaginaAnalsAureola == "Anal" & meta_orginal$motherorbaby_M_B == "M","MapDepth"])
mean_vaginal <- mean(meta_orginal[meta_orginal$SuperbodysiteOralSkinNoseVaginaAnalsAureola == "Vaginal" & meta_orginal$motherorbaby_M_B == "M","MapDepth"])

sd_skin <- sd(meta_orginal[meta_orginal$SuperbodysiteOralSkinNoseVaginaAnalsAureola == "Skin","MapDepth"])
sd_oral <- sd(meta_orginal[meta_orginal$SuperbodysiteOralSkinNoseVaginaAnalsAureola == "Oral","MapDepth"])
sd_analB <- sd(meta_orginal[meta_orginal$SuperbodysiteOralSkinNoseVaginaAnalsAureola == "Anal" & meta_orginal$motherorbaby_M_B == "B","MapDepth"])
sd_analM <- sd(meta_orginal[meta_orginal$SuperbodysiteOralSkinNoseVaginaAnalsAureola == "Anal" & meta_orginal$motherorbaby_M_B == "M","MapDepth"])
sd_vaginal <- sd(meta_orginal[meta_orginal$SuperbodysiteOralSkinNoseVaginaAnalsAureola == "Vaginal" & meta_orginal$motherorbaby_M_B == "M","MapDepth"])

#Number samples
total_s <- nrow(meta_orginal[meta_orginal$SuperbodysiteOralSkinNoseVaginaAnalsAureola == "Skin",])
total_o <- nrow(meta_orginal[meta_orginal$SuperbodysiteOralSkinNoseVaginaAnalsAureola == "Oral",])
total_ab <- nrow(meta_orginal[meta_orginal$SuperbodysiteOralSkinNoseVaginaAnalsAureola == "Anal" ,])
total_am <- nrow(meta_orginal[meta_orginal$SuperbodysiteOralSkinNoseVaginaAnalsAureola == "Anal",])
total_v <- nrow(meta_orginal[meta_orginal$SuperbodysiteOralSkinNoseVaginaAnalsAureola == "Vaginal",])

#Number dropped
below_50s <- nrow(meta_orginal[meta_orginal$SuperbodysiteOralSkinNoseVaginaAnalsAureola == "Skin" & meta_orginal$MapDepth < 50,])
below_50o <- nrow(meta_orginal[meta_orginal$SuperbodysiteOralSkinNoseVaginaAnalsAureola == "Oral" & meta_orginal$MapDepth < 50,])
below_50ab <- nrow(meta_orginal[meta_orginal$SuperbodysiteOralSkinNoseVaginaAnalsAureola == "Anal" & meta_orginal$MapDepth < 50 & meta_orginal$motherorbaby_M_B == "B",])
below_50am <- nrow(meta_orginal[meta_orginal$SuperbodysiteOralSkinNoseVaginaAnalsAureola == "Anal" & meta_orginal$MapDepth < 50 & meta_orginal$motherorbaby_M_B == "M",])
below_50v <- nrow(meta_orginal[meta_orginal$SuperbodysiteOralSkinNoseVaginaAnalsAureola == "Vaginal" & meta_orginal$MapDepth < 50,])

#Number OTUs
mapping$nOTUs <- as.numeric(as.character(mapping$nOTUs))
nOTUs_skin <- mean(mapping[mapping$SuperbodysiteOralSkinNoseVaginaAnalsAureola == "Skin", "nOTUs"])
nOTUs_oral <- mean(mapping[mapping$SuperbodysiteOralSkinNoseVaginaAnalsAureola == "Oral", "nOTUs"])
nOTUs_analB <- mean(mapping[mapping$SuperbodysiteOralSkinNoseVaginaAnalsAureola == "Anal" & mapping$motherorbaby_M_B=="B", "nOTUs"])
nOTUs_analM <- mean(mapping[mapping$SuperbodysiteOralSkinNoseVaginaAnalsAureola == "Anal" & mapping$motherorbaby_M_B== "M", "nOTUs"])
nOTUs_Vagina <- mean(mapping[mapping$SuperbodysiteOralSkinNoseVaginaAnalsAureola == "Vaginal", "nOTUs"])

sd_otus <- sd(mapping[mapping$SuperbodysiteOralSkinNoseVaginaAnalsAureola == "Skin", "nOTUs"])
sd_otuo <- sd(mapping[mapping$SuperbodysiteOralSkinNoseVaginaAnalsAureola == "Oral", "nOTUs"])
sd_otuab <- sd(mapping[mapping$SuperbodysiteOralSkinNoseVaginaAnalsAureola == "Anal" & mapping$motherorbaby_M_B=="B", "nOTUs"])
sd_otuam <- sd(mapping[mapping$SuperbodysiteOralSkinNoseVaginaAnalsAureola == "Anal" & mapping$motherorbaby_M_B== "M", "nOTUs"])
sd_otuv <- sd(mapping[mapping$SuperbodysiteOralSkinNoseVaginaAnalsAureola == "Vaginal", "nOTUs"])

#Number Taxa
mapping$nTaxa <- as.numeric(as.character(mapping$nTaxa))
taxa_skin <- mean(mapping[mapping$SuperbodysiteOralSkinNoseVaginaAnalsAureola == "Skin", "nTaxa"])
taxa_oral <- mean(mapping[mapping$SuperbodysiteOralSkinNoseVaginaAnalsAureola == "Oral", "nTaxa"])
taxa_analb <- mean(mapping[mapping$SuperbodysiteOralSkinNoseVaginaAnalsAureola == "Anal" & mapping$motherorbaby_M_B=="B", "nTaxa"])
taxa_analm <- mean(mapping[mapping$SuperbodysiteOralSkinNoseVaginaAnalsAureola == "Anal" & mapping$motherorbaby_M_B== "M", "nTaxa"])
taxa_vagina <- mean(mapping[mapping$SuperbodysiteOralSkinNoseVaginaAnalsAureola == "Vaginal", "nTaxa"])

sd_taxas <- sd(mapping[mapping$SuperbodysiteOralSkinNoseVaginaAnalsAureola == "Skin", "nTaxa"])
sd_taxao <- sd(mapping[mapping$SuperbodysiteOralSkinNoseVaginaAnalsAureola == "Oral", "nTaxa"])
sd_taxaab <- sd(mapping[mapping$SuperbodysiteOralSkinNoseVaginaAnalsAureola == "Anal" & mapping$motherorbaby_M_B=="B", "nTaxa"])
sd_taxaam <- sd(mapping[mapping$SuperbodysiteOralSkinNoseVaginaAnalsAureola == "Anal" & mapping$motherorbaby_M_B== "M", "nTaxa"])
sd_taxav <- sd(mapping[mapping$SuperbodysiteOralSkinNoseVaginaAnalsAureola == "Vaginal", "nTaxa"])





