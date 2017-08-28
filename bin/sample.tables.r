# Make Supplemental Tables
######################################################################
#Predominant (highest mean RA) taxa by sample type
sites <- list(Skin, Oral, Anal_B, Anal_M, Vagina)
names(sites) <- c("Skin", "Oral", "Anal_B", "Anal_M", "Vagina")
for(s in 1:length(sites)){
  site <- sites[[s]]
  table_now <- taxa_table[,site]
  table_now$mean <- rowMeans(table_now)
  table_now2 <- table_now[with(table_now, order(-mean)),]
  file_name <- paste(main_fp, "/", "MeanRA_", names(sites)[s], ".txt", sep="")
  write.table(table_now2[1:25, "mean", drop=F], file_name, sep="\t", quote=F, row.names = T)
}

#Most prevalent taxa by sample type
sites <- list(Skin, Oral, Anal_B, Anal_M, Vagina)
names(sites) <- c("Skin", "Oral", "Anal_B", "Anal_M", "Vagina")
for(s in 1:length(sites)){
  site <- sites[[s]]
  table_now <- taxa_table[,site]
  table_now$count <- rowSums(table_now > 0)
  table_now2 <- table_now[with(table_now, order(-count)),]
  file_name <- paste(main_fp, "/", "Prevalent_", names(sites)[s], ".txt", sep="")
  write.table(table_now2[1:25, "count", drop=F], file_name, sep="\t", quote=F, row.names = T)
}

#what percent is Malassezia?
site <- sites[[1]]
table_now <- taxa_table[,site]
table_now$mean <- rowMeans(table_now)
table_now2 <- table_now[with(table_now, order(-mean)),]
collapse_names <- str_split_fixed(rownames(table_now2), " ", 2)[,1]
genera <- as.matrix(table_now2)
rownames(genera) <- collapse_names
genera2 <- aggregate(genera, by=list(rownames(genera)), FUN = sum)
print("mean relative abundance of skin malassezia:")
print(genera2[genera2$Group.1 =="Malassezia","mean"])
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

#Subset to keep only samples in the meta data original --> no innoc/fec in here (line 43)
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
meta_orginal$SampleID <- rownames(meta_orginal) <- paste(rownames(meta_orginal), "s", sep="_")
total_s <- ncol(otutable[,meta_orginal[meta_orginal$SuperbodysiteOralSkinNoseVaginaAnalsAureola == "Skin" & meta_orginal$motherorbaby_M_B == "B", "SampleID"]])
total_o <- nncol(otutable[,meta_orginal[meta_orginal$SuperbodysiteOralSkinNoseVaginaAnalsAureola == "Oral" & meta_orginal$motherorbaby_M_B == "B", "SampleID"]])
total_ab <- ncol(otutable[,meta_orginal[meta_orginal$SuperbodysiteOralSkinNoseVaginaAnalsAureola == "Anal" & meta_orginal$motherorbaby_M_B == "B", "SampleID"]])
total_am <- ncol(otutable[,meta_orginal[meta_orginal$SuperbodysiteOralSkinNoseVaginaAnalsAureola == "Anal" & meta_orginal$motherorbaby_M_B == "M", "SampleID"]])
total_v <- ncol(otutable[,meta_orginal[meta_orginal$SuperbodysiteOralSkinNoseVaginaAnalsAureola == "Vaginal" & meta_orginal$motherorbaby_M_B == "M", "SampleID"]])

#Number dropped
below_50s <- length(which(colSums(otutable[,meta_orginal[meta_orginal$SuperbodysiteOralSkinNoseVaginaAnalsAureola == "Vaginal", "SampleID"]]) < 50 ))
below_50o <- length(which(colSums(otutable[,meta_orginal[meta_orginal$SuperbodysiteOralSkinNoseVaginaAnalsAureola == "Oral", "SampleID"]]) < 50 ))
below_50ab <- length(which(colSums(otutable[,meta_orginal[meta_orginal$SuperbodysiteOralSkinNoseVaginaAnalsAureola == "Anal" & meta_orginal$motherorbaby_M_B == "M", "SampleID"]]) < 50 ))
below_50am <- length(which(colSums(otutable[,meta_orginal[meta_orginal$SuperbodysiteOralSkinNoseVaginaAnalsAureola == "Anal" & meta_orginal$motherorbaby_M_B == "M", "SampleID"]]) < 50 ))
below_50v <- length(which(colSums(otutable[,meta_orginal[meta_orginal$SuperbodysiteOralSkinNoseVaginaAnalsAureola == "Vaginal", "SampleID"]]) < 50 ))

#Number OTUs
mapping$nOTUs <- as.numeric(as.character(mapping$nOTUs))
mean(colSums(otutable[,meta_orginal[meta_orginal$SuperbodysiteOralSkinNoseVaginaAnalsAureola == "Vaginal", "SampleID"]] > 0))
nOTUs_skin <- mean(colSums(otutable[,meta_orginal[meta_orginal$SuperbodysiteOralSkinNoseVaginaAnalsAureola == "Skin" & meta_orginal$motherorbaby_M_B == "B", "SampleID"]] > 0))
nOTUs_oral <- mean(colSums(otutable[,meta_orginal[meta_orginal$SuperbodysiteOralSkinNoseVaginaAnalsAureola == "Oral" & meta_orginal$motherorbaby_M_B == "B", "SampleID"]] > 0))
nOTUs_analB <- mean(colSums(otutable[,meta_orginal[meta_orginal$SuperbodysiteOralSkinNoseVaginaAnalsAureola == "Anal" & meta_orginal$motherorbaby_M_B == "B", "SampleID"]] > 0))
nOTUs_analM <- mean(colSums(otutable[,meta_orginal[meta_orginal$SuperbodysiteOralSkinNoseVaginaAnalsAureola == "Anal" & meta_orginal$motherorbaby_M_B == "M", "SampleID"]] > 0))
nOTUs_Vagina <- mean(colSums(otutable[,meta_orginal[meta_orginal$SuperbodysiteOralSkinNoseVaginaAnalsAureola == "Vaginal", "SampleID"]] > 0))

sd_otus <- sd(colSums(otutable[,meta_orginal[meta_orginal$SuperbodysiteOralSkinNoseVaginaAnalsAureola == "Skin" & meta_orginal$motherorbaby_M_B == "B", "SampleID"]] > 0))
sd_otuo <- sd(colSums(otutable[,meta_orginal[meta_orginal$SuperbodysiteOralSkinNoseVaginaAnalsAureola == "Oral" & meta_orginal$motherorbaby_M_B == "B", "SampleID"]] > 0))
sd_otuab <- sd(colSums(otutable[,meta_orginal[meta_orginal$SuperbodysiteOralSkinNoseVaginaAnalsAureola == "Anal" & meta_orginal$motherorbaby_M_B == "B", "SampleID"]] > 0))
sd_otuam <- sd(colSums(otutable[,meta_orginal[meta_orginal$SuperbodysiteOralSkinNoseVaginaAnalsAureola == "Anal" & meta_orginal$motherorbaby_M_B == "M", "SampleID"]] > 0))
sd_otuv <- sd(colSums(otutable[,meta_orginal[meta_orginal$SuperbodysiteOralSkinNoseVaginaAnalsAureola == "Vaginal", "SampleID"]] > 0))

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





