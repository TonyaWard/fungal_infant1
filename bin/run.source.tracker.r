# This runs SourceTracker

#otu table should be counts, mapping file should have only samples in otu
#full mapping = metadata, full otu table = otutable
otus <- t(otutable) 

#add column for source/sink
metadata$SourceSink <- metadata$SampleID
metadata$SourceSink <- gsub("_s", "", metadata$SourceSink)
metadata$SourceSink <- gsub("[0-9]+", "sink", metadata$SourceSink)
metadata$SourceSink <- gsub("Positive", "sink", metadata$SourceSink)
metadata$SourceSink <- gsub("Negative", "source", metadata$SourceSink)
metadata$Env <- metadata$SuperbodysiteOralSkinNoseVaginaAnalsAureola
metadata$Env[323] <- "Negative"
metadata$Env[324] <- "Postive"

source_info <-metadata[metadata$SourceSink == "source",]
metadata <- rbind(metadata, source_info)
metadata <- rbind(metadata, source_info)
rownames(metadata)[325] <- "Negative2_s"
rownames(metadata)[326] <- "Negative3_s"
source_info <- otus[rownames(otus) == "Negative_s",]
otus <- rbind(otus, source_info)
otus <- rbind(otus, source_info)
rownames(otus)[325] <- "Negative2_s"
rownames(otus)[326] <- "Negative3_s"

metadata <- metadata[rownames(otus),]

# extract the source environments and source/sink indices
train.ix <- which(metadata$SourceSink=='source')
test.ix <- which(metadata$SourceSink=='sink')

envs <- metadata$Env
if(is.element('Description',colnames(metadata))) desc <- metadata$Description

# load SourceTracker package
source('src/SourceTracker.r')

# tune the alpha values using cross-validation (this is slow!)
# tune.results <- tune.st(otus[train.ix,], envs[train.ix])
# alpha1 <- tune.results$best.alpha1
# alpha2 <- tune.results$best.alpha2
# note: to skip tuning, run this instead:
alpha1 <- alpha2 <- 0.001

# train SourceTracker object on training data
st <- sourcetracker(otus[train.ix,], envs[train.ix], rarefaction_depth=10)

# Estimate source proportions in test data
results <- predict(st,otus[test.ix,], alpha1=alpha1, alpha2=alpha2, rarefaction_depth=10, full=TRUE)

# get average of full results across restarts
res.mean <- apply(results$full.results,c(2,3,4),mean)

# Get depth of each sample for relative abundance calculation
sample.depths <- apply(results$full.results[1,,,,drop=F],4,sum)

# write each env separate file
for(i in 1:length(results$train.envs)){
  env.name <- results$train.envs[i]
  filename.fractions <- sprintf('%s_contributions.txt', env.name)
  res.mean.i <- res.mean[i,,]
  # handle the case where there is only one sink sample
  if(is.null(dim(res.mean.i))) res.mean.i <- matrix(res.mean.i,ncol=1)
  
  # make rows be samples, columns be features
  res.mean.i <- t(res.mean.i)
  
  # ensure proper names
  colnames(res.mean.i) <- colnames(otus)
  rownames(res.mean.i) <- results$samplenames
  
  # calculate and save relative abundance
  res.mean.i.ra <- sweep(res.mean.i,1,sample.depths,'/')
  sink(filename.fractions)
  cat('SampleID\t')
  write.table(res.mean.i.ra,quotes=F,sep='\t')
  sink(NULL)
}
