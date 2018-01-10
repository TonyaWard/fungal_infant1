#### Find and plot any differentiated taxa
source("bin/diff.test.r")

## Set new taxa table to be relative abundance
colSums(taxa_table)
taxa_table <- sweep(taxa_table, 2, colSums(taxa_table), FUN="/")

#set output dir
diff_dir <- paste(main_fp, "diff_taxa/Baby_VC/", sep='/')

#Set up tests to run
test.ixs <- list(Vaginal, Csection)
names(test.ixs) <- c("Vaginal", "Csection")

pvals<- c()
ALPHA <- 0.25
v_taxa <- taxa_table[,Vagina]
v_taxa <- v_taxa[rowSums(v_taxa) > 0,]

#For each bodysite and and timepoint run tests
for(i in 1:length(Bodysites_B)){
  Bodysite <- Bodysites_B[i]
  for(j in 1:length(Days)){
    Day <- Days[j]
    union1 <- intersect(Bodysite[[1]], Day[[1]])
    for(n in 1:(length(test.ixs)-1)){
      for(m in (n+1):length(test.ixs)){
        test.x <- test.ixs[n]
        test.y <- test.ixs[m]
        set1 <- intersect(union1, test.x[[1]])
        set2 <- intersect(union1, test.y[[1]])
        if(length(set1) > 2 && length(set2) > 2){
          full_set <- c(set1, set2)
          #keep vaginal taxa and the samples you are testing
          test_table <- t(taxa_table[rownames(v_taxa),full_set,drop=F])
          #Keep taxa that have at least one count
          test_table <- test_table[,colSums(test_table)>0, drop=F]
          map_test <- mapping[full_set,]
          map_test$Delivery_Vvaginal_Ccs_IcsInoc <- factor(map_test$Delivery_Vvaginal_Ccs_IcsInoc)
          
          #difftest <- differentiation.test(test_table, map_test$Delivery_Vvaginal_Ccs_IcsInoc, parametric=TRUE)
          difftest <- differentiation.test(test_table, map_test$Delivery_Vvaginal_Ccs_IcsInoc, parametric=FALSE)
          
          if(any(difftest$qvalues <= ALPHA)){
            signif.ix <- which(difftest$qvalues <= ALPHA)
            signif.ix <- signif.ix[order(difftest$pvalues[signif.ix])]
            for(k in 1:length(signif.ix)){
            #  if(!is.null(difftest$norm.test.pvals)){
            #    norm.test <- difftest$norm.test.pvals[k]
            #  } else {
            #    norm.test <- '0'
            #  }
            #  if(norm.test < 0.05){
            #    qval <- difftest$qvalues[k]
            #  } else {
            #    
              
              taxon <- names(signif.ix)[k]
              qval <- difftest$qvalues[[taxon]]
              pval <- difftest$pvalues[[taxon]]
              #name <- paste(names(Bodysite), names(Day), names(test.x), names(test.y), taxon, sep="-")
              name <- paste(names(Bodysite), names(Day), names(test.x), names(test.y), names(signif.ix)[k], sep="-")
              fp_name <- paste(diff_dir, name, sep="/")
              
              #Stats output  
              sink(paste(fp_name, "txt", sep="."))
              cat(paste('q=',qval,' taxon: ',taxon,'\n',sep=''))
              cat(paste('p=',pval,' taxon: ',taxon,'\n',sep=''))
              sink()
              
              #boxplots
              test_table <- as.data.frame(test_table)
              plot1 <- ggplot(test_table, aes(y=test_table[,names(signif.ix)[k]], x=map_test$Delivery_Vvaginal_Ccs_IcsInoc, fill=map_test$Delivery_Vvaginal_Ccs_IcsInoc)) +
                geom_boxplot(outlier.shape = NA) +
                geom_jitter(position=position_jitter(0.1), shape=1, size=3) +
                guides(fill=FALSE) +
                labs(x="", y = "Relative Abundance") +
                scale_fill_manual(values=c("#c3c823", "#053d58"))
              
              save_plot(paste(fp_name, "pdf", sep="."),plot1)
              #boxplot(test_table[,names(signif.ix)[k]] ~ map_test$Delivery_Vvaginal_Ccs_IcsInoc, 
              #        xlab='', ylab="Relative Abundance", main=name,
              #        col=c("#c3c823", "#053d58"))
              #dev.off()
            }
          } else {
            cat("not significant.")
            }
        } else {
          cat("Less than two samples in one group, skipping this test.")
        }
      }
    }
  }
}

######################################################################
#don't control for Day, test VC
#For each bodysite and and timepoint run tests
for(i in 1:length(Bodysites_B)){
  Bodysite <- Bodysites_B[i]
  union1 <- Bodysite[[1]]
    for(n in 1:(length(test.ixs)-1)){
      for(m in (n+1):length(test.ixs)){
        test.x <- test.ixs[n]
        test.y <- test.ixs[m]
        set1 <- intersect(union1, test.x[[1]])
        set2 <- intersect(union1, test.y[[1]])
        if(length(set1) > 2 && length(set2) > 2){
          full_set <- c(set1, set2)
          #keep vaginal taxa and the samples you are testing
          test_table <- t(taxa_table[rownames(v_taxa),full_set,drop=F])
          #Keep taxa that have at least one count
          test_table <- test_table[,colSums(test_table)>0, drop=F]
          map_test <- mapping[full_set,]
          map_test$Delivery_Vvaginal_Ccs_IcsInoc <- factor(map_test$Delivery_Vvaginal_Ccs_IcsInoc)
          
          #difftest <- differentiation.test(test_table, map_test$Delivery_Vvaginal_Ccs_IcsInoc, parametric=TRUE)
          difftest <- differentiation.test(test_table, map_test$Delivery_Vvaginal_Ccs_IcsInoc, parametric=FALSE)
          
          if(any(difftest$qvalues <= ALPHA)){
            signif.ix <- which(difftest$qvalues <= ALPHA)
            signif.ix <- signif.ix[order(difftest$pvalues[signif.ix])]
            for(k in 1:length(signif.ix)){
              #  if(!is.null(difftest$norm.test.pvals)){
              #    norm.test <- difftest$norm.test.pvals[k]
              #  } else {
              #    norm.test <- '0'
              #  }
              #  if(norm.test < 0.05){
              #    qval <- difftest$qvalues[k]
              #  } else {
              #    
              
              taxon <- names(signif.ix)[k]
              qval <- difftest$qvalues[[taxon]]
              pval <- difftest$pvalues[[taxon]]
              #name <- paste(names(Bodysite), names(Day), names(test.x), names(test.y), taxon, sep="-")
              name <- paste(names(Bodysite), names(test.x), names(test.y), names(signif.ix)[k], sep="-")
              fp_name <- paste(diff_dir, name, sep="/")
              
              #Stats output  
              sink(paste(fp_name, "txt", sep="."))
              cat(paste('q=',qval,' taxon: ',taxon,'\n',sep=''))
              cat(paste('p=',pval,' taxon: ',taxon,'\n',sep=''))
              sink()
              
              #boxplots
              test_table <- as.data.frame(test_table)
              plot1 <- ggplot(test_table, aes(y=test_table[,names(signif.ix)[k]], x=map_test$Delivery_Vvaginal_Ccs_IcsInoc, fill=map_test$Delivery_Vvaginal_Ccs_IcsInoc)) +
                geom_boxplot(outlier.shape = NA) +
                geom_jitter(position=position_jitter(0.1), shape=1, size=3) +
                guides(fill=FALSE) +
                labs(x="", y = "Relative Abundance") +
                scale_fill_manual(values=c("#c3c823", "#053d58"))
              
              save_plot(paste(fp_name, "pdf", sep="."),plot1)
              #boxplot(test_table[,names(signif.ix)[k]] ~ map_test$Delivery_Vvaginal_Ccs_IcsInoc, 
              #        xlab='', ylab="Relative Abundance", main=name,
              #        col=c("#c3c823", "#053d58"))
              #dev.off()
            }
          } else {
            cat("not significant.")
          }
        } else {
          cat("Less than two samples in one group, skipping this test.")
        }
      }
    }
}

##### Test bodysite in general ##########
###Test mom's sites####
diff_dir <- paste(main_fp, "diff_taxa/Mom_Sites/", sep='/')

pvals<- c()
ALPHA <- 0.25
mom <- c(Anal_M, Vagina)
test_table <- taxa_table[, mom]
test_table <- test_table[rowSums(test_table) > 0,]
test_table <- test_table[rowSums(test_table > 0 )/ncol(test_table) > 0.1, ]

test_table <- t(test_table)
map_test <- mapping[mom,]

difftest <- differentiation.test(test_table, map_test$Delivery_Vvaginal_Ccs_IcsInoc, parametric=TRUE)
difftest <- differentiation.test(test_table, map_test$SuperbodysiteOralSkinNoseVaginaAnalsAureola, parametric=FALSE)

if(any(difftest$qvalues <= ALPHA)){
  signif.ix <- which(difftest$qvalues <= ALPHA)
  signif.ix <- signif.ix[order(difftest$pvalues[signif.ix])]
  for(k in signif.ix){
      if(!is.null(difftest$norm.test.pvals)){
        norm.test <- difftest$norm.test.pvals[k]
      } else {
        norm.test <- '0'
      }
      if(norm.test < 0.05){
        qval <- difftest.np$qvalues[k]
        pval <- difftest.np$pvalues[k]
      } else {
        qval <- difftest$qvalues[k]
        pval <- difftest$pvalues[k]
      }
      taxon <- names(signif.ix)[k]
      #qval <- difftest$qvalues[[taxon]]
      pval <- difftest$pvalues[[taxon]]
      name <- paste("mom_V_A", names(signif.ix)[k], sep="-")
      fp_name <- paste(diff_dir, name, sep="/")
    
      #Stats output  
      sink(paste(fp_name, "txt", sep="."))
      cat(paste('q=',qval,'p=', pval , 'taxon: ',taxon,'\n',sep=''))
      sink()
    
      #boxplots
      test_table <- as.data.frame(test_table)
      plot1 <- ggplot(test_table, aes(y=test_table[,names(signif.ix)[k]], x=map_test$SuperbodysiteOralSkinNoseVaginaAnalsAureola, fill=map_test$SuperbodysiteOralSkinNoseVaginaAnalsAureola)) +
      geom_boxplot(outlier.shape = NA) +
      geom_jitter(position=position_jitter(0.1), shape=1, size=3) +
      guides(fill=FALSE) +
      labs(x="", y = "Relative Abundance") +
      scale_fill_manual(values=c("#b75f6d", "#99897e"))
    
    save_plot(paste(fp_name, "pdf", sep="."),plot1)
    }
} else {
  cat("No differentiated taxa for mom")
}


  

###Test Baby Sites, control for Day####
diff_dir <- paste(main_fp, "diff_taxa/Body_Sites/", sep='/')
test.ixs <- Bodysites_B

ALPHA <- 0.25
kids <- c(Bodysites_B[[1]], Bodysites_B[[2]], Bodysites_B[[3]])
test_table <- taxa_table[, kids]
test_table <- test_table[rowSums(test_table) > 0,]
test_table <- test_table[rowSums(test_table > 0 )/ncol(test_table) > 0.1, ]

for(i in 1:1){
  infants <- Infant
  for(j in 1:length(Days)){
    Day <- Days[j]
    union1 <- intersect(infants, Day[[1]])
    for(n in 1:(length(test.ixs)-1)){
      for(m in (n+1):length(test.ixs)){
        test.x <- test.ixs[n]
        test.y <- test.ixs[m]
        set1 <- intersect(union1, test.x[[1]])
        set2 <- intersect(union1, test.y[[1]])
        if(length(set1) > 2 && length(set2) > 2){
          full_set <- c(set1, set2)
          test_table2 <- t(test_table[,full_set])
          test_table2 <- test_table2[,colSums(test_table2)>0, drop=F]
          map_test <- mapping[full_set,]
          map_test$SuperbodysiteOralSkinNoseVaginaAnalsAureola <- factor(map_test$SuperbodysiteOralSkinNoseVaginaAnalsAureola)
          
          #difftest <- differentiation.test(test_table, map_test$Delivery_Vvaginal_Ccs_IcsInoc, parametric=TRUE)
          difftest <- differentiation.test(test_table2, map_test$SuperbodysiteOralSkinNoseVaginaAnalsAureola, parametric=FALSE)
          
          if(any(difftest$qvalues <= ALPHA)){
            signif.ix <- which(difftest$qvalues <= ALPHA)
            signif.ix <- signif.ix[order(difftest$pvalues[signif.ix])]
            for(k in 1:length(signif.ix)){
              #  if(!is.null(difftest$norm.test.pvals)){
              #    norm.test <- difftest$norm.test.pvals[k]
              #  } else {
              #    norm.test <- '0'
              #  }
              #  if(norm.test < 0.05){
              #    qval <- difftest$qvalues[k]
              #  } else {
              #    
              taxon <- names(signif.ix)[k]
              qval <- difftest$qvalues[[taxon]]
              pval <- difftest$pvalues[[taxon]]
              name <- paste(names(Day), names(test.x), names(test.y), names(signif.ix)[k], sep="-")
              fp_name <- paste(diff_dir, name, sep="/")
              file_name <- paste(fp_name, ".txt", sep="")
              
              #Stats output  
              sink(file = file_name)
              cat(paste('q=',qval,' taxon: ',taxon,'\n',sep=''))
              cat(paste('p=',pval,' taxon: ',taxon,'\n',sep=''))
              sink()
              
              #boxplots
              pdf(paste(fp_name, "pdf", sep="."),width=4,height=4)
              boxplot(test_table2[,names(signif.ix)[k]] ~ map_test$SuperbodysiteOralSkinNoseVaginaAnalsAureola, 
                      xlab='', ylab="Relative Abundance", main=name,
                      col=cols2(length(unique(map_test$SuperbodysiteOralSkinNoseVaginaAnalsAureola))))
              dev.off()
            }
          } else {
            cat("not significant.")
          }
        } else {
          cat("Less than two samples in one group, skipping this test.")
        }
      }
    }
  }
}
          
          
###Don't control for day #####

###Test Baby Sites####
diff_dir <- paste(main_fp, "diff_taxa/Body_Sites_gen/", sep='/')

test.ixs <- Bodysites_B
ALPHA <- 0.25
kids <- c(Bodysites_B[[1]], Bodysites_B[[2]], Bodysites_B[[3]])
test_table <- taxa_table[, kids]
test_table <- test_table[rowSums(test_table) > 0,]
test_table <- test_table[rowSums(test_table > 0 )/ncol(test_table) > 0.1, ]
map_test <- mapping[colnames(test_table),]
map_test$SuperbodysiteOralSkinNoseVaginaAnalsAureola <- factor(map_test$SuperbodysiteOralSkinNoseVaginaAnalsAureola, levels=c("Skin", "Oral", "Anal"))

difftest <- differentiation.test(t(test_table), map_test$SuperbodysiteOralSkinNoseVaginaAnalsAureola, parametric=FALSE)

if(any(difftest$qvalues <= ALPHA)){
  signif.ix <- which(difftest$qvalues <= ALPHA)
  signif.ix <- signif.ix[order(difftest$pvalues[signif.ix])]
  for(k in 1:length(signif.ix)){
    #  if(!is.null(difftest$norm.test.pvals)){
    #    norm.test <- difftest$norm.test.pvals[k]
    #  } else {
    #    norm.test <- '0'
    #  }
    #  if(norm.test < 0.05){
    #    qval <- difftest$qvalues[k]
    #  } else {
    #    
    taxon <- names(signif.ix)[k]
    qval <- difftest$qvalues[[taxon]]
    pval <- difftest$pvalues[[taxon]]
    name <- names(signif.ix)[k]
    fp_name <- paste(diff_dir, name, sep="/")
    file_name <- paste(fp_name, ".txt", sep="")
    
    #Stats output  
    sink(file = file_name)
    cat(paste('q=',qval,' taxon: ',taxon,'\n',sep=''))
    cat(paste('p=',pval,' taxon: ',taxon,'\n',sep=''))
    sink()
    working_table <- t(test_table)
    #boxplots
    working_table <- as.data.frame(working_table)
    plot1 <- ggplot(working_table, aes(y=working_table[,names(signif.ix)[k]], x=map_test$SuperbodysiteOralSkinNoseVaginaAnalsAureola, fill=map_test$SuperbodysiteOralSkinNoseVaginaAnalsAureola)) +
      geom_boxplot(outlier.shape = NA) +
      geom_jitter(position=position_jitter(0.1), shape=1, size=3) +
      guides(fill=FALSE) +
      labs(x="", y = "Relative Abundance") +
      scale_fill_manual(values=body_cols)
    save_plot(paste(fp_name, "pdf", sep="."), plot1)
    
    #pdf(paste(fp_name, "pdf", sep="."),width=4,height=4)
    #boxplot(working_table[,names(signif.ix)[k]] ~ map_test$SuperbodysiteOralSkinNoseVaginaAnalsAureola, 
    #        xlab='', ylab="Relative Abundance", main=name,
    #        col=cols2(length(unique(map_test$SuperbodysiteOralSkinNoseVaginaAnalsAureola))))
    #dev.off()
  }
}

diff_dir <- paste(main_fp, "diff_taxa/Body_Sites_gen/pairwise/", sep='/')
#Pairwise
for(i in 1:1){
  infants <- Infant
  union1 <- infants
    for(n in 1:(length(test.ixs)-1)){
      for(m in (n+1):length(test.ixs)){
        test.x <- test.ixs[n]
        test.y <- test.ixs[m]
        set1 <- intersect(union1, test.x[[1]])
        set2 <- intersect(union1, test.y[[1]])
        if(length(set1) > 2 && length(set2) > 2){
          full_set <- c(set1, set2)
          test_table2 <- t(test_table[,full_set])
          test_table2 <- test_table2[,colSums(test_table2)>0, drop=F]
          map_test <- mapping[full_set,]
          map_test$SuperbodysiteOralSkinNoseVaginaAnalsAureola <- factor(map_test$SuperbodysiteOralSkinNoseVaginaAnalsAureola)
          
          #difftest <- differentiation.test(test_table, map_test$Delivery_Vvaginal_Ccs_IcsInoc, parametric=TRUE)
          difftest <- differentiation.test(test_table2, map_test$SuperbodysiteOralSkinNoseVaginaAnalsAureola, parametric=FALSE)
          
          if(any(difftest$qvalues <= ALPHA)){
            signif.ix <- which(difftest$qvalues <= ALPHA)
            signif.ix <- signif.ix[order(difftest$pvalues[signif.ix])]
            for(k in 1:length(signif.ix)){
              #  if(!is.null(difftest$norm.test.pvals)){
              #    norm.test <- difftest$norm.test.pvals[k]
              #  } else {
              #    norm.test <- '0'
              #  }
              #  if(norm.test < 0.05){
              #    qval <- difftest$qvalues[k]
              #  } else {
              #    
              taxon <- names(signif.ix)[k]
              qval <- difftest$qvalues[[taxon]]
              pval <- difftest$pvalues[[taxon]]
              name <- paste(names(test.x), names(test.y), names(signif.ix)[k], sep="-")
              fp_name <- paste(diff_dir, name, sep="/")
              file_name <- paste(fp_name, ".txt", sep="")
              
              #Stats output  
              sink(file = file_name)
              cat(paste('q=',qval,' taxon: ',taxon,'\n',sep=''))
              cat(paste('p=',pval,' taxon: ',taxon,'\n',sep=''))
              sink()
              
              #boxplots
              pdf(paste(fp_name, "pdf", sep="."),width=4,height=4)
              boxplot(test_table2[,names(signif.ix)[k]] ~ map_test$SuperbodysiteOralSkinNoseVaginaAnalsAureola, 
                      xlab='', ylab="Relative Abundance", main=name,
                      col=cols2(length(unique(map_test$SuperbodysiteOralSkinNoseVaginaAnalsAureola))))
              dev.off()
            }
          } else {
            cat("not significant.")
          }
        } else {
          cat("Less than two samples in one group, skipping this test.")
        }
      }
    }
}

######################################################################
##Is Candida Parapsilosis higher in babies than moms?
diff_dir <- paste(main_fp, "diff_taxa/Baby_Mom_Parap/", sep='/')


for(b in 1:length(Bodysites_B)){
  test_table <- as.data.frame(t(taxa_table["Candida parapsilosis",]))
  baby_site <- Bodysites_B[[b]]
  baby_amount <- taxa_table["Candida parapsilosis",baby_site]
  for(m in 1:length(Bodysites_M)){
    mom_site <- Bodysites_M[[m]]
    mom_amount <- taxa_table["Candida parapsilosis",mom_site]
    test_out <- t.test(baby_amount, mom_amount)
    if(test_out$p.value < 0.05){
      keeps <- c(baby_site,mom_site)
      test_table <- test_table[keeps,,drop=F]
      test_map <- mapping[keeps,,drop=F]
      plot <- ggplot(test_table) +
        geom_boxplot(aes(x=test_map[,"motherorbaby_M_B"], y=test_table[,"Candida parapsilosis"], fill= factor(test_map[,"motherorbaby_M_B"])))+
        scale_fill_manual(values=c("#D9D3DF","#BF503F")) +
        guides(fill=FALSE) +
        labs(x="", y="Relative Abundance")
      fp <- paste(diff_dir, names(Bodysites_B)[b], names(Bodysites_M)[m], ".pdf", sep="")
      save_plot(fp, plot)
      fp_stats <- paste(diff_dir, names(Bodysites_B)[b], names(Bodysites_M)[m], ".txt", sep="")
      sink(fp_stats)
      print("p value:\n")
      cat(test_out$p.value)
      sink()
    }
  }
}
