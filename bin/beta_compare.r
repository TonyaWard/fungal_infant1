#Compare within and between distances for beta diversity
"beta_compare" <- function(beta_table, set_1, set_2, color_vals){
  for(i in 1:length(set_1)){
    test_this <- set_1[i]
    for(j in 1:length(Days)){
      test_this2 <- set_2[j]
      union1 <- intersect(test_this[[1]], test_this2[[1]])
      for(n in 1:(length(test.ixs)-1)){
        for(m in (n+1):length(test.ixs)){
          test.x <- test.ixs[n]
          test.y <- test.ixs[m]
          set1 <- intersect(union1, test.x[[1]])
          set2 <- intersect(union1, test.y[[1]])
          full_set <- c(set1, set2)
          ktest <- c()
          set1_within <- c()
          set2_within <- c()
          between_sets <- c()
          set1_mom <- c()
          set2_mom <- c()
          
          if(length(set1) > 2 && length(set2) > 2){
            for(c in 1:length(colnames(beta_table))){
              for(r in (c+1):length(rownames(beta_table))){
                if(rownames(beta_table)[r] %in% set1 && colnames(beta_table)[c] %in% set1){
                  set1_within <- c(set1_within, beta_table[r,c])
                } else {
                  if(rownames(beta_table)[r] %in% set1 && colnames(beta_table)[c] %in% set2){
                    between_sets <- c(between_sets, beta_table[r,c])
                  } else {
                    if(rownames(beta_table)[r] %in% set1 && colnames(beta_table)[c] %in% mom){
                      set1_mom <- c(set1_mom, beta_table[r,c])
                    } else {
                      if(rownames(beta_table)[r] %in% set2 && colnames(beta_table)[c] %in% set1){
                        between_sets <- c(between_sets, beta_table[r,c])
                      } else {
                        if(rownames(beta_table)[r] %in% set2 && colnames(beta_table)[c] %in% set2){
                          set2_within <- c(set2_within, beta_table[r,c])
                        } else {
                          if(rownames(beta_table)[r] %in% set2 && colnames(beta_table)[c] %in% mom){
                            set2_mom <- c(set2_mom, beta_table[r,c])
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
            sets_test <- c()
            sets_test <- list(set2_within, between_sets, set1_within)
            names(sets_test) <- c(names(test.y), paste(names(test.y), names(test.x), sep='-'), names(test.x)) 
            ktest <- kruskal.test(sets_test)
            pvals <- c(pvals,ktest$p.value)
            
            mom_tests <- c()
            mom_tests <- list(set1_mom, set2_mom)
            names(mom_tests) <- c(names(test.y), names(test.x))
            wtest <- wilcox.test(mom_tests[[1]], mom_tests[[2]])
            pvals_mom <- c(pvals_mom,wtest$p.value)
            
            #print stats to screen
            cat(sprintf('\n%s, %s:\n',names(test_this), names(test_this2)))
            print(ktest)
            print(wtest)
            
            #write stats to file
            sink(file_name, append =TRUE)
            cat(sprintf('\n%s, %s:\n',names(test_this), names(test_this2)))
            print(ktest)
            cat("V and C to mom")
            print(wtest)
            sink()
            
            #assign pdf name for plot
            name1 <- paste(names(test_this), names(test_this2), sep="-")
            name1 <- paste(name1, ".pdf", sep='')
            file_path <- paste(beta_dir, name1, sep='')
            pdf(file_path, height=4,width=6)
            
            #make beta div box plots
            title <- sprintf('%s, %s',names(test_this), names(test_this2))
            boxplot(sets_test,
                    xlab='',ylab='Bray Curtis', main=title,
                    col=cols2(3))
            dev.off()
            
            #assign pdf name for plot
            name2 <- paste(names(test_this), names(test_this2), "to_moms", sep="-")
            name2 <- paste(name2, ".pdf", sep='')
            file_path <- paste(beta_dir, name2, sep='')
            pdf(file_path, height=4,width=6)
            
            #make beta div box plots
            title <- sprintf('%s, %s, %s', names(test_this), names(test_this2), "to moms")
            boxplot(mom_tests,
                    xlab='',ylab='Bray Curtis', main=title,
                    col=cols2(2))
            dev.off()
            
            
            names(sets_test) <- c(paste(names(test.y), names(test_this2), names(test_this), sep="_"), paste(paste(names(test.x), names(test.y), sep='-'), names(test_this2), names(test_this), sep="_"), paste(names(test.x), names(test_this2), names(test_this), sep="_"))
            distances <- rbind(distances,melt(sets_test))
            
            names(mom_tests) <- c(paste(names(mom_tests[1]), names(test_this2), names(test_this), sep="_"), paste(names(mom_tests[2]), names(test_this2), names(test_this), sep="_"))
            distances2 <- rbind(distances2,melt(mom_tests))
            
          } else {
            ktest <- "Less than two samples in one group, skipping this test."
          }
        }
      }
    }
  }
  #print fdr corrected pvals to stats file
  print(pvals)
  fdr.pvals <- p.adjust(pvals, method="fdr")
  fdr.pvals.mom <- p.adjust(pvals_mom, method="fdr")
  sink(file_name, append =TRUE)
  cat("\nfdr adjusted pvals:")
  print(fdr.pvals)
  cat("\nfdr adjusted pvals to mom:")
  print(fdr.pvals.mom)
  sink()
  print(fdr.pvals)
  return(list(distances, distances2))
}