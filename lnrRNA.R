## By following Lao Wang instruction.
## Once having the selected genes, multiply the betas to the corresponding gene value in the dataset to get
## a total risk value. Then, it gets split into high and low or 1 and 0 for all the rows at the same time.
## Using this new data to run survdiff function. Find the p values for all the combinations of one-less selected
## genes in each cycle and take the model with smallest p value go on to the next cycle. Again with one-less gene to 
## find the model with smallest p value... till the last a couple of genes stay in the model.
## 2018-8-16

options(warn=-1)

setwd(choose.dir())

selectlnr <- function(){
  if (!require(tools)){install.packages('tools')}
  if (!require(survival)) {install.packages('survival')}
  library(tools)
  library(survival)

  f = file.choose()
  if (file_ext(f) == 'csv') {dat <- read.csv(f, header = T, check.names=F)}
  if (file_ext(f) == 'tsv') {dat <- read.table(f, header = T, sep = '\t', check.names=F)}
  else {
    'Please save your file in csv or tsv format before proceeding!'
  }

  # Make sure the colnames complying with R convention.
  colnames(dat) <- gsub("[-: ]", "_", colnames(dat))
  dat <- as.data.frame(dat[, -1])

  # Save betas in a dataframe.
  betas <- setNames(data.frame(matrix(NA,ncol=3)), c('gene', 'beta', 'P_value'))

  for (i in c(5:length(dat))){
    surv1 <- coxph(Surv(OS_Mon, OS_Statu) ~ as.matrix(dat[,i]), data = dat, method = 'efron')
    smy1 <- coef(summary(surv1))
    betas[(i-4), ] = c(colnames(dat)[i], smy1[1], smy1[5])
  }
  write.csv(betas, paste0(format(Sys.time(), "%Y%m%d_%H%M%S_"), 'All_single_betas.csv'))

  # Let's select the variables that gives p value <= 0.05.
  selected_ind <- which(betas$P_value <= 0.05)
  write.csv(betas$gene[selected_ind], paste0(format(Sys.time(), "%Y%m%d_%H%M%S_"),'Selected_genes.csv'))
  dat2 <- dat[, c('OS_Statu', 'OS_Mon', betas$gene[selected_ind])]
  selected_betas <- betas[selected_ind, c(1,2)]  
  dat3 <- dat2  

  all_dropped <- c()
  min_pvalues <- c()

  for (j in (1:(length(selected_betas[,1])-2))){
    dropped_genes <- setNames(data.frame(matrix(NA, ncol = 4)), c('cycle', 'loop', 'dropped_gene', 'P-value'))
    min_p_gene <- setNames(data.frame(matrix(NA, ncol = 3)), c('cycle', 'drop_gene', 'p_value'))
    pvalues <- c()

    for (k in (1:length(selected_betas[,1]))){
      bet <- selected_betas[, 2][-k]
      dat4 <- dat3[, -(k + 2)]
      mat4 <- as.matrix(dat4[, -c(1,2)]) %*% as.matrix(as.numeric(bet))
      m <- ifelse(mat4 >= median(mat4), 1, 0)

      sf1 = survdiff(Surv(OS_Mon, OS_Statu) ~ m, data = dat4, rho = 0)
      pvalue = pchisq(sf1$chisq, length(sf1$n)-1, lower.tail = FALSE)
      dropped_genes[k, 1] = paste0('cycle',j)
      dropped_genes[k, 2] = paste0('loop', k)
      dropped_genes[k, 3] = selected_betas[, 1][k]
      dropped_genes[k, 4] = pvalue
      pvalues <- c(pvalues, pvalue)
     }

    min_pvalue_ind <- which(pvalues == min(pvalues))[1]
    drop_gene <- selected_betas[,1][min_pvalue_ind]
    selected_betas <- selected_betas[- min_pvalue_ind, ]  
    dat3 <- dat3[, -(min_pvalue_ind + 2)]  
    min_p_gene[j, ] = c(j, drop_gene, min(pvalues))

    all_dropped <- rbind(all_dropped, dropped_genes)
    min_pvalues <- rbind(min_pvalues, na.omit(min_p_gene))
   }
   write.csv(all_dropped, paste0(format(Sys.time(), "%Y%m%d_%H%M%S_"),'All_Dropped_Genes.csv'))
   write.csv(min_pvalues, paste0(format(Sys.time(), "%Y%m%d_%H%M%S_"),'Min_p_gene.csv'))
   print("The process is done. You may find the results in four files in the working directory.")
   cat(paste('All_single_betas: coefs from coxph with all individual genes of the original dataset.',
            'Selected_genes: selected genes based on p_value <= 0.05 from coxph on single gene.',
            'Dropped_genes: genes dropped from modeling with total risk value with high and low.',
            'Min_p_gene: genes dropped from each cycle of screening and the smallest p value.', sep ='\n'))
}


## Double click the below code to run.

selectlnr()

