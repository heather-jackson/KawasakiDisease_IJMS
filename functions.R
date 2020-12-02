packages.load <- function(){
  library(limma)
  library(gprofiler2)
  library(tidyverse)
  library(NbClust)
  library(factoextra)
  library(ggExtra)
  library(glmnet)
  library(ggplot2)
  library(stringr)
  library(pROC)
}

# function to run Limma for differential abundance analysis 
limma.res <- function(data, p.thresh, comparisons, start.col, end.col, n.features, model.mat){
  # construct linear models
  l.m <- lmFit(t(data[,start.col:end.col]), model.mat)
  # for pairwise comparisons between groups, make a contrast matrix
  constrasts.design <- paste(unlist(comparisons), collapse = "-")
  constrasts.mat <- makeContrasts(contrasts = constrasts.design, levels = colnames(model.mat))
  cb.fit <- eBayes(contrasts.fit(l.m, constrasts.mat))
  # list of SDE features
  top.features <- topTable(cb.fit, coef = 1, adjust = "BH", number = n.features, p.value = p.thresh)
  top.features$Gene.Name <- feat$ILMN_Gene[match(rownames(top.features), feat$ID)]
  return(top.features)
}

# function to run pathway analysis 
pathways <- function(feature.list){
  # only include significant genes
  significant.genes <- feature.list[feature.list$adj.P.Val < 0.05,]
  # identify those up or down expressed
  genes.up <- significant.genes[significant.genes$logFC > 0,]
  genes.down <- significant.genes[significant.genes$logFC < 0,]
  # run pathway analysis
  up.pathways <- gost(query = genes.up$Gene.Name,
                      organism = 'hsapiens',
                      correction_method = 'fdr',
                      ordered_query = T,
                      multi_query = F,
                      source = "GO:BP")
  up.pathways <- up.pathways$result
  down.pathways <- gost(query = genes.down$Gene.Name,
                        organism = 'hsapiens',
                        correction_method = 'fdr',
                        ordered_query = T,
                        multi_query = F,
                        source = "GO:BP")
  down.pathways <- down.pathways$result
  # filter results 
  up.pathways <- up.pathways[up.pathways$p_value < 0.01,]
  down.pathways <- down.pathways[down.pathways$p_value < 0.01,]
  return(list(up.pathways = up.pathways, down.pathways = down.pathways))
}

# function to read in the REVIGO files and filter them into lists of pathways 
top.pathways <- function(gp1.up, gp1.down, gp2.up, gp2.down, gp3.up, gp3.down, 
                         gp1, gp2, gp3){
  # filter the lists,  order, and take top 10
  gp1.up <- gp1.up[gp1.up$eliminated==0,]
  gp1.up <- gp1.up[order(gp1.up$log10.p.value),][1:10,]
  gp1.down <- gp1.down[gp1.down$eliminated==0,]
  gp1.down <- gp1.down[order(gp1.down$log10.p.value),][1:10,]
  gp2.up <- gp2.up[gp2.up$eliminated==0,]
  gp2.up <- gp2.up[order(gp2.up$log10.p.value),][1:10,]
  gp2.down <- gp2.down[gp2.down$eliminated==0,]
  gp2.down <- gp2.down[order(gp2.down$log10.p.value),][1:10,]
  gp3.up <- gp3.up[gp3.up$eliminated==0,]
  gp3.up <- gp3.up[order(gp3.up$log10.p.value),][1:10,]
  gp3.down <- gp3.down[gp3.down$eliminated==0,]
  gp3.down <- gp3.down[order(gp3.down$log10.p.value),][1:10,]
  # bind pathways together
  top.pathways <- bind_rows(data.frame(gp1.up,
                                       Disease = gp1,
                                       Direction = "Up"),
                            data.frame(gp1.down,
                                       Disease = gp1,
                                       Direction = "Down"),
                            data.frame(gp2.up,
                                       Disease = gp2,
                                       Direction = "Up"),
                            data.frame(gp2.down,
                                       Disease = gp2,
                                       Direction = "Down"),
                            data.frame(gp3.up,
                                       Disease = gp3,
                                       Direction = "Up"),
                            data.frame(gp3.down,
                                       Disease = gp3,
                                       Direction = "Down"))
  top.pathways$Direction.P <- -top.pathways$log10.p.value
  top.pathways$Direction.P[top.pathways$Direction == "Down"] <- -top.pathways$Direction.P[top.pathways$Direction == "Down"]
  top.pathways.w <- as.data.frame(pivot_wider(data = top.pathways,
                                              names_from = Disease,
                                              values_from = Direction.P,
                                              id_cols = description))
  # order so the most similar results are next to each other
  order <- arrange( top.pathways.w, 2, 3, 4)
  order <- order$description
  top.pathways.w <- top.pathways.w[match(order, top.pathways.w$description),]
  top.pathways.w <- unnest(top.pathways.w, cols = c(gp1, gp2, gp3))
  top.pathways.w <- pivot_longer(top.pathways.w,
                                 cols = c(gp1, gp2, gp3),
                                 names_to = "Disease")
  top.pathways.w$num <- nrow(top.pathways.w):1
  top.pathways.w$direction <- "Up"
  top.pathways.w$direction[top.pathways.w$value < 0] <- "Down"
  return(top.pathways.w)
}

# function for preparing the dataset ready for clustering:
# remove contribution of other variables, identify most variable genes, determine optimal number of clusters
clustering.preparation <- function(model.mat, exp, meta, correct){
  if(correct == TRUE){
    limma.correct <- lmFit(as.matrix(exp), design = model.mat[,1:7])
    df.corrected <- residuals.MArrayLM(obj = limma.correct, y = as.matrix(exp))
  }
  else if(correct == FALSE){
    df.corrected <- exp
  }
  # df.corrected is now the corrected data for all samples
  # subset it so it only contains HC samples for downstream analysis
  meta.hc <- meta[!(meta$Group == "Control"),]
  df.corrected.all <- df.corrected
  df.corrected <- df.corrected[,colnames(df.corrected) %in% meta.hc$ID]
  c.meta <-  meta[meta$Group == "Control",]
  # filter the dataset so it only contains the most variable genes 
  var.genes <- apply(df.corrected, 1, function(x){
    var(x)
  })
  # only include those with variance over 0.25
  most.var <- var.genes[var.genes > 0.25]
  e.set.variable <- df.corrected[rownames(df.corrected) %in% names(most.var),]
  rm(df.corrected)
  # run NbClust to determine the optimal number of clusters 
  methods <- c("kl", "ch", "hartigan","mcclain",
               "dunn", "sdindex", "sdbw", "cindex",
               "silhouette", "ball","ptbiserial", 'ratkowsky')
  # run this on the dataset that has been corrected
  table.corrected = as.data.frame(matrix(ncol = 1, nrow = length(methods)))
  # determine optimal K 
  for(i in 1:length(methods)){
    #print(i)
    nb.cor = NbClust(t(e.set.variable), distance = 'euclidean',
                     min.nc = 2, max.nc = 10,
                     method = "kmeans", index =methods[i])
    table.corrected[i,] = nb.cor$Best.nc[1]
    rownames(table.corrected)[i] = methods[i]}
  getmode <- function(v) {
    uniqv <- unique(v)
    uniqv[which.max(tabulate(match(v, uniqv)))]
  }
  mode <- getmode(table.corrected$V1)
  print(paste("Optimal number of clusters for corrected dataset is:", mode, sep = ""))
  
  return(list(clusters.corrected = mode, 
              e.set.variable = e.set.variable, 
              df.corrected.all = df.corrected.all))
}


# function for training the classifier 
classifier <- function(group.interest, limma.results, expression.set){
  # change the group labels to either 0 vs 1 
  # 1 is the group to classify from the rest 
  expression.set$group <- as.character(expression.set$group)
  expression.set$group[!(expression.set$group == group.interest)] <- 0
  expression.set$group[expression.set$group== group.interest] <- 1
  expression.set$group <- as.numeric(expression.set$group)
  # set seed for reproducible results
  set.seed(1234)
  # train classifier 
  cvfit <- cv.glmnet(x = as.matrix(expression.set[,2:ncol(expression.set)]),
                     y = expression.set$group,
                     family= 'binomial',
                     type.measure = "mae",
                     alpha = 1,
                     folds = 10)
  # extract coefficients 
  coefs.fit <- coef(cvfit, s = "lambda.1se")
  coefs.fit <- t(t(as.matrix(coefs.fit)[-1, ]))
  coefs.fit <- coefs.fit[!(coefs.fit[,1]==0),]
  # check the fold change direction to see whether coefficient direction is concordant with limma results 
  lasso.limma <- data.frame(LFC = limma.results$logFC[match(names(coefs.fit), rownames(limma.results))],
                            Lasso = coefs.fit)
  lasso.limma$Concordant <- "Y"
  lasso.limma$Concordant[ which( lasso.limma$LFC > 0 & lasso.limma$Lasso < 0) ] <- "N"
  lasso.limma$Concordant[ which( lasso.limma$LFC < 0 & lasso.limma$Lasso > 0) ] <- "N"
  # add gene name for convenience 
  lasso.limma$Name <- feat$Symbol[match(rownames(lasso.limma), feat$ID)]
  print(table(lasso.limma$Concordant))
  return(lasso.limma)
}

# function for testing the classifiers 
classifier.test <- function(classifier.res, exp.test, kd.exp, group.interest){
  down.genes <- rownames(classifier.res)[classifier.res$Lasso < 0]
  up.genes <- rownames(classifier.res)[classifier.res$Lasso > 0]
  down.genes.exp <- exp.test[,match(down.genes, colnames(exp.test))]
  up.genes.exp <- exp.test[,match(up.genes, colnames(exp.test))]
  # calculate disease risk score
  # subtract the total expression of negatively regulated genes from the positively regulated genes
  up.genes.exp.sum <- apply(up.genes.exp, 1, function(x){sum(as.numeric(as.character(x)))})
  down.genes.exp.sum <- apply(down.genes.exp, 1, function(x){sum(as.numeric(as.character(x)))})
  drs <- up.genes.exp.sum - down.genes.exp.sum
  # scale between 0-1
  drs <- (drs-min(drs))/
    (max(drs)-min(drs))
  # using the DRS calculated on DB and DV samples and the known DB/DV classification, calculate thresholds above/below which a sample is DB/DV
  if(group.interest == "DV"){
    roc.res <- roc(fct_rev(droplevels(as.factor(exp.test$group))), drs, plot = T,ci = T, print.auc = T)
  }
  else if (group.interest == "DB"){
    roc.res <- roc(droplevels(as.factor(exp.test$group)), drs, plot = T,ci = T, print.auc = T)
  }
  my.coords <- coords(roc=roc.res, x = "all", transpose = FALSE)
  my.coords <- my.coords[my.coords$sensitivity >= .90, ]
  thresh <- my.coords[which(my.coords$sensitivity==min(my.coords$sensitivity)),]
  thresh <- thresh$threshold[thresh$specificity == max(thresh$specificity)]
  # calculate misclassification rate
  mis.class <- data.frame(DRS = drs,
                          truth = exp.test$group)
  mis.class$class <- ""
  mis.class$class[mis.class$DRS > thresh] <- group.interest
  print(table(mis.class$class, mis.class$truth))
  
  # now test the classifier on the KD samples 
  kd.exp.down <- kd.exp[rownames(kd.exp) %in% down.genes,]
  kd.exp.up <- kd.exp[rownames(kd.exp) %in% up.genes,]
  # calculate DRS 
  kd.up.sum <- apply(kd.exp.up, 2, function(x){sum(as.numeric(as.character(x)))})
  kd.down.sum <- apply(kd.exp.down, 2, function(x){sum(as.numeric(as.character(x)))})
  drs.kd <- kd.up.sum - kd.down.sum
  drs.kd <- (drs.kd-min(drs.kd))/
    (max(drs.kd)-min(drs.kd))
  #print(drs.kd)
  class.kd <- data.frame(DRS = drs.kd)
  class.kd$class <- paste('Not_', group.interest, sep = "")
  class.kd$class[class.kd$DRS > thresh] <- group.interest
  return(list(test.classes = mis.class, 
              kd.classes = class.kd, 
              threshold = thresh))
}
