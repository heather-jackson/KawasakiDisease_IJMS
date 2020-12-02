####################################################################

# load normalised data in addition to metadata and table to translate probe IDs to illumina IDs
# order of samples must match up
setwd('KawasakiDisease_IJMS/')

load(file = 'E-set-Trans.RData', verbose = T)
# e.set = normalised dataset 
# meta = metadata
# feat = feature information to match Illumina probe IDs to gene labels

# data has been normalised using RSN normalisation using the Lumi package


# source functions
source(file = 'functions.R')
packages.load()

####################################################################
# 1) Differential abundance analysis - 2.2.1
####################################################################

# create model matrix 
# cell proportions have been estimated using CibersortX
model.mat <- model.matrix(~0 + factor(Group) +
                           factor(Sex) +
                           as.numeric(age) +
                           as.numeric(Lymphocytes) +
                           as.numeric(Monocyte) +
                           as.numeric(Mast)+
                           as.numeric(Neutrophils),
                         data = meta)

colnames(model.mat) <- c("bacterial", "Control", "KD", "viral", "M", "age",
                        "Lymphocytes", "Monocyte", "Mast", "Neutrophils")

# rename features
feat <- feat[match(rownames(e.set), feat$Array_Address_Id),]
rownames(e.set) <- feat$ID


# identify genes differentially abundant between KD vs Controls
genes.sde.kd.control <- limma.res(data = t(e.set), 
                                  p.thresh = 1,
                                  comparisons = c("KD", "Control"),
                                  start.col = 1, 
                                  end.col = nrow(e.set),
                                  n.features = nrow(e.set), 
                                  model.mat = model.mat)

# identify genes differentially abundant between DB vs Controls
genes.sde.db.control <- limma.res(data = t(e.set),
                                  p.thresh = 1,
                                  comparisons = c("bacterial", "Control"),
                                  start.col = 1, 
                                  end.col = nrow(e.set),
                                  n.features = nrow(e.set), 
                                  model.mat = model.mat)

# identify genes differentially abundant between DV vs Controls
genes.sde.dv.control <- limma.res(data = t(e.set), 
                                  p.thresh = 1,
                                  comparisons = c("viral", "Control"),
                                  start.col = 1, 
                                  end.col = nrow(e.set),
                                  n.features = nrow(e.set), 
                                  model.mat = model.mat)

####################################################################
# 2) Pathway analysis - 2.2.2
####################################################################

kd.paths.up <- pathways(genes.sde.kd.control)[[1]]
kd.paths.down <- pathways(genes.sde.kd.control)[[2]]
db.paths.up <- pathways(genes.sde.db.control)[[1]]
db.paths.down <- pathways(genes.sde.db.control)[[2]]
dv.paths.up <- pathways(genes.sde.dv.control)[[1]]
dv.paths.down <- pathways(genes.sde.dv.control)[[2]]

# the lists of pathways must be input into the REVIGO web server to reduce redundancy
# the 'medium' parameter was selected in REVIGO to determine the extend of filtering of significant pathways 
# note these lists are available in the supplementary materials
# save each file as a CSV 

# read in each file 
# update /path/to/file.csv with your file name and path
kd.up.revigo <- read.csv(file = '/path/to/file.csv', header = T)
kd.down.revigo <- read.csv(file = '/path/to/file.csv', header = T)
db.up.revigo <- read.csv(file = '/path/to/file.csv', header = T)
db.down.revigo <- read.csv(file = '/path/to/file.csv', header = T)
dv.up.revigo <- read.csv(file = '/path/to/file.csv', header = T)
dv.down.revigo <- read.csv(file = '/path/to/file.csv', header = T)

# condense the pathway analysis results into one dataframe with top 10 results for each run of revigo
top.paths <- top.pathways(gp1.up = kd.up.revigo, 
                          gp1.down = kd.down.revigo, 
                          gp2.up = db.up.revigo, 
                          gp2.down = db.down.revigo, 
                          gp3.up = dv.up.revigo,
                          gp3.down = dv.down.revigo, 
                          gp1 = 'KD',
                          gp2 = "DB", 
                          gp3 = "DV")

# some manual adaptations may be needed
# i.e. cell activation is up and down regulated in the KD group
top.paths$direction[c(7,10)]<- "Both"

ggplot(top.paths[!(is.na(top.paths$value)),],
       aes(y = fct_reorder(description, num),
           x = fct_relevel(Disease, "DV", "KD", "DB"),
           size = abs(value),
           color = fct_rev(direction)))+
  theme_bw()+
  geom_point()+
  labs(size = "-log10 P-value",
       color = "Direction",
       x = "",
       y = "")+
  theme(axis.text.x = element_text(angle =35, hjust = 1),
        plot.margin = margin(10, 10, 10, 10),
        text = element_text(size = 12))+
  scale_color_manual(breaks = c("Up", "Down", "Both"), values = c("#d73027", "#313695", "#fee090"))+
  scale_size_continuous(range = c(2,6))+
  scale_x_discrete(labels = c("Definite Viral vs Healthy Controls", "Kawasaki Disease vs Healthy Controls",
                              "Definite Bacterial vs Healthy Controls" ))+
  guides(colour = guide_legend(override.aes = list(size=4)))

####################################################################
# 3) Clustering analysis - 2.2.3
####################################################################

# before running clustering analysis, first remove the contributions of age, sex and immune cell proportions (estimated by CIBERSORTx)
# isolate these variables
cells <- meta[,c(4, 5, 10:14)]
rownames(cells) <- meta$ID
cells <- cells[match(colnames(e.set), rownames(cells)),]

# create model matrix from the variables we are correcting for 
correct.matrix <- model.matrix(~0 + as.numeric(age) +
                                 as.numeric(Lymphocytes) +
                                 as.numeric(Monocyte)+
                                 as.numeric(Mast)+
                                 as.numeric(Eosinophils)+
                                 as.numeric(Neutrophils)+
                                 as.factor(Sex),
                               data = cells)
colnames(correct.matrix) <- c('age', 'lymphocytes', "monocyte", "mast", "eosinophil", "neutrophil", "female", "male")

# run the function to prepare the data and identify optimal number of clusters
clusters.prep <- clustering.preparation(correct.matrix, e.set, meta, correct = TRUE)

# run K-Means on the corrected dataset using the number of clusters determined by NbClust
set.seed(2020)
clusters <- kmeans(t(clusters.prep$e.set.variable),
                             centers = clusters.prep$clusters.corrected,
                             nstart = 15,
                             iter.max = 15)

# create plot with clusters split up by disease group 
clusters <- data.frame(clusters = clusters$cluster, 
                       group = meta$Group[match(names(clusters$cluster), meta$ID)])

ggplot(clusters, aes(y = fct_rev(as.factor(clusters)), 
                       fill = group))+
  geom_bar(position = 'fill')+
  theme_bw()+
  scale_x_continuous(breaks=seq(0,1,0.2), 
                     labels = c("0%", "20%", '40%', '60%',  '80%',  '100%'))+
  scale_fill_grey(start = 0.8, end = 0.3)+
  labs(x = "Proportion of each cluster", 
       fill = "Disease groups", 
       y = "")+
  scale_y_discrete(labels = c("Cluster 3", "Cluster 2", "Cluster 1"))

# run pathway analysis on the KD samples that fall in the different clusters 
# first of all, need to identify genes significantly differentially expressed 
clusters.kd <- clusters[clusters$group == "KD", ]
rm(clusters)

e.set.variable.kd <- clusters.prep$e.set.variable[,match(rownames(clusters.kd), colnames(clusters.prep$e.set.variable))]
model.mat.clusters <- model.matrix(~0 + as.factor(clusters), data = clusters.kd)
colnames(model.mat.clusters) <- c('Cl.1', 'Cl.2', 'Cl.3')


# run limma 
# compare cluster 1 to clusters 2 and 3
cluster1.sde.genes <- limma.res(data = t(e.set.variable.kd), 
          p.thresh = 0.05,
          comparisons = c("Cl.1"),
          start.col = 1, 
          end.col = nrow(e.set.variable.kd),
          n.features = nrow(e.set.variable.kd), 
          model.mat = model.mat.clusters)

# compare cluster 2 to clusters 1 and 3
cluster2.sde.genes <- limma.res(data = t(e.set.variable.kd),
                       p.thresh = 0.05,
                       comparisons = c('Cl.2'),
                       start.col = 1,
                       end.col = nrow(e.set.variable.kd),
                       n.features = nrow(e.set.variable.kd),
                       model.mat = model.mat.clusters)

# compare cluster 3 to clusters 2 and 3
cluster3.sde.genes <- limma.res(data = t(e.set.variable.kd),
                       p.thresh = 0.05,
                       comparisons = c('Cl.3'),
                       start.col = 1,
                       end.col = nrow(e.set.variable.kd),
                       n.features = nrow(e.set.variable.kd),
                       model.mat = model.mat.clusters)

# run pathway analysis
cl.paths.up <- pathways(cluster1.sde.genes)[[1]]
cl.paths.down <- pathways(cluster1.sde.genes)[[2]]
c2.paths.up <- pathways(cluster2.sde.genes)[[1]]
c2.paths.down <- pathways(cluster2.sde.genes)[[2]]
c3.paths.up <- pathways(cluster3.sde.genes)[[1]]
c3.paths.down <- pathways(cluster3.sde.genes)[[2]]

# the lists of pathways must be input into the REVIGO web server to reduce redundancy
# the 'medium' parameter was selected in REVIGO to determine the extend of filtering of significant pathways 
# note these lists are available in the supplementary materials
# save each file as a CSV 
# read in each file 

# update /path/to/file.csv with your file name and path
c1.up.revigo <- read.csv(file = '/path/to/file.csv', header = T)
c1.down.revigo <- read.csv(file = '/path/to/file.csv', header = T)
c2.up.revigo <- read.csv(file = '/path/to/file.csv', header = T)
c2.down.revigo <- read.csv(file = '/path/to/file.csv', header = T)
c3.up.revigo <- read.csv(file = '/path/to/file.csv', header = T)
c3.down.revigo <- read.csv(file = '/path/to/file.csv', header = T)

# identify top 10 pathways up and down regulated in each cluster
cluster.top.pathways <- top.pathways(gp1.up = c1.up.revigo, 
             gp1.down = c1.down.revigo, 
             gp2.up = c2.up.revigo, 
             gp2.down = c2.down.revigo, 
             gp3.up = c3.up.revigo,
             gp3.down = c3.down.revigo, 
             gp1 = 'Cluster1',
             gp2 = "Cluster2", 
             gp3 = "Cluster3")

# again, some may need to be manually edited if the pathway is up and down regulated in the same group
cluster.top.pathways$direction[c(12, 15, 18, 21, 24, 27)] <- 'Both'

ggplot(cluster.top.pathways[!is.na(cluster.top.pathways$value),], 
       aes(y = fct_reorder(description, num),
           x = as.factor(Disease),
           size = abs(value),
           color = direction))+
  theme_bw()+
  geom_point()+
  labs(x = "", y = "", color = "Direction", size = "-log10 p.value")+
  theme(axis.text.x = element_text(angle = 35, hjust = 1))+
  scale_color_manual(breaks = c("Up", "Down", "Both"), values = c("#d73027", "#313695", "#fee090"))+
  theme(axis.text.x = element_text(angle = 35, hjust = 1),
        plot.margin = margin(10, 10, 10, 10),
        text = element_text(size = 12))+
  scale_size_continuous(range = c(2,6))+
  scale_x_discrete(labels = c("KD in cluster 1",
                              "KD in cluster 2",
                              "KD in cluster 3"))+
  guides(colour = guide_legend(override.aes = list(size=4)))

####################################################################
# 4) Classification analysis - 2.2.4
####################################################################

# remove any bacterial samples that might not be pure, i.e. those with known coinfections
not.pure <- c('bacterialgneg_8_SMH',
              'bacterialgpos_15_SMH',
              'bacterialgneg_5_SMH',
              'bacterialgpos_13_SMH',
              'bacterialgneg_18_SMH',
              'bacterialgneg_1_SMH',
              'bacterialgpos_23_SOT',
              'bacterialgpos_14_SMH',
              'bacterialgpos_10_SMH',
              'bacterialgpos_12_SMH',
              'bacterialgneg_11_SMH',
              'bacterialgneg_2_SMH')

meta <- meta[!(meta$ID) %in% not.pure,]

# use only genes SDE between DB vs (HC+DV) and DV vs (DB+HC)
# run limma to identify these genes 
# subset the data to remove KD 
meta.bv.hc <- meta[!(meta$Group=="KD"),]
meta.bv.hc$DB <- "Yes"
meta.bv.hc$DB[!(meta.bv.hc$Group == "bacterial")] <- "No"
meta.bv.hc$DV <- "Yes"
meta.bv.hc$DV[!(meta.bv.hc$Group == "viral")] <- "No"

# we are using the NON corrected data to identify the SDE genes, therefore need to correct for age, sex and immune cells 
# create the model matrices to use for limma
model.mat.db <- model.matrix(~0 + factor(DB) +
                                     factor(Sex) +
                                     as.numeric(age) +
                                     as.numeric(Lymphocytes) +
                                     as.numeric(Monocyte) +
                                     as.numeric(Mast)+
                                     as.numeric(Neutrophils),
                                   data = meta.bv.hc)

colnames(model.mat.db) <- c("DB.No", "DB.Yes", "M", "age",
                                  "Lymphocytes", "Monocyte", "Mast", "Neutrophils")
rownames(model.mat.db) <- meta.bv.hc$ID

model.mat.dv <- model.matrix(~0 + factor(DV) +
                                     factor(Sex) +
                                     as.numeric(age) +
                                     as.numeric(Lymphocytes) +
                                     as.numeric(Monocyte) +
                                     as.numeric(Mast)+
                                     as.numeric(Neutrophils),
                                   data = meta.bv.hc)

colnames(model.mat.dv) <- c("DV.No", "DV.Yes", "M", "age",
                                  "Lymphocytes", "Monocyte", "Mast", "Neutrophils")
rownames(model.mat.dv) <- meta.bv.hc$ID

# filter expression set to remove KD 
e.set.bv.hc <- e.set[,match(meta.bv.hc$ID, colnames(e.set))]
model.mat.db <- model.mat.db[match(colnames(e.set.bv.hc), rownames(model.mat.db)),]
model.mat.dv <- model.mat.dv[match(colnames(e.set.bv.hc), rownames(model.mat.dv)),]

# identify genes SDE between DB vs (DV+HC)
db.sde.genes <- limma.res(data = t(e.set.bv.hc),
                                p.thresh = 0.01,
                                comparisons =  c("DB.Yes", "DB.No"),
                                start.col = 1,
                                end.col = nrow(e.set.bv.hc),
                                n.features = nrow(e.set.bv.hc),
                                model.mat = model.mat.db)

# identify genes SDE between DB vs (DV+HC)
dv.sde.genes <- limma.res(data = t(e.set.bv.hc),
                          p.thresh = 0.01,
                          comparisons =  c("DV.Yes", "DV.No"),
                          start.col = 1,
                          end.col = nrow(e.set.bv.hc),
                          n.features = nrow(e.set.bv.hc),
                          model.mat = model.mat.dv)


# create 2 new data frames with these features but this time using the corrected dataset
exp.db <- clusters.prep$df.corrected.all[rownames(clusters.prep$df.corrected.all) %in% rownames(db.sde.genes),
                                         match(meta.bv.hc$ID, colnames(clusters.prep$df.corrected.all))]
exp.dv <- clusters.prep$df.corrected.all[rownames(clusters.prep$df.corrected.all) %in% rownames(dv.sde.genes),
                                         match(meta.bv.hc$ID, colnames(clusters.prep$df.corrected.all))]
# extract the KD samples with the same set of features- these are to test the classifiers on 
exp.db.kd <- clusters.prep$df.corrected.all[rownames(clusters.prep$df.corrected.all) %in% rownames(db.sde.genes),
                                            !(colnames(clusters.prep$df.corrected.all) %in% meta.bv.hc$ID) & 
                                              !(colnames(clusters.prep$df.corrected.all) %in% not.pure)]
exp.dv.kd <- clusters.prep$df.corrected.all[rownames(clusters.prep$df.corrected.all) %in% rownames(dv.sde.genes),
                                            !(colnames(clusters.prep$df.corrected.all) %in% meta.bv.hc$ID) & 
                                              !(colnames(clusters.prep$df.corrected.all) %in% not.pure)]

# add on class labels to traning data 
# this is the training data for the bacterial classifier 
exp.db <- data.frame(Group = as.factor(meta.bv.hc$Group), 
                     t(exp.db))
# this is the training data for the viral classifier 
exp.dv <- data.frame(Group = as.factor(meta.bv.hc$Group), 
                     t(exp.dv))

# read in the validation dataset used to test the classifier performance 
load(file = 'validation_transcriptomics.RData', 
     verbose = T)


# remove the contribution of age, sex and immune cell proportions (determined using CIBERSORTx)
model.mat.val  <- model.matrix(~0 + as.numeric(age) +
                                     as.numeric(Lymphocytes) +
                                     as.numeric(Monocyte)+
                                     as.numeric(Mast)+
                                     as.numeric(Eosinophils)+
                                     as.numeric(Neutrophils)+
                                     as.factor(sex),
                                   data = meta.val)
colnames(model.mat.val) <- c('age', 'lymphocytes', "monocyte", "mast", "eosinophil", "neutrophil", "female", "male")

# model residuals of variables
limma.correct.val <- lmFit(as.matrix(e.set.val), design = model.mat.val[,1:7])
df.corrected.val <- residuals.MArrayLM(obj = limma.correct.val, y = as.matrix(e.set.val))
rm(limma.correct.val, e.set.val)

# ensure the same set of genes are present in the discovery and validation datasets for each classifier 
db.intersect <- intersect(rownames(df.corrected.val), colnames(exp.db))
dv.intersect <- intersect(rownames(df.corrected.val), colnames(exp.dv))

exp.db.val <- df.corrected.val[match(db.intersect, rownames(df.corrected.val)),]
exp.dv.val <- df.corrected.val[match(dv.intersect, rownames(df.corrected.val)),]
exp.db <- data.frame(group= exp.db$Group, 
                     exp.db[,match(db.intersect, colnames(exp.db))])
exp.dv <- data.frame(group = exp.dv$Group, 
                     exp.dv[,match(dv.intersect, colnames(exp.dv))])


# and the KD test discovery
exp.db.kd <- exp.db.kd[match(db.intersect, rownames(exp.db.kd)),]
exp.dv.kd <- exp.dv.kd[match(dv.intersect, rownames(exp.dv.kd)),]

exp.db.val <- data.frame(group = meta.val$group[match(str_remove_all(colnames(exp.db.val), "X"), meta.val$ID)],
                         t(exp.db.val))
exp.dv.val <- data.frame(group = meta.val$group[match(str_remove_all(colnames(exp.dv.val), "X"), meta.val$ID)],
                         t(exp.dv.val))
# extract the KD samples for one test dataset
exp.db.val_KD <- exp.db.val[exp.db.val$group == "Acute", ]
exp.dv.val_KD <- exp.dv.val[exp.dv.val$group == "Acute", ]

# keep only bacterial and viral samples in another test dataset
exp.db.val <- exp.db.val[exp.db.val$group == "bacterial" |
                           exp.db.val$group == "viral",]
exp.dv.val <- exp.dv.val[exp.dv.val$group == "bacterial" |
                           exp.dv.val$group == "viral",]


### train the classifiers 
# bacterial 
bact.classifier <- classifier('bacterial', db.sde.genes, exp.db)
# viral 
vir.classifier <- classifier('viral', dv.sde.genes, exp.dv)

# check concordance
table(bact.classifier$Concordant)
table(vir.classifier$Concordant)

# test the classifiers, first of all on the bacterial and viral samples from the validation dataset 
# this will return thresholds that can be applied to the KD samples 
# first of all make one large dataframe with all kawasaki samples included 

##### now test on the KD samples
# make one big df for KD samples
exp.db.val_KD <- exp.db.val_KD[,match(rownames(exp.db.kd), colnames(exp.db.val_KD))]
exp.dv.val_KD <- exp.dv.val_KD[,match(rownames(exp.dv.kd), colnames(exp.dv.val_KD))]

all.kd.db <- data.frame(exp.db.kd,
                        t(exp.db.val_KD))
all.kd.dv <- data.frame(exp.dv.kd,
                        t(exp.dv.val_KD))

# test the classifiers on the DB/DV datasets and also the merged KD dataset 
db.classifier.res <- classifier.test(bact.classifier, exp.db.val, all.kd.db, 'DB')
dv.classifier.res <- classifier.test(vir.classifier, exp.dv.val, all.kd.dv, 'DV')

# join the results 
classifier.results.kd <- data.frame(DB.DRS = db.classifier.res$kd.classes, 
                                    DV.DRS = dv.classifier.res$kd.classes)

classifier.results.kd$class.joined <- NA
classifier.results.kd$class.joined[classifier.results.kd$DB.DRS.class == "DB"] <- "DB"
classifier.results.kd$class.joined[classifier.results.kd$DV.DRS.class == "DV"] <- "DV"

classifier.results.kd$class.joined[classifier.results.kd$DB.DRS.class == "DB" & 
                                     classifier.results.kd$DV.DRS.class == "DV"] <- "Both"
classifier.results.kd$class.joined[classifier.results.kd$DB.DRS.class == "Not_DB" & 
                                     classifier.results.kd$DV.DRS.class == "Not_DV"] <- "Neither"

# How many were classified in each group?
table(classifier.results.kd$class.joined)

# join the results with the DB/DV test dataset to visualisse
all.results.classifier <- rbind(data.frame(DB.DRS = db.classifier.res$test.classes$DRS, 
                     DV.DRS = dv.classifier.res$test.classes$DRS, 
                     Truth = db.classifier.res$test.classes$truth), 
                data.frame(DB.DRS = classifier.results.kd$DB.DRS.DRS, 
                           DV.DRS = classifier.results.kd$DV.DRS.DRS, 
                           Truth = "KD"))

# plot 
p <- ggplot(all.results.classifier,
                aes(x = DB.DRS,
                    y = DV.DRS))+
  geom_jitter(aes(col=as.factor(Truth)),
              size = 2)+
  theme_bw()+
  labs(color = "Disease Group",
       y = "Viral DRS",
       x = "Bacterial DRS")+
  scale_color_manual(breaks = c("bacterial", "KD", "viral"),
                     values=c("#d73027",'darkgrey',  "#313695"),
                     labels = c("DB", "KD", "DV"))+
  geom_vline(xintercept = db.classifier.res$threshold, linetype = 'dashed')+
  geom_hline(yintercept = dv.classifier.res$threshold, linetype = 'dashed')

ggMarginal(p,
           type = "boxplot",
           groupFill = T,
           size = 7,
           groupColour = T)

####################################################################
# 5) Clustering analysis with KD - 2.2.4
####################################################################

# remove the KD samples from the corrected expression set 
meta.kd <- meta[match(colnames(e.set.variable.kd), meta$ID),]
correct.matrix.kd <- correct.matrix[match(meta.kd$ID, rownames(correct.matrix)),]
# remove the contribution of age, sex and immune cell proportions (estimated using CIBERSORTx)
# determine the optimal number of clusters
clusters.prep.kd <- clustering.preparation(correct.matrix.kd, e.set.variable.kd, meta.kd, correct = FALSE)


# run K-Means on the corrected dataset using the number of clusters determined by NbClust
set.seed(2020)
clusters.kd.alone <- kmeans(t(clusters.prep.kd$e.set.variable),
                   centers = clusters.prep.kd$clusters.corrected,
                   nstart = 15,
                   iter.max = 15)

# run pathway analysis on the samples that fall in the different clusters 
# first of all, need to identify genes significantly differentially expressed 
model.mat.clusters.kd <- model.matrix(~0 + as.factor(clusters.kd.alone$cluster))
colnames(model.mat.clusters.kd) <- c('Cl.1', 'Cl.2', 'Cl.3')

# run limma 
# compare cluster 1 to clusters 2 and 3
cluster1.sde.genes.kd <- limma.res(data = t(clusters.prep.kd$e.set.variable), 
                                p.thresh = 0.05,
                                comparisons = c("Cl.1"),
                                start.col = 1, 
                                end.col = nrow(clusters.prep.kd$e.set.variable),
                                n.features = nrow(clusters.prep.kd$e.set.variable), 
                                model.mat = model.mat.clusters.kd)

# compare cluster 2 to clusters 1 and 3
cluster2.sde.genes.kd <- limma.res(data = t(clusters.prep.kd$e.set.variable),
                                p.thresh = 0.05,
                                comparisons = c('Cl.2'),
                                start.col = 1,
                                end.col = nrow(clusters.prep.kd$e.set.variable),
                                n.features = nrow(clusters.prep.kd$e.set.variable),
                                model.mat = model.mat.clusters.kd)

# compare cluster 3 to clusters 2 and 3
cluster3.sde.genes.kd <- limma.res(data = t(clusters.prep.kd$e.set.variable),
                                p.thresh = 0.05,
                                comparisons = c('Cl.3'),
                                start.col = 1,
                                end.col = nrow(clusters.prep.kd$e.set.variable),
                                n.features = nrow(clusters.prep.kd$e.set.variable),
                                model.mat = model.mat.clusters.kd)

# run pathway analysis
cl.paths.up.kd <- pathways(cluster1.sde.genes.kd)[[1]]
cl.paths.down.kd <- pathways(cluster1.sde.genes.kd)[[2]]
c2.paths.up.kd <- pathways(cluster2.sde.genes.kd)[[1]]
c2.paths.down.kd <- pathways(cluster2.sde.genes.kd)[[2]]
c3.paths.up.kd <- pathways(cluster3.sde.genes.kd)[[1]]
c3.paths.down.kd <- pathways(cluster3.sde.genes.kd)[[2]]

# the lists of pathways must be input into the REVIGO web server to reduce redundancy
# the 'medium' parameter was selected in REVIGO to determine the extend of filtering of significant pathways 
# note these lists are available in the supplementary materials
# save each file as a CSV 
# read in each file 
# replace path/to/file.csv with the path to your files from REVIGO
c1.up.revigo <- read.csv(file = 'path/to/file.csv', header = T)
c1.down.revigo <- read.csv(file = 'path/to/file.csv', header = T)
c2.up.revigo <- read.csv(file = 'path/to/file.csv', header = T)
c2.down.revigo <- read.csv(file = 'path/to/file.csv', header = T)
c3.up.revigo <- read.csv(file = 'path/to/file.csv', header = T)
c3.down.revigo <- read.csv(file = 'path/to/file.csv', header = T)

# identify top 10 pathways up and down regulated in each cluster
cluster.top.pathways.kd <- top.pathways(gp1.up = c1.up.revigo, 
                                     gp1.down = c1.down.revigo, 
                                     gp2.up = c2.up.revigo, 
                                     gp2.down = c2.down.revigo, 
                                     gp3.up = c3.up.revigo,
                                     gp3.down = c3.down.revigo, 
                                     gp1 = 'Cluster1',
                                     gp2 = "Cluster2", 
                                     gp3 = "Cluster3")

# again, some may need to be manually edited if the pathway is up and down regulated in the same group
cluster.top.pathways.kd$direction[c(22, 25, 32, 35)] <- 'Both'

ggplot(cluster.top.pathways.kd[!is.na(cluster.top.pathways.kd$value),], 
       aes(y = fct_reorder(description, num),
           x = as.factor(Disease),
           size = abs(value),
           color = direction))+
  theme_bw()+
  geom_point()+
  labs(x = "", y = "", color = "Direction", size = "-log10 p.value")+
  theme(axis.text.x = element_text(angle = 35, hjust = 1))+
  scale_color_manual(breaks = c("Up", "Down", "Both"), values = c("#d73027", "#313695", "#fee090"))+
  theme(axis.text.x = element_text(angle = 35, hjust = 1),
        plot.margin = margin(10, 10, 10, 10),
        text = element_text(size = 12))+
  scale_size_continuous(range = c(2,6))+
  scale_x_discrete(labels = c("Cluster 1",
                              "Cluster 2",
                              "Cluster 3"))+
  guides(colour = guide_legend(override.aes = list(size=4)))