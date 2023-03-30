#This script is for analysis of the connectivity matrices generated in MatLab for 
#each clinical subgroup. It also creates a 'stats_table' dataframe with graph measures 
#and clinical information to then perform linear regression.

#---------------------------------------------------------------------------------------
library(DescTools)
library(R.matlab)
library(dgof)
library(stringr)
library(ggplot2)
library(ggpubr)
library(readxl)
library(ggsignif)

df <- read.table('./anxiety_0_thresholded.txt', sep = ",")
x <- df[upper.tri(df)]
df <- read.table('./anxiety_1_thresholded.txt', sep = ",")
y <- df[upper.tri(df)]
ks.test(x,y)

df <- read.table('./mania_0_thresholded.txt', sep = ",")
x <- df[upper.tri(df)]
df <- read.table('./mania_1_thresholded.txt', sep = ",")
y <- df[upper.tri(df)]
ks.test(x,y)

df <- read.table('./mania_0_thresholded.txt', sep = ",")
x <- df[upper.tri(df)]
df <- read.table('./mania_2_thresholded.txt', sep = ",")
y <- df[upper.tri(df)]
ks.test(x,y)

df <- read.table('./mania_0_thresholded.txt', sep = ",")
x <- df[upper.tri(df)]
df <- read.table('./mania_3_thresholded.txt', sep = ",")
y <- df[upper.tri(df)]
ks.test(x,y)

df <- read.table('./mania_1_thresholded.txt', sep = ",")
x <- df[upper.tri(df)]
df <- read.table('./mania_2_thresholded.txt', sep = ",")
y <- df[upper.tri(df)]
ks.test(x,y)

df <- read.table('./mania_1_thresholded.txt', sep = ",")
x <- df[upper.tri(df)]
df <- read.table('./mania_3_thresholded.txt', sep = ",")
y <- df[upper.tri(df)]
ks.test(x,y)

df <- read.table('./mania_2_thresholded.txt', sep = ",")
x <- df[upper.tri(df)]
df <- read.table('./mania_3_thresholded.txt', sep = ",")
y <- df[upper.tri(df)]
ks.test(x,y)

#----------------------------------Network measures---------------------------------------

data <- readMat('YMRS_RUN10_OCT_22_structure.mat')

subjects <- as.data.frame(data$subjects)

motions <- read.csv(file = 'motion_parameters.txt')
motions_subject <- as.data.frame(motions[,1])
motions_subject[] <- lapply(motions_subject, as.character)
motions[,1] <- motions_subject

stats_table <- matrix(nrow = nrow(subjects), ncol = 15)
stats_table <- as.data.frame(stats_table)
colnames(stats_table) <- c("subject", "age","sex", "cluster", "strength", "cpl", "transitivity", "assortativity", "clustering_coefficient", "participation_coeff", "mean_con", "framewise_displacement", "medication", "medication_binary", "dataset")
stats_table$subject <- data$subjects
stats_table$subject[] <- lapply(stats_table$subject, as.character)
class(stats_table$subject) <- "character"
stats_table$cluster <- data$clusterlabels
stats_table$cluster[] <- lapply(stats_table$cluster, as.character)
class(stats_table$cluster) <- "character"
strength <- rowMeans(data$strength)
stats_table$strength <- strength
stats_table$cpl <- data$cpl
stats_table$transitivity <- data$trans
stats_table$assortativity <- data$assor
stats_table$clustering_coefficient <- data$clustcoeff
part_coeff <- rowMeans(data$part)
stats_table$participation_coeff <- part_coeff
stats_table$mean_con <- data$meancon

for (i in 1:nrow(motions))
{
  motions[i,1] <- str_remove(as.character(motions[i,1]), "sub-")
}

for (subject in stats_table[,1])
{
  c <- stats_table[,1] == subject
  d <- motions[,1] == subject
  A <- as.numeric(motions[d,2])
  stats_table[c,12] <- A
}

demographic <- read.table(file = 'ymrs_all_FINAL_SELECTED.txt', sep = "|")
demo_subject <- as.data.frame(demographic[,1])
demo_subject[] <- lapply(demo_subject, as.character)
demographic[,1] <- demo_subject

for (subject in stats_table[,1])
{
  c <- stats_table[,1] == subject
  d <- demographic[,1] == subject
  A <- (as.numeric(demographic[d,2])/12)
  B <- as.character(demographic[d,3])
  stats_table[c,2] <- A
  stats_table[c,3] <- B
}

meds <- read.csv(file = 'meds_BD_HC.csv', header = TRUE)
meds <- meds[,2:3]
meds_subject <- as.data.frame(meds[,1])
meds_subject[] <- lapply(meds_subject, as.character)
meds[,1] <- meds_subject

for (subject in stats_table[,1])
{
  c <- stats_table[,1] == subject
  d <- meds[,1] == subject
  A <- as.character(meds[d,2])
  stats_table[c,13] <- A
}

meds_bin <- read.csv(file = 'meds_bin.csv', header = TRUE)
meds_subject <- as.data.frame(meds_bin[,1])
meds_subject[] <- lapply(meds_subject, as.character)
meds_bin[,1] <- meds_subject

for (subject in stats_table[,1])
{
  c <- stats_table[,1] == subject
  d <- meds_bin[,1] == subject
  A <- as.character(meds_bin[d,2])
  stats_table[c,14] <- A
}

#import image03_fmri file
dataset_info <- read_xlsx('image03_rsfmri.xlsx')
dataset_info2 <- dataset_info

for (i in 1:nrow(dataset_info))
{
  if ((dataset_info[i,5] %in% stats_table[,1]) == F)
  {
    dataset_info[i,] <- NA
  }
}

dup_subs <- duplicated(dataset_info$src_subject_id)
sum(dup_subs) #if sum is not zero:

# remove duplicates
dup_subs <- as.data.frame(dup_subs)

for (i in 1:nrow(dup_subs))
{
  if (dup_subs[i,] == TRUE)
  {
    dataset_info[i,] <- NA
  }
}

# remove NAs
dataset_info <- dataset_info[rowSums(is.na(dataset_info)) != ncol(dataset_info), ]
dataset_info <- as.data.frame(dataset_info)

for (subject in stats_table[,1])
{
  c <- stats_table[,1] == subject
  d <- dataset_info[,5] == subject
  A <- dataset_info[d,3]
  stats_table[c,15] <- A
}

#remove subjects with abnormally high connectivity values
for (i in 1:nrow(stats_table))
{
  if (stats_table[i,11] > 0.8)
  {
    stats_table[i,] <- NA
  }
}

# remove NAs
stats_table <- na.omit(stats_table)

#--------------------------------------add healthy controls--------------------------------------
data_HC <- readMat('HC_RUN10_OCT_22_structure.mat')

subjects_HC <- as.data.frame(data_HC$subjects)

stats_table_HC <- matrix(nrow = nrow(subjects_HC), ncol = 15)
stats_table_HC <- as.data.frame(stats_table_HC)
colnames(stats_table_HC) <- c("subject", "age","sex", "cluster", "strength", "cpl", "transitivity", "assortativity", "clustering_coefficient", "participation_coeff", "mean_con", "framewise_displacement", "medication", "medication_binary", "dataset")
stats_table_HC$subject <- data_HC$subjects
stats_table_HC$subject[] <- lapply(stats_table_HC$subject, as.character)
class(stats_table_HC$subject) <- "character"
stats_table_HC$cluster <- "HC"
stats_table_HC$cluster[] <- lapply(stats_table_HC$cluster, as.character)
class(stats_table_HC$cluster) <- "character"
strength <- rowMeans(data_HC$strength)
stats_table_HC$strength <- strength
stats_table_HC$cpl <- data_HC$cpl
stats_table_HC$transitivity <- data_HC$trans
stats_table_HC$assortativity <- data_HC$assor
stats_table_HC$clustering_coefficient <- data_HC$clustcoeff
part_coeff <- rowMeans(data_HC$part)
stats_table_HC$participation_coeff <- part_coeff
stats_table_HC$mean_con <- data_HC$meancon

for (subject in stats_table_HC[,1])
{
  c <- stats_table_HC[,1] == subject
  d <- motions[,1] == subject
  A <- as.numeric(motions[d,2])
  stats_table_HC[c,12] <- A
}

demographic_HC <- read.table(file = 'dataset_demo_info.txt', sep = "|")
demo_subject <- as.data.frame(demographic_HC[,5])
demo_subject[] <- lapply(demo_subject, as.character)
demographic_HC[,5] <- demo_subject

for (subject in stats_table_HC[,1])
{
  c <- stats_table_HC[,1] == subject
  d <- demographic_HC[,5] == subject
  A <- (as.numeric(demographic_HC[d,7])/12)
  B <- as.character(demographic_HC[d,8])
  stats_table_HC[c,2] <- A
  stats_table_HC[c,3] <- B
}

for (subject in stats_table_HC[,1])
{
  c <- stats_table_HC[,1] == subject
  d <- meds[,1] == subject
  A <- as.character(meds[d,2])
  stats_table_HC[c,13] <- A
}

meds_bin <- read.csv(file = 'meds_bin.csv', header = TRUE)
meds_subject <- as.data.frame(meds_bin[,1])
meds_subject[] <- lapply(meds_subject, as.character)
meds_bin[,1] <- meds_subject

for (subject in stats_table_HC[,1])
{
  c <- stats_table_HC[,1] == subject
  d <- meds_bin[,1] == subject
  A <- as.character(meds_bin[d,2])
  stats_table_HC[c,14] <- A
}

for (i in 1:nrow(dataset_info2))
{
  if ((dataset_info2[i,5] %in% stats_table_HC[,1]) == F)
  {
    dataset_info2[i,] <- NA
  }
}

dup_subs <- duplicated(dataset_info2$src_subject_id)
sum(dup_subs) #if sum is not zero:

# remove duplicates
dup_subs <- as.data.frame(dup_subs)

for (i in 1:nrow(dup_subs))
{
  if (dup_subs[i,] == TRUE)
  {
    dataset_info2[i,] <- NA
  }
}

# remove NAs
dataset_info2 <- dataset_info2[rowSums(is.na(dataset_info2)) != ncol(dataset_info2), ]
dataset_info2 <- as.data.frame(dataset_info2)

for (subject in stats_table_HC[,1])
{
  c <- stats_table_HC[,1] == subject
  d <- dataset_info2[,5] == subject
  A <- dataset_info2[d,3]
  stats_table_HC[c,15] <- A
}

#remove subjects with abnormally high connectivity values
for (i in 1:nrow(stats_table_HC))
{
  if (stats_table_HC[i,11] > 0.8)
  {
    stats_table_HC[i,] <- NA
  }
}

# remove NAs
stats_table_HC <- na.omit(stats_table_HC)

combined_stats <- rbind(stats_table, stats_table_HC)
combined_stats$sex <- as.factor(combined_stats$sex)
combined_stats$cluster <- as.factor(combined_stats$cluster)
combined_stats$medication <- as.factor(combined_stats$medication)
combined_stats$medication_binary <- as.factor(combined_stats$medication_binary)
combined_stats$dataset <- as.factor(combined_stats$dataset)

#--------------------------------------CAS clusters--------------------------------------
data <- readMat('CAS_RUN10_OCT_22_structure.mat')

subjects <- as.data.frame(data$subjects)

stats_table <- matrix(nrow = nrow(subjects), ncol = 15)
stats_table <- as.data.frame(stats_table)
colnames(stats_table) <- c("subject", "age","sex", "cluster", "strength", "cpl", "transitivity", "assortativity", "clustering_coefficient", "participation_coeff", "mean_con", "framewise_displacement", "medication", "medication_binary", "dataset")
stats_table$subject <- data$subjects
stats_table$subject[] <- lapply(stats_table$subject, as.character)
class(stats_table$subject) <- "character"
stats_table$cluster <- data$clusterlabels
stats_table$cluster[] <- lapply(stats_table$cluster, as.character)
class(stats_table$cluster) <- "character"
strength <- rowMeans(data$strength)
stats_table$strength <- strength
stats_table$cpl <- data$cpl
stats_table$transitivity <- data$trans
stats_table$assortativity <- data$assor
stats_table$clustering_coefficient <- data$clustcoeff
part_coeff <- rowMeans(data$part)
stats_table$participation_coeff <- part_coeff
stats_table$mean_con <- data$meancon

for (i in 1:nrow(motions))
{
  motions[i,1] <- str_remove(as.character(motions[i,1]), "sub-")
}

for (subject in stats_table[,1])
{
  c <- stats_table[,1] == subject
  d <- motions[,1] == subject
  A <- as.numeric(motions[d,2])
  stats_table[c,12] <- A
}

demographic <- read.table(file = 'cas_all_FINAL_SELECTED.txt', sep = "|")
demo_subject <- as.data.frame(demographic[,1])
demo_subject[] <- lapply(demo_subject, as.character)
demographic[,1] <- demo_subject

for (subject in stats_table[,1])
{
  c <- stats_table[,1] == subject
  d <- demographic[,1] == subject
  A <- (as.numeric(demographic[d,2])/12)
  B <- as.character(demographic[d,3])
  stats_table[c,2] <- A
  stats_table[c,3] <- B
}

for (subject in stats_table[,1])
{
  c <- stats_table[,1] == subject
  d <- meds[,1] == subject
  A <- as.character(meds[d,2])
  stats_table[c,13] <- A
}

meds_bin <- read.csv(file = 'meds_bin.csv', header = TRUE)
meds_subject <- as.data.frame(meds_bin[,1])
meds_subject[] <- lapply(meds_subject, as.character)
meds_bin[,1] <- meds_subject

for (subject in stats_table[,1])
{
  c <- stats_table[,1] == subject
  d <- meds_bin[,1] == subject
  A <- as.character(meds_bin[d,2])
  stats_table[c,14] <- A
}

#import image03_fmri file
dataset_info <- read_xlsx('image03_rsfmri.xlsx')
dataset_info2 <- dataset_info

for (i in 1:nrow(dataset_info))
{
  if ((dataset_info[i,5] %in% stats_table[,1]) == F)
  {
    dataset_info[i,] <- NA
  }
}

dup_subs <- duplicated(dataset_info$src_subject_id)
sum(dup_subs) #if sum is not zero:

# remove duplicates
dup_subs <- as.data.frame(dup_subs)

for (i in 1:nrow(dup_subs))
{
  if (dup_subs[i,] == TRUE)
  {
    dataset_info[i,] <- NA
  }
}

# remove NAs
dataset_info <- dataset_info[rowSums(is.na(dataset_info)) != ncol(dataset_info), ]
dataset_info <- as.data.frame(dataset_info)

for (subject in stats_table[,1])
{
  c <- stats_table[,1] == subject
  d <- dataset_info[,5] == subject
  A <- dataset_info[d,3]
  stats_table[c,15] <- A
}

#remove subjects with abnormally high connectivity values
for (i in 1:nrow(stats_table))
{
  if (stats_table[i,11] > 0.8)
  {
    stats_table[i,] <- NA
  }
}

# remove NAs
stats_table <- na.omit(stats_table)

combined_stats <- rbind(stats_table, stats_table_HC)
combined_stats$sex <- as.factor(combined_stats$sex)
combined_stats$cluster <- as.factor(combined_stats$cluster)
combined_stats$medication <- as.factor(combined_stats$medication)
combined_stats$medication_binary <- as.factor(combined_stats$medication_binary)
combined_stats$dataset <- as.factor(combined_stats$dataset)
#-------------------------------------Linear Models-------------------------------------
strength.model <- lm((strength - mean(strength)) ~ cluster + age + sex + framewise_displacement + medication_binary + dataset,  data = combined_stats)
cpl.model <- lm((cpl - mean(cpl)) ~ cluster + age + sex + framewise_displacement + medication_binary + dataset,  data = combined_stats)
transitivity.model <- lm((transitivity - mean(transitivity)) ~ cluster + age + sex + framewise_displacement + medication_binary + dataset,  data = combined_stats)
assortativity.model <- lm((assortativity - mean(assortativity)) ~ cluster + age + sex + framewise_displacement + medication_binary + dataset,  data = combined_stats)
clustering_coefficient.model <- lm((clustering_coefficient - mean(clustering_coefficient)) ~ cluster + age + sex + framewise_displacement + medication_binary + dataset,  data = combined_stats)
participation_coeff.model <- lm((participation_coeff - mean(participation_coeff)) ~ cluster + age + sex + framewise_displacement + medication_binary + dataset,  data = combined_stats)

#---------------------Pair-wise comparison with linear regression-------------------------
library(ggplot2)
library(tidyverse)
library(viridis)
library(ggstatsplot)
library(stats)

combined_stats$cluster <- as.factor(combined_stats$cluster)
variables <- colnames(combined_stats)
variables <- variables[5:10]

perm_comps <- as.data.frame(t(combn(as.character(unique(combined_stats$cluster)),2)))
colnames(perm_comps) <- c("group_1","group_2")

output <- list()
for(j in 1:length(variables)){
  var = variables[j]
  tempout <- perm_comps
  tempout$p <- tempout$r <- tempout$dendf <- tempout$numdf <- tempout$f <- tempout$size_g1 <- tempout$size_g2 <- tempout$median_g1 <-tempout$median_g2 <- NA
  
  for (comp in 1:nrow(perm_comps)) {
    group1 = perm_comps[comp,"group_1"]
    group1 <- factor(group1, levels=c(as.character(unique(combined_stats$cluster))))
    group2 = perm_comps[comp,"group_2"]
    group2 <- factor(group2, levels=c(as.character(unique(combined_stats$cluster))))
    
    x <- combined_stats[combined_stats$cluster == group1, ]
    y <- combined_stats[combined_stats$cluster == group2, ]
    xy <- rbind(x, y)
    
    tempout$size_g1[comp] <- length(x[,var])
    tempout$size_g2[comp] <- length(y[,var])
    tempout$median_g1[comp] <- median(x[,var],na.rm = TRUE)
    tempout$median_g2[comp] <- median(y[,var],na.rm = TRUE)
    
    test <- lm((xy[,var] - mean(xy[,var])) ~ xy$cluster + xy$age + xy$sex + xy$framewise_displacement + xy$medication_binary + xy$dataset)
    fstat <- summary(test)$fstatistic
    tempout$f[comp] <- fstat[1]
    tempout$numdf[comp] <- fstat[2]
    tempout$dendf[comp] <- fstat[3]
    tempout$r[comp] <- summary(test)$adj.r.squared
    tempout$p[comp] <- pf(fstat[1], fstat[2], fstat[3], lower.tail=FALSE)
  }
  tempout$variable <- var
  tempout$p.adjusted <- p.adjust(tempout$p,method = "BH")
  output[[var]] <- tempout
}

output <- bind_rows(output)

#CAS figure
#PLOT
sign_vars<- unique(output$variable)

plot_output <- list()
for (var in sign_vars)
{
  limit <- max(combined_stats[,var])
  pdata <- as.data.frame(matrix(nrow = nrow(combined_stats), ncol = 2))
  colnames(pdata) <- c("variable", "cluster")
  pdata$variable <- combined_stats[,var]
  pdata$cluster <- combined_stats$cluster
  if (var == "clustering_coefficient") {varname <- "cc"}
  else if (var == "participation_coeff") {varname <- "pc"}
  else {varname <- var}
  temp_plot <- ggplot(pdata, aes(x = cluster, y = variable, fill=cluster)) + 
    geom_boxplot(width=0.4) +
    theme_classic() +
    labs(y= varname, x = "CAS cluster") +
    theme(legend.position="none", legend.title=element_blank()) +
    theme(axis.text=element_text(size=10))
  for (comp in 1:nrow(perm_comps)) {
    group1 = perm_comps[comp,"group_1"]
    group2 = perm_comps[comp,"group_2"]
    new <- output[grep(var, output$variable), ]
    new2 <- new[new$group_1 %like% as.character(group1), ]
    new2 <- new2[new2$group_2 %like% as.character(group2), ]
    if (new2$p.adjusted < .05) { 
      temp_plot <- temp_plot + geom_signif(comparisons=list(c(as.character(group1), as.character(group2))), annotations="*", y_position = limit, tip_length = 0.05, vjust=0.4)
      limit <- limit+(limit/15)
    }
  }
  plot_output[[var]] <- temp_plot
  rm(temp_plot)
}

ggarrange(plot_output[[1]], plot_output[[2]], plot_output[[3]], plot_output[[4]], plot_output[[5]], plot_output[[6]], labels = c("A", "B", "C", "D", "E", "F"), ncol=3, nrow=2, common.legend = TRUE, legend="right") +
  ggtitle("Global network properties per each anxiety cluster and controls \n") + 
  theme(plot.title = element_text(face="bold", hjust = 0.5, size = 14)) 

ggsave(filename = "CAS_global_figure.tiff", width = 10, height = 6, device='tiff', dpi=700)
dev.off()

#YMRS figure
#PLOT
sign_vars<- unique(output$variable)

plot_output <- list()
for (var in sign_vars)
{
  limit <- max(combined_stats[,var])
  pdata <- as.data.frame(matrix(nrow = nrow(combined_stats), ncol = 2))
  colnames(pdata) <- c("variable", "cluster")
  pdata$variable <- combined_stats[,var]
  pdata$cluster <- combined_stats$cluster
  if (var == "clustering_coefficient") {varname <- "cc"}
  else if (var == "participation_coeff") {varname <- "pc"}
  else {varname <- var}
  temp_plot <- ggplot(pdata, aes(x = cluster, y = variable, fill=cluster)) + 
    geom_boxplot(width=0.4) +
    theme_classic() +
    labs(y= varname, x = "YMRS cluster") +
    theme(legend.position="none", legend.title=element_blank()) +
    theme(axis.text=element_text(size=10))
  for (comp in 1:nrow(perm_comps)) {
    group1 = perm_comps[comp,"group_1"]
    group2 = perm_comps[comp,"group_2"]
    new <- output[grep(var, output$variable), ]
    new2 <- new[new$group_1 %like% as.character(group1), ]
    new2 <- new2[new2$group_2 %like% as.character(group2), ]
    if (new2$p.adjusted < .05) { 
      temp_plot <- temp_plot + geom_signif(comparisons=list(c(as.character(group1), as.character(group2))), annotations="*", y_position = limit, tip_length = 0.05, vjust=0.4)
      limit <- limit+(limit/15)
    }
  }
  plot_output[[var]] <- temp_plot
  rm(temp_plot)
}

ggarrange(plot_output[[1]], plot_output[[2]], plot_output[[3]], plot_output[[4]], plot_output[[5]], plot_output[[6]], labels = c("A", "B", "C", "D", "E", "F"), ncol=3, nrow=2, common.legend = TRUE, legend="right") +
  ggtitle("Global network properties per each mania cluster and controls \n") + 
  theme(plot.title = element_text(face="bold", hjust = 0.5, size = 14)) 

ggsave(filename = "YMRS_global_figure.tiff", width = 10, height = 6, device='tiff', dpi=700)
dev.off()

#get means and SD
table(combined_stats$cluster)

means1 <- matrix(nrow=71, ncol=2)
means1 <- as.data.frame(means1)

n = 1
for (i in 1:nrow(combined_stats))
{
  if (combined_stats[i,4] == 0)
  {
    means1[n,1] <- combined_stats[i,4]
    means1[n,2] <- combined_stats[i,8]
    n = n+1
  }
}

means2 <- matrix(nrow=33, ncol=2)
means2 <- as.data.frame(means2)

n = 1
for (i in 1:nrow(combined_stats))
{
  if (combined_stats[i,4] == 1)
  {
    means2[n,1] <- combined_stats[i,4]
    means2[n,2] <- combined_stats[i,8]
    n = n+1
  }
}

means3 <- matrix(nrow=70, ncol=2)
means3 <- as.data.frame(means1)

n = 1
for (i in 1:nrow(combined_stats))
{
  if (combined_stats[i,4] == 'HC')
  {
    means3[n,1] <- combined_stats[i,4]
    means3[n,2] <- combined_stats[i,8]
    n = n+1
  }
}

