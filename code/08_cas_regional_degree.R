library(DescTools)
library(R.matlab)
library(dgof)
library(stringr)
library(ggplot2)
library(ggpubr)
library(readxl)
library(ggsignif)
library(ggpubr)
library(tibble)

#----------------------------------Network measures---------------------------------------
data2 <- readMat('HC_RUN10_OCT_22_structure.mat')

subjects <- as.data.frame(data2$subjects)

motions <- read.csv(file = 'motion_parameters.txt')
motions_subject <- as.data.frame(motions[,1])
motions_subject[] <- lapply(motions_subject, as.character)
motions[,1] <- motions_subject

HC_stats_table <- matrix(nrow = nrow(subjects), ncol = 15)
HC_stats_table <- as.data.frame(HC_stats_table)
colnames(HC_stats_table) <- c("subject", "age","sex", "cluster", "strength", "cpl", "transitivity", "assortativity", "clustering_coefficient", "participation_coeff", "mean_con", "framewise_displacement", "medication", "medication_binary", "dataset")
HC_stats_table$subject <- data2$subjects
HC_stats_table$subject[] <- lapply(HC_stats_table$subject, as.character)
class(HC_stats_table$subject) <- "character"
HC_stats_table$cluster <- 'HC'
HC_stats_table$cluster[] <- lapply(HC_stats_table$cluster, as.character)
class(HC_stats_table$cluster) <- "character"
strength <- rowMeans(data2$strength)
HC_stats_table$strength <- strength
HC_stats_table$cpl <- data2$cpl
HC_stats_table$transitivity <- data2$trans
HC_stats_table$assortativity <- data2$assor
HC_stats_table$clustering_coefficient <- data2$clustcoeff
part_coeff <- rowMeans(data2$part)
HC_stats_table$participation_coeff <- part_coeff
HC_stats_table$mean_con <- data2$meancon

for (i in 1:nrow(motions))
{
  motions[i,1] <- str_remove(as.character(motions[i,1]), "sub-")
}

for (subject in HC_stats_table[,1])
{
  c <- HC_stats_table[,1] == subject
  d <- motions[,1] == subject
  A <- as.numeric(motions[d,2])
  HC_stats_table[c,12] <- A
}

demographic_HC <- read.table(file = 'dataset_demo_info.txt', sep = "|")
demo_subject <- as.data.frame(demographic_HC[,5])
demo_subject[] <- lapply(demo_subject, as.character)
demographic_HC[,5] <- demo_subject

for (subject in HC_stats_table[,1])
{
  c <- HC_stats_table[,1] == subject
  d <- demographic_HC[,5] == subject
  A <- (as.numeric(demographic_HC[d,7])/12)
  B <- as.character(demographic_HC[d,8])
  HC_stats_table[c,2] <- A
  HC_stats_table[c,3] <- B
}

meds <- read.csv(file = 'meds_BD_HC.csv', header = TRUE)
meds <- meds[,2:3]
meds_subject <- as.data.frame(meds[,1])
meds_subject[] <- lapply(meds_subject, as.character)
meds[,1] <- meds_subject

for (subject in HC_stats_table[,1])
{
  c <- HC_stats_table[,1] == subject
  d <- meds[,1] == subject
  A <- as.character(meds[d,2])
  HC_stats_table[c,13] <- A
}

meds_bin <- read.csv(file = 'meds_bin.csv', header = TRUE)
meds_subject <- as.data.frame(meds_bin[,1])
meds_subject[] <- lapply(meds_subject, as.character)
meds_bin[,1] <- meds_subject

for (subject in HC_stats_table[,1])
{
  c <- HC_stats_table[,1] == subject
  d <- meds_bin[,1] == subject
  A <- as.character(meds_bin[d,2])
  HC_stats_table[c,14] <- A
}

#import image03_fmri file
dataset_info <- read_xlsx('image03_rsfmri.xlsx')
dataset_info2 <- dataset_info

for (i in 1:nrow(dataset_info))
{
  if ((dataset_info[i,5] %in% HC_stats_table[,1]) == F)
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

for (subject in HC_stats_table[,1])
{
  c <- HC_stats_table[,1] == subject
  d <- dataset_info[,5] == subject
  A <- dataset_info[d,3]
  HC_stats_table[c,15] <- A
}

#remove subjects with abnormally high connectivity values
for (i in 1:nrow(HC_stats_table))
{
  if (HC_stats_table[i,11] > 0.8)
  {
    HC_stats_table[i,] <- NA
  }
}

# remove NAs
HC_stats_table <- na.omit(HC_stats_table)

#-------------------------------------------------------------------------
degree_regional <- matrix(nrow = nrow(subjects), ncol = 367)
degree_regional <- as.data.frame(degree_regional)
degree_regional[,8:367] <- data2$deg
names(degree_regional)[names(degree_regional) == 'V1'] <- 'cluster'
degree_regional$cluster <- 'HC'
degree_regional$cluster[] <- lapply(degree_regional$cluster, as.character)
class(degree_regional$cluster) <- "character"
names(degree_regional)[names(degree_regional) == 'V2'] <- 'age'
degree_regional$age <- HC_stats_table$age
names(degree_regional)[names(degree_regional) == 'V3'] <- 'sex'
degree_regional$sex <- HC_stats_table$sex
names(degree_regional)[names(degree_regional) == 'V4'] <- 'framewise_displacement'
degree_regional$framewise_displacement <- HC_stats_table$framewise_displacement
names(degree_regional)[names(degree_regional) == 'V5'] <- 'medication'
degree_regional$medication <- HC_stats_table$medication
names(degree_regional)[names(degree_regional) == 'V6'] <- 'medication_binary'
degree_regional$medication_binary <- HC_stats_table$medication_binary
names(degree_regional)[names(degree_regional) == 'V7'] <- 'dataset'
degree_regional$dataset <- HC_stats_table$dataset
degree_regional2 <- degree_regional

#--------------------------add CAS cluster data----------------------------
data <- readMat('CAS_RUN10_OCT_22_structure.mat')

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

#-----------------------regional degree----------------------------
degree_regional <- matrix(nrow = nrow(subjects), ncol = 367)
degree_regional <- as.data.frame(degree_regional)
degree_regional[,8:367] <- data$deg
names(degree_regional)[names(degree_regional) == 'V1'] <- 'cluster'
degree_regional$cluster <- data$clusterlabels
degree_regional$cluster[] <- lapply(degree_regional$cluster, as.character)
class(degree_regional$cluster) <- "character"
names(degree_regional)[names(degree_regional) == 'V2'] <- 'age'
degree_regional$age <- stats_table$age
names(degree_regional)[names(degree_regional) == 'V3'] <- 'sex'
degree_regional$sex <- stats_table$sex
names(degree_regional)[names(degree_regional) == 'V4'] <- 'framewise_displacement'
degree_regional$framewise_displacement <- stats_table$framewise_displacement
names(degree_regional)[names(degree_regional) == 'V5'] <- 'medication'
degree_regional$medication <- stats_table$medication
names(degree_regional)[names(degree_regional) == 'V6'] <- 'medication_binary'
degree_regional$medication_binary <- stats_table$medication_binary
names(degree_regional)[names(degree_regional) == 'V7'] <- 'dataset'
degree_regional$dataset <- stats_table$dataset

degree_regional <- rbind(degree_regional, degree_regional2) #add healthy controls

#-----------stats-----------
my_lms <- lapply(8:367, function(i) lm((degree_regional[,i] - mean(degree_regional[,i])) ~ degree_regional$cluster + degree_regional$age + degree_regional$sex + degree_regional$framewise_displacement + degree_regional$medication_binary + degree_regional$dataset))

#get regions labels
regionlabels <- read.table('HCP-MMP1_on_MNI152_ICBM2009a_nlin.txt', sep = "")
regionlabels <- as.data.frame(regionlabels[,2])
regionlabels[] <- lapply(regionlabels, as.character)
regionlabels <- as.data.frame(regionlabels)
colnames(regionlabels) <- c("REG")
extra_labels <- regionlabels$REG
regionlabels <- gsub('^.', '', regionlabels$REG)
regionlabels <- gsub('^.', '', regionlabels)
regionlabels <- gsub('_ROI','', regionlabels)
extra_labels <- gsub('_ROI','', extra_labels)
regionlabels <- gsub('ROI','', regionlabels)
extra_labels <- gsub('ROI','', extra_labels)
r_extra_labels <- gsub('^.','R', extra_labels)
regionlabels <- as.data.frame(regionlabels)
extra_labels <- as.data.frame(extra_labels)
r_extra_labels <- as.data.frame(r_extra_labels)
all_regions <- matrix(nrow = 360, ncol = 1)
all_regions <- as.data.frame(all_regions)
regionlabels[] <- lapply(regionlabels, as.character)
extra_labels[] <- lapply(extra_labels, as.character)
r_extra_labels[] <- lapply(r_extra_labels, as.character)
all_regions[] <- lapply(all_regions, as.character)
all_regions[1:180,] <- regionlabels
#regionlabels <- gsub('L_', 'R_', regionlabels$REG)
all_regions[181:360,] <- regionlabels
names(degree_regional) <- c(colnames(degree_regional[,1:7]),extra_labels$extra_labels,r_extra_labels$r_extra_labels)

rstat <- matrix(nrow = 360, ncol = 6)
rstat <- as.data.frame(rstat)
rstat[,1] <- all_regions
rstat[1:180, 2] <- 'left'
rstat[181:360, 2] <- 'right'
for (i in 1:360)
{
  sum1 <- summary(my_lms[[i]])
  lm_p <- pf(sum1$fstatistic[1], sum1$fstatistic[2], sum1$fstatistic[3], lower.tail = FALSE)
  rstat[i,3] <- lm_p
  rstat[i,4] <- sum1$adj.r.squared
  rstat[i,5] <- sum1$fstatistic[1]
}
rstat_adj <- rstat
rstat_adj[,3] <- p.adjust(rstat[,3], method = 'fdr')
rstat_adj <- as.data.frame(rstat_adj)

names(rstat_adj) <- c('region', 'hemi', 'pval','rstat', 'fstat', 'sign')

#----------------------Brain plot--------------------------
options(repos = c(
  ggseg = 'https://ggseg.r-universe.dev',
  CRAN = 'https://cloud.r-project.org'))

library(ggsegGlasser)
library(ggseg)
library(ggplot2)
library(scico)

glass_object <- ggsegGlasser::glasser$data

glass_object$label[170] <- "lh_L_10pp"
glass_object$region[170] <- "10pp"
glass_object$side[40] <- "medial"

glasser$data <- glass_object

rstat_adj$sign <- ifelse(rstat_adj$pval <0.05,"Yes","No")

#--------------------------------Figures-----------------------------------
Fig1 <- rstat_adj %>% ggseg(mapping=aes(fill=(fstat),color=sign),  atlas = "glasser") +
  scale_fill_scico(palette = "bilbao",direction = 1) +
  scale_color_manual(values = c("white","red"), na.value = "white") +
  theme(legend.position = "bottom", axis.title.x = element_blank(), 
        axis.title.y = element_blank(), axis.text.y = element_blank(),
        strip.text.x = element_text(size = 14)) +
  ggtitle("Thresholded at 10% \n")
rstat_adj_10 <- rstat_adj

Fig2 <- rstat_adj %>% ggseg(mapping=aes(fill=(fstat),color=sign),  atlas = "glasser") +
  scale_fill_scico(palette = "bilbao",direction = 1) +
  scale_color_manual(values = c("white","red"), na.value = "white") +
  theme(legend.position = "bottom", axis.title.x = element_blank(), 
        axis.title.y = element_blank(), axis.text.y = element_blank(),
        strip.text.x = element_text(size = 14)) +
  ggtitle("Thresholded at 2% \n")
rstat_adj_2 <- rstat_adj

Fig3 <- rstat_adj %>% ggseg(mapping=aes(fill=(fstat),color=sign),  atlas = "glasser") +
  scale_fill_scico(palette = "bilbao",direction = 1) +
  scale_color_manual(values = c("white","red"), na.value = "white") +
  theme(legend.position = "bottom", axis.title.x = element_blank(), 
        axis.title.y = element_blank(), axis.text.y = element_blank(),
        strip.text.x = element_text(size = 14)) +
  ggtitle("Thresholded at 5% \n")
rstat_adj_5 <- rstat_adj

Fig4 <- rstat_adj %>% ggseg(mapping=aes(fill=(fstat),color=sign),  atlas = "glasser") +
  scale_fill_scico(palette = "bilbao",direction = 1) +
  scale_color_manual(values = c("white","red"), na.value = "white") +
  theme(legend.position = "bottom", axis.title.x = element_blank(), 
        axis.title.y = element_blank(), axis.text.y = element_blank(),
        strip.text.x = element_text(size = 14)) +
  ggtitle("Thresholded at 20% \n")
rstat_adj_20 <- rstat_adj

ggarrange(Fig1, Fig2, Fig3, Fig4, labels = c("A", "B", "C", "D"), ncol=2, nrow=2, common.legend = TRUE, legend="right") +
  ggtitle("Regional differences in degree between anxiety subgroups for different thresholds \n") + 
  theme(plot.title = element_text(face="bold", hjust = 0.5, size = 16))

rstat_adj_10[is.na(rstat_adj_10)] <- 0
rstat_adj_2[is.na(rstat_adj_2)] <- 0
rstat_adj_5[is.na(rstat_adj_5)] <- 0
rstat_adj_20[is.na(rstat_adj_20)] <- 0

new_data=data.frame(col10=rstat_adj_10$rstat,col2=rstat_adj_2$rstat,
                    col5=rstat_adj_5$rstat,col20=rstat_adj_20$rstat)
cor(new_data, method = c("spearman"))

