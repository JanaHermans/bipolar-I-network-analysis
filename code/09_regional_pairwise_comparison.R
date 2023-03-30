library(ggplot2)
library(tidyverse)
library(viridis)
library(ggstatsplot)
library(stats)

#--------------------------linear model STRENGTH-----------------------------------------
strength_regional$cluster <- as.factor(strength_regional$cluster)
rstat_adj_pair <- rstat_adj[rstat_adj$sign == "Yes",]
rstat_adj_pair$region <- paste0(rstat_adj_pair$hemi,"_",rstat_adj_pair$region)
rstat_adj_pair$region <- gsub("right_", "R_", rstat_adj_pair$region)
rstat_adj_pair$region <- gsub("left_", "L_", rstat_adj_pair$region)

perm_comps <- as.data.frame(t(combn(as.character(unique(strength_regional$cluster)),2)))
colnames(perm_comps) <- c("group_1","group_2")

output <- list()
for(region in 1:length(unique(rstat_adj_pair$region))){
  i = rstat_adj_pair$region[region]
  tempout <- perm_comps
  tempout$p <- tempout$r <- tempout$dendf <- tempout$numdf <- tempout$f <- tempout$size_g1 <- tempout$size_g2 <- tempout$median_g1 <-tempout$median_g2 <- NA
  
  for (comp in 1:nrow(perm_comps)) {
    group1 = perm_comps[comp,"group_1"]
    group1 <- factor(group1, levels=c(as.character(unique(strength_regional$cluster))))
    group2 = perm_comps[comp,"group_2"]
    group2 <- factor(group2, levels=c(as.character(unique(strength_regional$cluster))))
    
    x <- strength_regional[strength_regional$cluster == group1, ]
    y <- strength_regional[strength_regional$cluster == group2, ]
    xy <- rbind(x, y)
    
    tempout$size_g1[comp] <- length(x[,i])
    tempout$size_g2[comp] <- length(y[,i])
    tempout$median_g1[comp] <- median(x[,i],na.rm = TRUE)
    tempout$median_g2[comp] <- median(y[,i],na.rm = TRUE)
    
    test <- lm((xy[,i] - mean(xy[,i])) ~ xy$cluster + xy$age + xy$sex + xy$framewise_displacement + xy$medication_binary + xy$dataset)
    fstat <- summary(test)$fstatistic
    tempout$f[comp] <- fstat[1]
    tempout$numdf[comp] <- fstat[2]
    tempout$dendf[comp] <- fstat[3]
    tempout$r[comp] <- summary(test)$adj.r.squared
    tempout$p[comp] <- pf(fstat[1], fstat[2], fstat[3], lower.tail=FALSE)
  }
  tempout$region <- i
  output[[i]] <- tempout
}

output <- bind_rows(output)
output$p_adjusted <- p.adjust(output$p,method = "BH")

#-------------------------plotting--------------------------------
sign_regions <- unique(output$region)

plot_output <- list()
for (i in sign_regions)
{
  if (grepl("R_", i))
  {
    region_hemi <- "right "
    region_name <- gsub("R_", "", i)
  }
  else if (grepl("L_", i))
  {
    region_hemi <- "left "
    region_name <- gsub("L_", "", i)
  }
  pdata <- as.data.frame(matrix(nrow = nrow(strength_regional), ncol = 2))
  colnames(pdata) <- c("region", "cluster")
  pdata$region <- strength_regional[,i]
  pdata$cluster <- strength_regional$cluster
  temp_plot <- ggplot(pdata, aes(x = cluster, y = region, fill=cluster)) + 
    geom_boxplot(width=0.4) +
    theme_classic() +
    theme(axis.title.x=element_blank(), axis.title.y = element_blank()) +
    ggtitle(paste0(region_hemi, region_name)) +
    theme(legend.position="none", legend.title=element_blank()) +
    theme(axis.text=element_text(size=10))
  for (comp in 1:nrow(perm_comps)) {
    group1 = perm_comps[comp,"group_1"]
    group2 = perm_comps[comp,"group_2"]
    new <- output[grep(i, output$region), ]
    new2 <- new[new$group_1 %like% as.character(group1), ]
    new2 <- new2[new2$group_2 %like% as.character(group2), ]
    if (new2$p_adjusted < .05) { temp_plot <- temp_plot + geom_signif(comparisons=list(c(as.character(group1), as.character(group2))), annotations="*", tip_length = 0.05, vjust=0.4)}
  }
  plot_output[[i]] <- temp_plot
  rm(temp_plot)
}

#CAS
figure <- ggarrange(plot_output[[1]],plot_output[[2]],plot_output[[3]],plot_output[[4]],plot_output[[5]],plot_output[[6]],plot_output[[7]],plot_output[[8]],plot_output[[9]], plot_output[[10]],
                    plot_output[[11]],plot_output[[12]],plot_output[[13]],plot_output[[14]],plot_output[[15]],
                    plot_output[[16]],plot_output[[17]],plot_output[[18]],plot_output[[19]],plot_output[[20]],
                    plot_output[[21]],plot_output[[22]],plot_output[[23]],plot_output[[24]],plot_output[[25]], 
                    labels = c("A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T","U","V","W","X","Y"), 
                    ncol=5, nrow=5, common.legend = TRUE, legend="right") +
  ggtitle("Pairwise comparison of regional strength for anxiety clusters\n") + 
  theme(plot.title = element_text(face="bold", hjust = 0.5, size = 14)) 

annotate_figure(figure,
                bottom = text_grob("CAS cluster", face = "bold", size = 12),
                left = text_grob("Strength", face = "bold", rot = 90, size = 12),
)

ggsave(filename = "CAS_suppl_regional_strength.tiff", width = 12, height = 12, device='tiff', dpi=700)
dev.off()

#YMRS
figure <- ggarrange(plot_output[[1]],plot_output[[2]],plot_output[[3]],plot_output[[4]],plot_output[[5]],
                    plot_output[[6]],plot_output[[7]],plot_output[[8]],plot_output[[9]], plot_output[[10]],
                    plot_output[[11]],plot_output[[12]],plot_output[[13]],plot_output[[14]],plot_output[[15]],
                    plot_output[[16]],plot_output[[17]],plot_output[[18]],plot_output[[19]],plot_output[[20]],
                    plot_output[[21]],plot_output[[22]],plot_output[[23]],plot_output[[24]],plot_output[[25]],plot_output[[26]],
                    labels = c("A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T","U","V","W","X","Y","Z"), 
                    ncol=5, nrow=6, common.legend = TRUE, legend="right") +
  ggtitle("Pairwise comparison of regional strength for mania clusters\n") + 
  theme(plot.title = element_text(face="bold", hjust = 0.5, size = 14)) 

annotate_figure(figure,
                bottom = text_grob("YMRS cluster", face = "bold", size = 12),
                left = text_grob("Strength", face = "bold", rot = 90, size = 12),
)

ggsave(filename = "YMRS_suppl_regional_strength.tiff", width = 12, height = 15, device='tiff', dpi=700)
dev.off()

#--------------------------linear model DEGREE-----------------------------------------
degree_regional$cluster <- as.factor(degree_regional$cluster)
rstat_adj_pair <- rstat_adj[rstat_adj$sign == "Yes",]
rstat_adj_pair$region <- paste0(rstat_adj_pair$hemi,"_",rstat_adj_pair$region)
rstat_adj_pair$region <- gsub("right_", "R_", rstat_adj_pair$region)
rstat_adj_pair$region <- gsub("left_", "L_", rstat_adj_pair$region)

perm_comps <- as.data.frame(t(combn(as.character(unique(degree_regional$cluster)),2)))
colnames(perm_comps) <- c("group_1","group_2")

output <- list()
for(region in 1:length(unique(rstat_adj_pair$region))){
  i = rstat_adj_pair$region[region]
  tempout <- perm_comps
  tempout$p <- tempout$r <- tempout$dendf <- tempout$numdf <- tempout$f <- tempout$size_g1 <- tempout$size_g2 <- tempout$median_g1 <-tempout$median_g2 <- NA
  
  for (comp in 1:nrow(perm_comps)) {
    group1 = perm_comps[comp,"group_1"]
    group1 <- factor(group1, levels=c(as.character(unique(degree_regional$cluster))))
    group2 = perm_comps[comp,"group_2"]
    group2 <- factor(group2, levels=c(as.character(unique(degree_regional$cluster))))
    
    x <- degree_regional[degree_regional$cluster == group1, ]
    y <- degree_regional[degree_regional$cluster == group2, ]
    xy <- rbind(x, y)
    
    tempout$size_g1[comp] <- length(x[,i])
    tempout$size_g2[comp] <- length(y[,i])
    tempout$median_g1[comp] <- median(x[,i],na.rm = TRUE)
    tempout$median_g2[comp] <- median(y[,i],na.rm = TRUE)
    
    test <- lm((xy[,i] - mean(xy[,i])) ~ xy$cluster + xy$age + xy$sex + xy$framewise_displacement + xy$medication_binary + xy$dataset)
    fstat <- summary(test)$fstatistic
    tempout$f[comp] <- fstat[1]
    tempout$numdf[comp] <- fstat[2]
    tempout$dendf[comp] <- fstat[3]
    tempout$r[comp] <- summary(test)$adj.r.squared
    tempout$p[comp] <- pf(fstat[1], fstat[2], fstat[3], lower.tail=FALSE)
  }
  tempout$region <- i
  output[[i]] <- tempout
}

output <- bind_rows(output)
output$p_adjusted <- p.adjust(output$p,method = "BH")

#-------------------------plotting--------------------------------
sign_regions <- unique(output$region)

plot_output <- list()
for (i in sign_regions)
{
  if (grepl("R_", i))
  {
    region_hemi <- "right "
    region_name <- gsub("R_", "", i)
  }
  else if (grepl("L_", i))
  {
    region_hemi <- "left "
    region_name <- gsub("L_", "", i)
  }
  pdata <- as.data.frame(matrix(nrow = nrow(degree_regional), ncol = 2))
  colnames(pdata) <- c("region", "cluster")
  pdata$region <- degree_regional[,i]
  pdata$cluster <- degree_regional$cluster
  temp_plot <- ggplot(pdata, aes(x = cluster, y = region, fill=cluster)) + 
    geom_boxplot(width=0.4) +
    theme_classic() +
    theme(axis.title.x=element_blank(), axis.title.y = element_blank()) +
    ggtitle(paste0(region_hemi, region_name)) +
    theme(legend.position="none", legend.title=element_blank()) +
    theme(axis.text=element_text(size=10))
  for (comp in 1:nrow(perm_comps)) {
    group1 = perm_comps[comp,"group_1"]
    group2 = perm_comps[comp,"group_2"]
    new <- output[grep(i, output$region), ]
    new2 <- new[new$group_1 %like% as.character(group1), ]
    new2 <- new2[new2$group_2 %like% as.character(group2), ]
    if (new2$p_adjusted < .05) { temp_plot <- temp_plot + geom_signif(comparisons=list(c(as.character(group1), as.character(group2))), annotations="*", tip_length = 0.05, vjust=0.4)}
  }
  plot_output[[i]] <- temp_plot
  rm(temp_plot)
}

#CAS
figure <- ggarrange(plot_output[[1]], plot_output[[2]], plot_output[[3]], plot_output[[4]], plot_output[[5]], plot_output[[6]], plot_output[[7]], plot_output[[8]], plot_output[[9]], plot_output[[10]], plot_output[[11]], plot_output[[12]], plot_output[[13]], plot_output[[14]], plot_output[[15]], plot_output[[16]], labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P"), ncol=4, nrow=4, common.legend = TRUE, legend="right") +
  ggtitle("Pairwise comparison of regional degree for anxiety clusters\n") + 
  theme(plot.title = element_text(face="bold", hjust = 0.5, size = 14)) 

annotate_figure(figure,
                bottom = text_grob("CAS cluster", face = "bold", size = 12),
                left = text_grob("Degree", face = "bold", rot = 90, size = 12),
)

ggsave(filename = "CAS_suppl_regional_degree.tiff", width = 10, height = 10, device='tiff', dpi=700)
dev.off()

#YMRS
figure <- ggarrange(plot_output[[1]], plot_output[[2]], plot_output[[3]], plot_output[[4]], plot_output[[5]], plot_output[[6]], plot_output[[7]], plot_output[[8]], plot_output[[9]], plot_output[[10]], plot_output[[11]], plot_output[[12]], plot_output[[13]], plot_output[[14]], labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N"), ncol=4, nrow=4, common.legend = TRUE, legend="right") +
  ggtitle("Pairwise comparison of regional degree for mania clusters\n") + 
  theme(plot.title = element_text(face="bold", hjust = 0.5, size = 14)) 

annotate_figure(figure,
                bottom = text_grob("YMRS cluster", face = "bold", size = 12),
                left = text_grob("Degree", face = "bold", rot = 90, size = 12),
)

ggsave(filename = "YMRS_suppl_regional_degree.tiff", width = 10, height = 10, device='tiff', dpi=700)
dev.off()

#--------------------------linear model-----------------------------------------
participation_regional$cluster <- as.factor(participation_regional$cluster)
rstat_adj_pair <- rstat_adj[rstat_adj$sign == "Yes",]
rstat_adj_pair$region <- paste0(rstat_adj_pair$hemi,"_",rstat_adj_pair$region)
rstat_adj_pair$region <- gsub("right_", "R_", rstat_adj_pair$region)
rstat_adj_pair$region <- gsub("left_", "L_", rstat_adj_pair$region)

perm_comps <- as.data.frame(t(combn(as.character(unique(participation_regional$cluster)),2)))
colnames(perm_comps) <- c("group_1","group_2")

output <- list()
for(region in 1:length(unique(rstat_adj_pair$region))){
  i = rstat_adj_pair$region[region]
  tempout <- perm_comps
  tempout$p <- tempout$r <- tempout$dendf <- tempout$numdf <- tempout$f <- tempout$size_g1 <- tempout$size_g2 <- tempout$median_g1 <-tempout$median_g2 <- NA
  
  for (comp in 1:nrow(perm_comps)) {
    group1 = perm_comps[comp,"group_1"]
    group1 <- factor(group1, levels=c(as.character(unique(participation_regional$cluster))))
    group2 = perm_comps[comp,"group_2"]
    group2 <- factor(group2, levels=c(as.character(unique(participation_regional$cluster))))
    
    x <- participation_regional[participation_regional$cluster == group1, ]
    y <- participation_regional[participation_regional$cluster == group2, ]
    xy <- rbind(x, y)
    
    tempout$size_g1[comp] <- length(x[,i])
    tempout$size_g2[comp] <- length(y[,i])
    tempout$median_g1[comp] <- median(x[,i],na.rm = TRUE)
    tempout$median_g2[comp] <- median(y[,i],na.rm = TRUE)
    
    test <- lm((xy[,i] - mean(xy[,i])) ~ xy$cluster + xy$age + xy$sex + xy$framewise_displacement + xy$medication_binary + xy$dataset)
    fstat <- summary(test)$fstatistic
    tempout$f[comp] <- fstat[1]
    tempout$numdf[comp] <- fstat[2]
    tempout$dendf[comp] <- fstat[3]
    tempout$r[comp] <- summary(test)$adj.r.squared
    tempout$p[comp] <- pf(fstat[1], fstat[2], fstat[3], lower.tail=FALSE)
  }
  tempout$region <- i
  output[[i]] <- tempout
}

output <- bind_rows(output)
output$p_adjusted <- p.adjust(output$p,method = "BH")

#-------------------------plotting--------------------------------
sign_regions <- unique(output$region)

plot_output <- list()
for (i in sign_regions)
{
  if (grepl("R_", i))
  {
    region_hemi <- "right "
    region_name <- gsub("R_", "", i)
  }
  else if (grepl("L_", i))
  {
    region_hemi <- "left "
    region_name <- gsub("L_", "", i)
  }
  pdata <- as.data.frame(matrix(nrow = nrow(participation_regional), ncol = 2))
  colnames(pdata) <- c("region", "cluster")
  pdata$region <- participation_regional[,i]
  pdata$cluster <- participation_regional$cluster
  temp_plot <- ggplot(pdata, aes(x = cluster, y = region, fill=cluster)) + 
    geom_boxplot(width=0.4) +
    theme_classic() +
    theme(axis.title.x=element_blank(), axis.title.y = element_blank()) +
    ggtitle(paste0(region_hemi, region_name)) +
    theme(legend.position="none", legend.title=element_blank()) +
    theme(axis.text=element_text(size=10))
  for (comp in 1:nrow(perm_comps)) {
    group1 = perm_comps[comp,"group_1"]
    group2 = perm_comps[comp,"group_2"]
    new <- output[grep(i, output$region), ]
    new2 <- new[new$group_1 %like% as.character(group1), ]
    new2 <- new2[new2$group_2 %like% as.character(group2), ]
    if (new2$p_adjusted < .05) { temp_plot <- temp_plot + geom_signif(comparisons=list(c(as.character(group1), as.character(group2))), annotations="*", tip_length = 0.05, vjust=0.4)}
  }
  plot_output[[i]] <- temp_plot
  rm(temp_plot)
}

#CAS
figure <- ggarrange(plot_output[[1]], plot_output[[2]], plot_output[[3]], plot_output[[4]], plot_output[[5]], labels = c("A", "B", "C", "D", "E"), ncol=3, nrow=2, common.legend = TRUE, legend="right") +
  ggtitle("Pairwise comparison of regional participation coefficients for anxiety clusters\n") + 
  theme(plot.title = element_text(face="bold", hjust = 0.5, size = 14)) 

annotate_figure(figure,
                bottom = text_grob("CAS cluster", face = "bold", size = 12),
                left = text_grob("Particpation coefficient", face = "bold", rot = 90, size = 12),
)

ggsave(filename = "CAS_suppl_regional_pc.tiff", width = 10, height = 7, device='tiff', dpi=700)
dev.off()

