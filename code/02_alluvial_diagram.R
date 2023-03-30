library(ggplot2)
#install.packages('ggalluvial')
library(ggalluvial)

setwd("/Users/janahermans/Documents")
df <- read.table('newcom_label_ymrs_cas_FINAL.txt', sep = "|", header = T, row.names = 1)
df <- df[,-c(1,3,4)]

figure1 <- ggplot(as.data.frame(df),
       aes(axis1 = CAS_cluster, axis2 = YMRS_cluster)) +
  geom_alluvium(aes(fill = as.character(YMRS_cluster)), width = 1/12) +
  geom_stratum(width = 1/12, fill = "white", color = "grey") +
  geom_label(stat = "stratum", aes(label = after_stat(stratum)), min.y = 200) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), color = "grey35") +
  scale_x_discrete(limits = c("Anxiety cluster", "Mania cluster")) +
  theme_classic() +
  theme_void() +
  theme(legend.position = "none", axis.ticks.y=element_blank()) +
  theme(axis.text.x = element_text(face = "plain"), plot.margin = margin(0, 0, 10, 0))

figure1 + ggtitle('Correspondence between mania and anxiety cluster labels') +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
