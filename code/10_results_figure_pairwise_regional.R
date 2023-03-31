# This script creates a results table for regional pairwise comparisons
# for strength, degree, and participation coefficient (added in Supplement).

table_strength <- rstat_adj
table_strength$sign <- NULL
table_strength$hemi <- NULL
table_strength <- table_strength[, c(1, 3, 4, 2)]
colnames(table_strength)[which(names(table_strength) == "pval")] <- "pval_adjusted"
print(sum1$fstatistic[2]) #add to table manually
print(sum1$fstatistic[3]) #add to table manually

table_degree <- rstat_adj
table_degree$sign <- NULL
table_degree$hemi <- NULL
table_degree <- table_degree[, c(1, 3, 4, 2)]
colnames(table_degree)[which(names(table_degree) == "pval")] <- "pval_adjusted"
print(sum1$fstatistic[2]) #add to table manually
print(sum1$fstatistic[3]) #add to table manually

table_participation <- rstat_adj
table_participation$sign <- NULL
table_participation$hemi <- NULL
table_participation <- table_participation[, c(1, 3, 4, 2)]
colnames(table_participation)[which(names(table_participation) == "pval")] <- "pval_adjusted"
print(sum1$fstatistic[2]) #add to table manually
print(sum1$fstatistic[3]) #add to table manually

#save
names <- list('Regional strength' = table_strength, 'Regional degree' = table_degree, 'Regional participation' = table_participation)
openxlsx::write.xlsx(names, file = 'supplementary_table_01_anxiety_subgroups.xlsx')

names <- list('Regional strength' = table_strength, 'Regional degree' = table_degree)
openxlsx::write.xlsx(names, file = 'supplementary_table_02_mania_subgroups.xlsx')
