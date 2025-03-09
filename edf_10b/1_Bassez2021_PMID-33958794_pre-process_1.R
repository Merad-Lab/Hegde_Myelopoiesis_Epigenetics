# Load Seurat object

# Cohort 1
mypath <- '../yourpath/1863-counts_cells_cohort1.rds'
df <- readRDS(mypath)
df <- as.matrix(df)
write.csv(df, "../yourpath/cells_cohort1_matrix.csv")

# Cohort 2
mypath <- '../yourpath/1867-counts_cells_cohort2.rds'
df <- readRDS(mypath)
df <- as.matrix(df)
write.csv(df, "../yourpath/cells_cohort2_matrix.csv")
