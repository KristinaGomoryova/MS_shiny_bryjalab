# Libraries required
library(dplyr)
library(HGNChelper)
library(tidyr)
library(here)
library(stringr)

########################################## Prepare the data ###############################################

# 2006 RNF43
data_2006.RNF43 <- read.csv(here('data', 'shinyapp_2006-RNF43.csv'))
long_2006.RNF43 <- data_2006.RNF43 %>% select(Gene.names, starts_with("logFC"), starts_with("adjpval"))
long_2006.RNF43$contrast <- "RNF43_vs_ctrl"
long_2006.RNF43$dataset <- "2006.RNF43"
colnames(long_2006.RNF43) <- c("GeneID", "logFC", "padj", "contrast", "dataset")

# 2434_2660
data_2434.2660 <- read.csv(here('data', 'shinyapp_2434-2660.csv'))
long_2434.2660 <- data_2434.2660 %>% select(Gene.name, starts_with("logFC"), starts_with("adjpval"))
long_2434.2660$contrast <- "RIK_vs_ctrl"
long_2434.2660$dataset <- "2434.2660"
colnames(long_2434.2660) <- c("GeneID", "logFC", "padj", "contrast", "dataset")

# 2006 RNF122
data_2006.RNF122 <- read.csv(here('data', 'shinyapp_2006-RNF122.csv'))
long_2006.RNF122 <- data_2006.RNF122 %>% select(Gene.names, starts_with("logFC"), starts_with("adjpval"))
long_2006.RNF122$contrast <- "RNF122_vs_ctrl"
long_2006.RNF122$dataset <- "2006.RNF122"
colnames(long_2006.RNF122) <- c("GeneID", "logFC", "padj", "contrast", "dataset")

# 2006 RNFT2
data_2006.RNFT2 <- read.csv(here('data', 'shinyapp_2006-RNFT2.csv'))
long_2006.RNFT2 <- data_2006.RNFT2 %>% select(Gene.names, starts_with("logFC"), starts_with("adjpval"))
long_2006.RNFT2$contrast <- "RNFT2_vs_ctrl"
long_2006.RNFT2$dataset <- "2006.RNFT2"
colnames(long_2006.RNFT2) <- c("GeneID", "logFC", "padj", "contrast", "dataset")

# 2999
data_2999 <- read.csv(here('data', 'shinyapp_2999.csv'))
long_2999 <- data_2999 %>% select(Gene.names, starts_with("logFC"), starts_with("adjpval"))
colnames(long_2999) <- gsub("adjpval_", "adjpvalX", colnames(long_2999))
colnames(long_2999) <- gsub("logFC_", "logFCX", colnames(long_2999))

long_2999 <- pivot_longer(long_2999, 
             cols = !Gene.names,
             names_to = c(".value", "contrast"),
             names_sep = "X")
long_2999 <- long_2999[, c(1,3,4,2)]
long_2999$dataset <- "2999_TTLL"
colnames(long_2999) <- c("GeneID", "logFC", "padj", "contrast", "dataset")
long_2999$contrast <- stringr::str_replace_all(long_2999$contrast, "\\.", "_vs_") #Escape the dot!

# 3019
data_3019 <- read.csv(here('data', 'shinyapp_3019.csv'))
long_3019 <- data_3019 %>% select(Gene.name, starts_with("logFC"), starts_with("adjpval"))
colnames(long_3019) <- gsub("adjpval_", "adjpvalX", colnames(long_3019))
colnames(long_3019) <- gsub("logFC_", "logFCX", colnames(long_3019))

long_3019 <- pivot_longer(long_3019, 
                          cols = !Gene.name,
                          names_to = c(".value", "contrast"),
                          names_sep = "X")
long_3019 <- long_3019[, c(1,3,4,2)]
long_3019$dataset <- "3019_WntPCP"
colnames(long_3019) <- c("GeneID", "logFC", "padj", "contrast", "dataset")
long_3019$contrast <- stringr::str_replace_all(long_3019$contrast, "\\.", "_vs_") #Escape the dot!

# 3676 has no statistical evaluation, so will not be used here

# 4310 CK1s
data_4310.CK1 <- read.csv(here('data', 'shinyapp_4310-CK1s.csv'))
long_4310.CK1 <- data_4310.CK1 %>% select(Gene.name, starts_with("logFC"), starts_with("adjpval"))
colnames(long_4310.CK1) <- gsub("adjpval_", "adjpvalX", colnames(long_4310.CK1))
colnames(long_4310.CK1) <- gsub("logFC_", "logFCX", colnames(long_4310.CK1))

long_4310.CK1 <- pivot_longer(long_4310.CK1, 
                          cols = !Gene.name,
                          names_to = c(".value", "contrast"),
                          names_sep = "X")
long_4310.CK1 <- long_4310.CK1[, c(1,3,4,2)]
long_4310.CK1$dataset <- "4310_CK1s"
colnames(long_4310.CK1) <- c("GeneID", "logFC", "padj", "contrast", "dataset")
long_4310.CK1$contrast <- stringr::str_replace_all(long_4310.CK1$contrast, "\\.", "_vs_") #Escape the dot!

# 4310 FZDs
data_4310.FZD <- read.csv(here('data', 'shinyapp_4310-FZDs.csv'))
long_4310.FZD <- data_4310.FZD %>% select(name, starts_with("logFC"), starts_with("adjpval"))
colnames(long_4310.FZD) <- gsub("adjpval_", "adjpvalX", colnames(long_4310.FZD))
colnames(long_4310.FZD) <- gsub("logFC_", "logFCX", colnames(long_4310.FZD))

long_4310.FZD <- pivot_longer(long_4310.FZD, 
                              cols = !name,
                              names_to = c(".value", "contrast"),
                              names_sep = "X")
long_4310.FZD <- long_4310.FZD[, c(1,3,4,2)]
long_4310.FZD$dataset <- "4310_FZDs"
colnames(long_4310.FZD) <- c("GeneID", "logFC", "padj", "contrast", "dataset")
long_4310.FZD$contrast <- stringr::str_replace_all(long_4310.FZD$contrast, "\\.", "_vs_") #Escape the dot!

# 4453 DDA
data_4453.DDA <- read.csv(here('data', 'shinyapp_4453-DDA.csv'))
long_4453.DDA <- data_4453.DDA %>% select(Gene.name, starts_with("logFC"), starts_with("padj"))
colnames(long_4453.DDA) <- gsub("padj_", "adjpvalX", colnames(long_4453.DDA))
colnames(long_4453.DDA) <- gsub("logFC_", "logFCX", colnames(long_4453.DDA))

long_4453.DDA <- pivot_longer(long_4453.DDA, 
                              cols = !Gene.name,
                              names_to = c(".value", "contrast"),
                              names_sep = "X")
long_4453.DDA <- long_4453.DDA[, c(1,3,4,2)]
long_4453.DDA$dataset <- "4453_DDA"
colnames(long_4453.DDA) <- c("GeneID", "logFC", "padj", "contrast", "dataset")
long_4453.DDA$contrast <- stringr::str_replace_all(long_4453.DDA$contrast, "\\.", "_vs_") #Escape the dot!

# 4453 DIA
data_4453.DIA <- read.csv(here('data', 'shinyapp_4453-DIA.csv'))
long_4453.DIA <- data_4453.DIA %>% select(Gene.names, starts_with("logFC"), starts_with("padj"))
colnames(long_4453.DIA) <- gsub("padj_", "adjpvalX", colnames(long_4453.DIA))
colnames(long_4453.DIA) <- gsub("logFC_", "logFCX", colnames(long_4453.DIA))

long_4453.DIA <- pivot_longer(long_4453.DIA, 
                              cols = !Gene.names,
                              names_to = c(".value", "contrast"),
                              names_sep = "X")
long_4453.DIA <- long_4453.DIA[, c(1,3,4,2)]
long_4453.DIA$dataset <- "4453_DIA"
colnames(long_4453.DIA) <- c("GeneID", "logFC", "padj", "contrast", "dataset")
long_4453.DIA$contrast <- stringr::str_replace_all(long_4453.DIA$contrast, "\\.", "_vs_") #Escape the dot!

##################################### Merge into long dataset ############################################
# Merge into long dataset
datasets.long <- rbind(long_2006.RNF43,
                       long_2434.2660,
                       long_2006.RNF122,
                       long_2006.RNFT2,
                       long_2999,
                       long_3019,
                       long_4310.CK1,
                       long_4310.FZD,
                       long_4453.DDA,
                       long_4453.DIA)

# In case there are multiple gene names, take the 1st one
datasets.long$GeneID <- sapply(strsplit(datasets.long$GeneID,";"), `[`, 1)


str(datasets.long)

# 'datasets.long <- datasets.long %>%
#   mutate(change = case_when(
#     logFC > 1 & padj < 0.05 ~ "upregulated",
#     logFC < -1 & padj < 0.05 ~ "downregulated",
#     TRUE ~ "ns"
#   ))

##################################### Update the gene names ###############################################
# Get the most recent mapping of gene names
genenames_newest <- getCurrentHumanMap()
save(genenames_newest, file = ("HGNC_genenames_20230920.RData")) #adjust the date

genenames_update <- checkGeneSymbols(datasets.long$GeneID, species = "human", map = genenames_newest, unmapped.as.na = FALSE)
datasets.long$GeneID <- genenames_update$Suggested.Symbol #replace the original with updated IDs

# After adding new datasets, don't forget to always create a new table!
write.table(datasets.long, sep = "\t", quote = FALSE, row.names = FALSE, file = "datasets_long.txt")
