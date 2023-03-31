#  生存ESCA ESCC ESAD的表达临床信息混合标准表
# 标准： 1.csv格式 2. 第一列为样本名，不设行名
#        3. 第2-3列为OS, OS.time, 第4-12列为临床信息，若不满12构建空列补足
#        4. 第13开始为基因表达矩阵，未log化
# sample_id        OS OS.time Age Gender Stage  T  N  M Hist sample_type batch_num     RAB4B
# TCGA-2H-A9GF-01A  1     784  67   Male   III T3 N1 M0 ESCA       Tumor  382.51.0  4.332118
# TCGA-2H-A9GG-01A  1     610  66   Male   III T3 N1 M0 ESCA       Tumor  382.51.0 10.743245
# TCGA-2H-A9GH-01A  1     951  44   Male    II T1 N1 M0 ESCA       Tumor  382.51.0 11.941258
# TCGA-2H-A9GI-01A  1     435  68   Male   III T3 N1 M0 ESCA       Tumor  382.51.0  6.266653
rm(list = ls());gc()
library(data.table)
library(tidyverse)
source('E:/OneDrive/Desktop/DM.R1/v5/0.fun/p1.fun.R')

# 注意导入 ##
#load('E:/OneDrive/Desktop/DM.R1/v5/0.fun/data.flow/xenaLUAD.env.RData')
self <- new.env(parent = baseenv())
load("E:/BaiduSyncdisk/Public_TCGA/UCSCxena/annot.2/gencode.v22.annot.RData", 
     envir = self)

ESCAfpkm <- fread('E:/Public_data/TCGA_XENA/ESCA/eset/TCGA-ESCA.htseq_fpkm.tsv.gz',
                  data.table = F)
ESCAfpkm[1:4,1:4]
survdf <- fread('E:/Public_data/TCGA_XENA/ESCA/clin/TCGA-ESCA.survival.tsv',
                data.table = F)
phen <- fread('E:/Public_data/TCGA_XENA/ESCA/clin/xena.esca.csv', data.table = F)
phendf <- phen[,1:10]
# surv
survdf <- Xena.surv(survdf)
survdf %>% head
survdf$OS.time %>% range()
phendf[,1:4] %>% head

survdf <- merge(survdf, phendf, by='sample_id', all = T)

# # #
# =========================
# fpkm, 提取mRNA,并将FPKM转为TPM
TPMdf <- Xena.process(ESCAfpkm, annot = self$all_anot, TPM = T)
# TPMdf.lnc <- Xena.process(fpkmdf, annot = self$all_anot, type = 'lncRNA', TPM = T)
TPMdf[1:4,1:4]
survdf[1:4,]
eset_clin <- Xena.mergeEC(df = TPMdf, clic = survdf) 
eset_clin[1:5,1:5]
PR.checkna(eset_clin[,1:15])
eset_clin[1:4,1:15]

# ========================
FPKMdf <- Xena.process(ESCAfpkm, annot = self$all_anot, TPM = F)
eset_clin.fpkm <- Xena.mergeEC(df = FPKMdf, clic = survdf)
write.csv(eset_clin, file = 'DataClean_ESCA.tpm&clic.csv',
          row.names = F)
write.csv(eset_clin.fpkm, file = 'DataClean_ESCA.fpkm&clic.csv',
          row.names = F)

###
eset_clin.escc <- eset_clin %>% filter(Hist %in% c("ESCC", "ESCC_keratinizing"))
eset_clin.esad <- eset_clin %>% filter(Hist == "ESCA")

eset_clin.fpkm.escc <- eset_clin.fpkm %>% filter(Hist %in% c("ESCC", "ESCC_keratinizing"))
eset_clin.fpkm.esad <- eset_clin.fpkm %>% filter(Hist == "ESCA")
write.csv(eset_clin.escc, file = 'DataClean_ESCC.tpm&clic.csv',
          row.names = F)
write.csv(eset_clin.fpkm.escc, file = 'DataClean_ESCC.fpkm&clic.csv',
          row.names = F)

write.csv(eset_clin.esad, file = 'DataClean_ESAD.tpm&clic.csv',
          row.names = F)
write.csv(eset_clin.fpkm.esad, file = 'DataClean_ESAD.fpkm&clic.csv',
          row.names = F)