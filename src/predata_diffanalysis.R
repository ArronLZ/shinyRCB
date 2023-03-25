fpkmToTpm <- function(fpkm) {
  # fpkm必须为非log化数据
  exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}

Xena.process <- function(df, annot, df.log=T, TPM=FALSE, type='mRNA') {
  # df 可为counts或fpkm, 默认未counts, 但若TPM若为T，则df必须是FPKM ！！！！
  # df.log表示df是否为已log2化后的数据
  #df: 为数据框 
  #          Ensembl_ID TCGA-97-7938-01A TCGA-55-7574-01A TCGA-05-4250-01A
  #1 ENSG00000000003.13        10.989394         9.967226        12.386940
  #2  ENSG00000000005.5         4.000000         0.000000         2.584963
  #3 ENSG00000000419.11        10.253847         9.541097        11.501340
  #annot: 为csv： (前三列必须为id, symbol：注释名称, type：选择类型)
  #      ensembl_gene_id   symbol  type
  #    1 ENSG00000121410   A1BG   lncRNA
  #    2 ENSG00000148584   A1CF   mRNA
  #    3 ENSG00000175899   A2M    mRNA
  #type='mRNA' or 'lncRNA'
  #TPM若为T，则需要df是FPKM
  #output:
  # 注释symbol, 均数去重复基因，并挑选出mRNA, 返回列名为样本，行名为symbol.
  # dplyr, limma

  # 1. 数据是否为log数据 
  row.index <- df[,1]
  if (df.log) {
    df  <- apply(df[, 2:ncol(df)], 2, function(x) 2^x - 1)
  } else {
    df  <- as.matrix(df[, 2:ncol(df)])
  }
  # 2. 是否将fpkm转为tpm
  if (TPM) {
    df <- apply(df, 2, fpkmToTpm)
  }
  rownames(df)  <- row.index
  df[1:4,1:3]
  # 3. 注释
  df  <- df[intersect(rownames(df), annot[,1]), ]
  rownames(df) <- annot[match(rownames(df), annot[,1]), 2]

  # 4. 去重
  df <- avereps(df)
  df <- data.frame(df, check.names = F)
  # 返回结果
  if (type == 'all') {
    df  <- df
  } else if (type == 'mRNA') {
    df  <- df[rownames(df) %in% annot[annot[,3] == "mRNA", 2], ]
  } else if (type == 'lncRNA') {
    df  <- df[rownames(df) %in% annot[annot[,3] == "lncRNA", 2], ]
  }
  return(df)
}

DEG.edgeR <- function(exprset.group, pval=0.05,fdr=0.1, logfc=1) {
  # 此函数为edegR差异分析，返回两个表格组成的list：
  #   1.全部的分析表格 + 2.差异表格(P<0.05&FDR<0.1&abs(et$log2FC)>1)
  
  # exprset.group要求为list (主要为了方便lapply函数做批量差异分析)
  # exprset.group 由下列三个数据打包而成，可用代码示例：
  #   exprset.group <- list(eset=exprset, group=phen, f_mark=c())
  
  #** exprset.group$eset：矩阵数据格式(数值型，整型)
  #       | row1 | row2 | row3  | row4
  # gene1 |  34  |  23  |  56   |  23
  # gene2 |  35  |  23  |  12   |  23
  # gene3 |  12  |  78  |  78   |  78
  # ** exprset.group$group：分组数据格式：
  #   需要组的行名=表达谱的列名 rownames(group) == colname(eset)
  #            type
  # rowname1 | tumor
  # rowname2 | tumor
  # rowname3 | normal
  # rowname4 | normal
  # ** exprset.group$f_mark：list标记项，可以为空，但是如果是批量差异分析，
  #     建议必须要有此参数，方便后续对个结果进行保存和导出
  
  # 1.0 count，group
  cat('\n', exprset.group$f_mark, ' ================\n')
  exprset <- exprset.group$eset
  pheno <- exprset.group$group
  pheno$type
  
  # 1.1 构建DEGList
  y <- DGEList(counts=exprset, group=pheno$type)
  keep <- filterByExpr(y)
  y <- y[keep, , keep.lib.sizes=FALSE]
  # 1.2 计算离散因子
  y <- calcNormFactors(y)
  design <- model.matrix(~pheno$type)
  y <- estimateDisp(y, design)
  #To perform quasi-likelihood F-tests:
  fit <- glmQLFit(y, design)
  qlf <- glmQLFTest(fit, coef=2)
  #topTags(qlf)
  
  # 获取结果
  # 获取排名靠前的基因，这里设置n=80000是为了输出所以基因
  et <- topTags(qlf, n=80000)
  et <- as.data.frame(et) # 转换为数据框类型
  et <- cbind(rownames(et), et)  # 将行名粘贴为数据框的第一列
  colnames(et) <- c("Gene", "log2FC", "log2CPM", "F", "PValue", "FDR")
  et <- et %>% dplyr::arrange(desc(log2FC), PValue)
  
  # 差异基因筛选
  etSig <- et[which(et$PValue < pval & et$FDR < fdr & abs(et$log2FC) > logfc),]
  return(list(qlf=qlf, resdf=et, deg=etSig))
}

DEG.DESeq2 <- function(exprset.group) {
  # ** 矩阵数据格式(数值型，整型)
  #       | row1 | row2 | row3  | row4
  # gene1 |  34  |  23  |  56   |  23
  # gene2 |  35  |  23  |  12   |  23
  # gene3 |  12  |  78  |  78   |  78
  # ** 分组数据格式：需要组的行名=表达谱的列名 rownames(group) == colname(eset)
  #            type
  # rowname1 | tumor
  # rowname2 | tumor
  # rowname3 | normal
  # rowname4 | normal
  
  # 1.0 count
  #keepGene <- rowSums(edgeR::cpm(exprset)>0) >=2
  #exprset <- exprset[keepGene,]
  # 1.1 group
  cat('\n', exprset.group$f_mark, ' ================\n')
  exprset <- exprset.group$eset
  pheno <- exprset.group$group
  design_type <- as.formula(paste('~ ', colnames(pheno)[1]))
  # 1.2 dds
  dds <- DESeqDataSetFromMatrix(countData = exprset,
                                colData = pheno,
                                design = design_type)
  # 1.3
  dds <- DESeq(dds, parallel = T)
  

  ## 数据转换
  keep <- rowSums(counts(dds) >= 15) >= 4  #过滤低表达基因，至少有3个样品都满足10个以上的reads数  
  dds <- dds[keep, ]
  # resultsNames(dds)                # 查看结果的名称。
  #type1 <- levels(dds$type)[1]
  #type2 <- levels(dds$type)[2]
  # 如不指定contrast，默认第二因子比第一因子。建议将对照的因子水平放在第一个。
  # 也可以通过results(dds, contrast=c("type", type2, type1)): type2 比 type 1
  res <- results(dds) 
  # summary(res)      #看一下结果的概要信息，p值默认小于0.1。
  
  ### 差异分析总表
  res_df <- data.frame(res) 
  eset_norma <- data.frame(counts(dds, normalized=TRUE), check.names = F)
  res_df <- merge(res_df, eset_norma, by="row.names", sort=FALSE)
  res_df <- res_df %>% 
    dplyr::rename(Gene = 'Row.names',
                  log2FC = log2FoldChange,
                  PValue = pvalue,
                  FDR=padj) %>% 
    dplyr::arrange(desc(log2FC), PValue)
  rownames(res_df) <- res_df$Gene
  
  # 差异基因筛选
  etSig <- res_df %>% filter(PValue < 0.05 & FDR < 0.1 & abs(log2FC) > 1)
  return(list(dds=dds, resdf=res_df, deg=etSig))
}

DEG.voom <- function(exprset.group) {
  cat('\n', exprset.group$f_mark, ' ================\n')
  exprset <- exprset.group$eset
  pheno <- exprset.group$group
  
  #pheno$type
  design <- model.matrix(~0 + pheno$type)
  colnames(design) <- levels(pheno$type)
  colnames(design) <- c('con','trt')  # 第二因子比第一因子
  #contrast.matrix <- makeContrasts(group2-group1, group3-group2, group3-group1, levels=design)
  contrast.matrix <- makeContrasts(trt-con, levels=design)  # 第二因子比第一因子
  
  # edgeR
  y <- DGEList(counts=exprset)
  keep <- filterByExpr(y, design, min.count = 15, min.total.count = 20, 
                       large.n = 10, min.prop = 1)
  y <- y[keep, , keep.lib.sizes=FALSE]
  y <- calcNormFactors(y)
  # voom
  v <- voom(y, design, plot=TRUE)
  fit <- lmFit(v, design)
  fit <- contrasts.fit(fit, contrast.matrix)
  fit <- eBayes(fit)
  #topTable(fit, coef=ncol(design))
  tab <- topTable(fit, sort.by = "P", n = Inf)
  tab <- cbind(rownames(tab), tab)  # 将行名粘贴为数据框的第一列
  colnames(tab) <- c("Gene", "log2FC", "AveExp", "T", "PValue", "FDR", "B")
  tab <- tab %>% dplyr::arrange(desc(log2FC), PValue)
  
  deg <- tab[which(tab$PValue < 0.05 & tab$FDR < 0.1 & abs(tab$log2FC) > 1),]
  return(list(resdf=tab, deg=deg))
}

## 自定义函数，go分析，并画go关系图，保存golist对象
DEG_GO <- function(genelist, orgdb="org.Hs.eg.db", sigNodes=20, resultdir, filemark) {
  #genelist为数据框，需有一列为entrzid
  #resultdir <- 'result_stringtie/p005fc15',输出文件夹
  #filemark <- 'sh2118'，特定文件名
  go_list <- list()
  for (ont in c('ALL', 'BP', 'CC', 'MF')) {
    GO <- enrichGO(genelist[, "ENTREZID"],#GO富集分析BP模块
                   OrgDb = orgdb,
                   keyType = "ENTREZID",
                   ont = ont,
                   minGSSize = 1,
                   pvalueCutoff = 1,
                   pAdjustMethod = "BH",
                   qvalueCutoff = 1,
                   readable = T)
    if (ont=='ALL') {
      go_list[[ont]] <- GO
      go_all_df <- data.frame(GO)
      write.table(go_all_df, sep = '\t', quote = F, 
                  file = paste0(resultdir, '/GOall', '_',filemark, '.txt'))
    } else {
      go_list[[ont]] <- GO
      # plot关系图
      fname <- paste0(resultdir, '/', ont, '_',filemark, '.pdf')
      cat(fname)
      pdf(fname, height = 8, width = 12)
      plotGOgraph(GO, firstSigNodes = sigNodes)
      dev.off()
    }
  }
  save(go_list, file = paste0(resultdir, '/GOall', '_',filemark, '.Rdata'), compress = T)
  return(go_list)
}

## 自定义函数，kegg分析，保存kegg对象（原y叔函数实为p.adjust）,人为取出p<0.05的数据
DEG_KEGG <- function(genelist, orgdb="org.Hs.eg.db", org_kegg='hsa',resultdir, filemark) {
  #genelist为数据框，需有一列为entrzid
  #resultdir <- 'result_stringtie/p005fc15',输出文件夹
  #filemark <- 'sh2118'，特定文件名
  KEGG <- enrichKEGG(genelist[, "ENTREZID"],#KEGG富集分析
                     organism = org_kegg, minGSSize = 1,
                     pvalueCutoff = 1, qvalueCutoff = 1)
  KEGG <- setReadable(KEGG, OrgDb = orgdb, keyType="ENTREZID")
  keggdf <- filter(KEGG@'result', pvalue < 1)
  save(KEGG, file = paste0(resultdir, '/KEGGpathway', '_', filemark, '.Rdata'), compress = T)
  write.table(keggdf, sep = '\t', quote = F, 
              file = paste0(resultdir, '/KEGG', '_', filemark, '.txt'))
  return(list(KEGG=KEGG, pSigDF=keggdf))
}
