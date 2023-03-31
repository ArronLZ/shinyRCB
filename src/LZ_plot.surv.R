shiny.readeset <- function(filename, gene1.col) {
  # xena的LZST2023数据
  eset_clin <- fread(filename, data.table = F) # "data/tcga/DataClean_LUAD.tpm&os.csv.gz"
  # 只保留01A数据
  rindex <- str_sub(eset_clin[,1], 14, 16) == "01A"
  eset_clin <- eset_clin[rindex, ]
  # 根据病人去重
  eset_clin[, 1] <- str_sub(eset_clin[,1], 1, 12)
  eset_clin <- eset_clin %>% distinct(., sample_id, .keep_all = T)
  # log化
  eset_clin[, gene1.col:ncol(eset_clin)] <- log2(eset_clin[, gene1.col:ncol(eset_clin)] + 1)
  return(eset_clin)
}

xena.surv.cut <- function(eset_os.df, tagetGene="SIN3B", method=2) {
  eset_clin_p <- eset_os.df
  tagetGene.five <- fivenum(eset_clin_p[,tagetGene])
  if (method == "mean") {
    # mean
    eset_clin_p$Group <- as.vector(ifelse(eset_clin_p[,tagetGene] > mean(eset_clin_p[,tagetGene]), "high.mean","low.mean"))
  } else if (method == "upper75 vs low25") {
    # upper 75 vs low 25
    eset_clin_p$Group <- as.vector(ifelse(eset_clin_p[,tagetGene] > tagetGene.five[2], "high75","low25"))
  } else if (method == "median") {
    # median
    eset_clin_p$Group <- as.vector(ifelse(eset_clin_p[,tagetGene] > tagetGene.five[3], "high.median","low.median"))
  } else if (method == "upper25 vs low75") {
    # upper 25 vs low 75
    eset_clin_p$Group <- as.vector(ifelse(eset_clin_p[,tagetGene] > tagetGene.five[4], "high25","low75"))
  } else if (method == "auto cut point") {
    # cut poion
    res.cut <- surv_cutpoint(eset_clin_p, #数据集
                             time = "OS.time", #生存状态
                             event = "OS", #生存时间
                             variables = tagetGene) #需要计算的数据列名
    eset_clin_p$Group <- as.vector(ifelse(eset_clin_p[,tagetGene] > res.cut$cutpoint$cutpoint, "high.point","low.point"))
  } else if (method == "upper25 vs low25") {
    eset_clin_p.low25 <- eset_clin_p[eset_clin_p[,tagetGene] <= tagetGene.five[2], ]
    eset_clin_p.low25$Group <- "low25"
    eset_clin_p.high25 <- eset_clin_p[eset_clin_p[,tagetGene] > tagetGene.five[4], ]
    eset_clin_p.high25$Group <- "high25"
    eset_clin_p <- rbind(eset_clin_p.low25, eset_clin_p.high25) %>% data.frame(check.names = F)
  }
  return(eset_clin_p)
}

xena.surv.getPvale <- function(rt, num.tran=365, main.text='LUAD') {
  # 指定group.name , 默认为"risk"
  rt$OS.time <- rt$OS.time/num.tran
  diff <- survdiff(Surv(OS.time, OS) ~ Group, data = rt)
  fit <- survfit(Surv(OS.time, OS) ~ Group, data = rt)
  pValue <- 1 - pchisq(diff$chisq, 
                       df = (length(levels(as.factor(rt$Group))) - 1)
  )
  pValue <- signif(pValue, 4)
  # plot
  jco <- ggsci::pal_npg()(9)[c(1,4,9)]
  ggs <- ggsurvplot(fit, # 创建的拟合对象
                    data = rt,  # 指定变量数据来源
                    palette = c(jco[1], jco[2]),
                    conf.int = F, # 显示置信区间
                    pval = T, # 添加P值
                    surv.median.line = "hv",  # 添加中位生存时间线
                    risk.table = TRUE, # 添加风险表
                    xlab = "Time in Years", # 指定x轴标签
                    legend = c(0.8,0.8), # 指定图例位置
                    legend.title = main.text, # 设置图例标题
                    #legend.labs = 1:5, # 指定图例分组标签
                    break.x.by = 2, # 设置x轴刻度间距
                    ggtheme = theme_survminer(font.legend = c(11, "plain", "black")),
                    tables.theme = theme_survminer(font.main = 11) )
  ggs$plot <- ggs$plot + 
    scale_y_continuous(expand = c(0, 0)) + 
    scale_x_continuous(expand = c(0, 0)) 
  return(list(pValue=pValue, fit=fit, pic=ggs))
}
