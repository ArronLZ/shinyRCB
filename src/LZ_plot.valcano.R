suppressMessages(library(scales))
suppressMessages(library(ggsci))
suppressMessages(library(ggrepel))
suppressMessages(library(ggpubr))
# show_col(pal_npg()(9))

DEGplot.volcano <- function(result, logFC = 1.5, adj_P = 0.1, label_geneset = NULL) {
# result必须为数据框，必须具有logFolodChange, adj列，如果需要在图上标记，请设置label_geneset.
  result$Type = "NONE"
  result$Type[which(result$log2FoldChange > logFC & result$padj < adj_P)] = "UP"
  result$Type[which(result$log2FoldChange < (-logFC) & result$padj < adj_P)] = "DOWN"
  xlim = max(abs(result$log2FoldChange))
  if(is.null(label_geneset)) {
    p = ggplot(result, aes(x = log2FoldChange, y = -log10(padj)))+
      geom_point(data = result, aes(x = log2FoldChange, y = -log10(padj), color = Type)) +
          theme_bw()+
    geom_vline(xintercept = c(-logFC, logFC), lty = 2)+
    geom_hline(yintercept = c(-log10(adj_P)), lty = 2)+
    scale_x_continuous(limits = c(-xlim, xlim))+
    coord_fixed(ratio = ( 2*xlim )/(max(-log10(result$padj), na.rm = T)))+
    xlab("log10FoldChange")+
    ylab("-log10P-value")+
    scale_color_manual(values = c("UP" = "#E64B35FF", "DOWN" = "#3C5488FF", "NONE" = "grey"))+
    theme(panel.grid = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text = element_text(color = "black"))
  } else {
    p = ggplot(result, aes(x = log2FoldChange, y = -log10(padj)))+
      geom_point(data = result, 
                 aes(x = log2FoldChange, y = -log10(padj), color = Type), alpha=0.9)+
      geom_point(data = result[which(result$Gene %in% label_geneset),],
                 aes(x = log2FoldChange, y = -log10(padj)),color = "black",size = 4)+
      geom_point(data = result[which(result$Gene %in% label_geneset),],
                 aes(x = log2FoldChange, y = -log10(padj)),color = "white",size = 2.5)+
      geom_point(data = result[which(result$Gene %in% label_geneset),],
                 aes(x = log2FoldChange, y = -log10(padj),color = Type),size = 1.5)+
      geom_text_repel(data = result[which(result$Gene %in% label_geneset),],
                      aes(x = log2FoldChange, y = -log10(padj), label = Gene), 
                      fontface = "bold") + #fontface = "italic"
          theme_bw()+
    geom_vline(xintercept = c(-logFC, logFC), lty = 2)+
    geom_hline(yintercept = c(-log10(adj_P)), lty = 2)+
    scale_x_continuous(limits = c(-xlim, xlim))+
    coord_fixed(ratio = ( 2*xlim )/(max(-log10(result$padj), na.rm = T)))+
    xlab("log10FoldChange")+
    ylab("-log10P-value")+
    scale_color_manual(values = c("UP" = "#E64B35FF", "DOWN" = "#3C5488FF", "NONE" = "grey"))+
    theme(panel.grid = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text = element_text(color = "black"))  
  }
  return(p)
}


