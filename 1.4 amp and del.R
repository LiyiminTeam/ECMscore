rm(list = ls())
library(ComplexHeatmap) # 用于绘制热图
library(circlize) # 用于热图颜色设置
library(data.table) # 用于读取大文件
Sys.setenv(LANGUAGE = "en") #显示英文报错信息
options(stringsAsFactors = FALSE) #禁止chr转成factor
###导入数据
load(file = "data/TCGA_OV_TPM.Rdata")
expr <- ov.exprSet
###
load(file = "temp/ndegs.Rdata")
gene <- data.table::fread("resource/ECM.txt",data.table = F)
gene <- gene[gene$`Gene Symbol`%in% ndegs$Gene_ID,]

gene$Category <- dplyr::recode(gene$Category, 
                               'ECM Glycoproteins'="Glycoproteins",
                               'ECM-affiliated Proteins'="Affiliated",
                               'ECM Regulators'="Regulators" )

gene$Category <- factor(gene$Category,
                        levels = unique(gene$Category),  
                        ordered = F )
###
expr <- expr[gene$`Gene Symbol`,]

cna <- read.table("data/TCGA.OV.sampleMap_Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes.gz",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
cna$gene <- sapply(strsplit(rownames(cna),"|",fixed = T),"[",1)
cna <- cna[!duplicated(cna$gene),]; cna <- cna[,setdiff(colnames(cna),"gene")]
is.element(gene$`Gene Symbol`,rownames(cna)) # 有些基因没有对应的拷贝数结果

cna <- cna[intersect(rownames(cna),gene$`Gene Symbol`),]
cna[cna > 1] <- 1 # 统一扩增
cna[cna < -1] <- -1 # 统一缺失

## 提取共同样本
comsam <- intersect(colnames(expr), colnames(cna))

expr <- expr[,comsam]
cna <- cna[,comsam]
##################################################
ampFreq <- delFreq <- 
  as.data.frame(matrix(NA,
                       nrow = nrow(gene),
                       ncol = 1, 
                       dimnames = list(gene$`Gene Symbol`, 
                                       "OV")))


## 扩增/缺失频率
for (i in gene$`Gene Symbol`) {
  if(!is.element(i, rownames(cna))) { # 同理，如果存在拷贝数中缺失某基因，则保持NA
    ampFreq[i,] <- NA 
    delFreq[i,] <- NA
  } else { # 否则
    # 计算i在总样本中的频率
    ampFreqInAll <- sum(as.numeric(cna[i,]) == 1)/ncol(cna) # 总样本中扩增的数目除以总样本数
    delFreqInAll <- sum(as.numeric(cna[i,]) == -1)/ncol(cna) # 总样本中缺失的数目除以总样本数
    ampFreq[i, 1] <- ampFreqInAll
    delFreq[i, 1] <- delFreqInAll
  }
}
##################################画图
my.col= c("#D53E4FFF","#F46D43FF","#FDAE61FF","#9ECAE1","#4292C6","#08519C")

library(gplots)
ha1 <- HeatmapAnnotation(Category=gene$Category,
                         col = list(Category = c("Glycoproteins"= my.col[1],
                                                 "Collagens"=my.col[2],
                                                 "Proteoglycans" = my.col[3],
                                                 "Affiliated" = my.col[4],
                                                 "Regulators"=my.col[5],
                                                 "Secreted Factors"=my.col[6])
                         ),
                         simple_anno_size = unit(0.25, "cm"),
                         show_annotation_name =F, ## 是否显示注释的标题. 
                         show_legend          = F, # 不显示亚型的图例，因为一目了然
                         which = "col")

############
range(na.omit(ampFreq))
my_color.2=c(colorpanel(80,low="white",high="#0FCFC0"))

hm.ampFreq <- pheatmap(as.matrix(t(ampFreq)),
                          border_color = "grey80", # 热图单元格无边框
                          #rect_gp = gpar(col = "grey80"), 
                          top_annotation = ha1, 
                          #bottom_annotation = ha2,
                          #left_annotation = ha1,
                          #row_names_side="left",
                          color = my_color.2,
                          #row_names_side     = "left",
                          column_names_side  = "top",
                          show_rownames = F, # 显示行名
                          show_colnames = T, # 不显示列名
                          cellheight = 10, # 热图高度固定
                          cellwidth = 10, # 热图宽度固定
                          name = "Amplification\nFrequency", # 图例名字
                          #angle_col = 180,
                          column_split = gene$Category,
                          #heatmap_legend_param = list(direction = "horizontal"),
                          #row_split= pv$Type,
                          #gaps_col = cumsum(table(pv$Type)),#743, # 列分割cumsum(table(cafs.group$CAF_Scores))
                          #gaps_row = cumsum(table(pv$Type)),#cumsum(table(pathway.group$Group))
                          cluster_rows = F, # 行不聚类
                          cluster_cols = F) # 列不聚类

range(na.omit(delFreq))
my_color.3=c(colorpanel(80,low="white",high="#C75DAA"))
hm.delFreq <- pheatmap(as.matrix(t(delFreq)),
                       border_color = "grey80", # 热图单元格无边框
                       #rect_gp = gpar(col = "grey80"), 
                       #top_annotation = ha1, 
                       #bottom_annotation = ha2,
                       #left_annotation = ha1,
                       #row_names_side="left",
                       color = my_color.3,
                       #row_names_side     = "left",
                       column_names_side  = "top",
                       show_rownames = F, # 显示行名
                       show_colnames = F, # 不显示列名
                       cellheight = 10, # 热图高度固定
                       cellwidth = 10, # 热图宽度固定
                       name = "Deletion\nFrequency", # 图例名字
                       #angle_col = 180,
                       column_split = gene$Category,
                       #heatmap_legend_param = list(direction = "horizontal"),
                       #row_split= pv$Type,
                       #gaps_col = cumsum(table(pv$Type)),#743, # 列分割cumsum(table(cafs.group$CAF_Scores))
                       #gaps_row = cumsum(table(pv$Type)),#cumsum(table(pathway.group$Group))
                       cluster_rows = F, # 行不聚类
                       cluster_cols = F) # 列不聚类
#pdf(file = "complexheatmap of immunomodulator.pdf", width = 8,height = 12)
draw(hm.ampFreq %v% hm.delFreq, # 水平衔接各个子热图
     heatmap_legend_side = "right") # 热图颜色图例显示在下方
#invisible(dev.off())
library(export)
graph2ppt(file="out/1.4 amp and del.pptx")
