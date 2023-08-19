rm(list = ls())
library(data.table) #载入数据用
library(tidyverse)
#表达矩阵
gtex.tpm=fread("data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz",header = T, sep = '\t',data.table = F)
rownames(gtex.tpm)=gtex.tpm[,1]
df <- gtex.tpm[1:5,1:5]
#########
annotat= gtex.tpm[,c(1,2),drop=F]   ###注释
rownames(annotat) <- annotat$Name
colnames(annotat)[2] <- "gene"
###############
exp_gtex.tpm=gtex.tpm[,-1]

#样本信息
data_cl=fread("data/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt",header = T, sep = '\t',data.table = F)
data_cl=data_cl[,c(1,6)]
names(data_cl)=c('Barcode','Tissue')
data_cl=data_cl[data_cl$Tissue == 'Ovary',] #筛选出ovarian的数据

#筛选，筛选之后还剩100个barcode
exp_gtex.tpm=exp_gtex.tpm[,colnames(exp_gtex.tpm) %in% data_cl$Barcode]
#还原为TPM
range(exp_gtex.tpm)
exp_gtex.tpm= log2(exp_gtex.tpm+1)

#基因注释
exp_gtex.tpm=as.matrix(exp_gtex.tpm)
t_index=intersect(rownames(exp_gtex.tpm),rownames(annotat)) #行名取交集，t_index中是能够进行注释的probe_id
exp_gtex.tpm=exp_gtex.tpm[t_index,]
annotat=annotat[t_index,]
rownames(exp_gtex.tpm)=annotat$gene

#去除重复基因名
t_index1=order(rowMeans(exp_gtex.tpm),decreasing = T)
t_data_order=exp_gtex.tpm[t_index1,]
keep=!duplicated(rownames(t_data_order))#对于有重复的基因，保留第一次出现的那个，即行平均值大的那个
exp_gtex.tpm=as.data.frame(t_data_order[keep,]) #得到最后处理之后的表达谱矩阵
range(exp_gtex.tpm)
##################################################################################
load(file = "data/TCGA_OV_TPM.1.Rdata")
range(ov.exprSet)
ov.exprSet <- log2(ov.exprSet+1)
comgene <- intersect(rownames(exp_gtex.tpm),rownames(ov.exprSet))
dataset=c("GTEx",
          "TCGA")
combined.expr <- cbind.data.frame(exp_gtex.tpm[comgene,],
                                  ov.exprSet[comgene,])
library(cluster)
library(oompaBase)
library(sva)
source("resource/batchPCA.R")
library(paletteer)
# 绘制PCA散点图，检查批次效应
batchPCA(indata = t(scale(t(combined.expr))),
         batch= rep(dataset, times = c(ncol(exp_gtex.tpm),
                                       ncol(ov.exprSet))),
         fig.dir = ".",
         PCA.fig.title = "temp/TCGA",
         cols =paletteer_c("grDevices::Spectral", 2),
         showID = F,
         cex = 0.7,
         showLegend = T) # 

# 去除批次效应
batch <- data.frame(batch = rep(c("GTEx",
                                  "TCGA"), 
                                times = c(ncol(exp_gtex.tpm),ncol(ov.exprSet))))

modcombat = model.matrix(~1, data = batch)
combined.expr.combat <- as.data.frame(ComBat(dat=as.matrix(combined.expr), batch=batch$batch, mod=modcombat))

batchPCA(indata = t(scale(t(combined.expr.combat))),
         batch= rep(dataset, 
                    times = c(ncol(exp_gtex.tpm),ncol(ov.exprSet))),
         fig.dir = ".",
         PCA.fig.title = "temp/combined TCGA",
         cols =paletteer_c("grDevices::Spectral", 2),
         showID = F,
         cex = 0.7,
         showLegend = T) 
####################
exprSet <- combined.expr.combat
load(file = "temp/ndegs.Rdata")
gene <- data.table::fread("resource/ECM.txt",data.table = F)
gene <- gene[gene$`Gene Symbol`%in% ndegs$Gene_ID,]

co.gene <- intersect(gene$`Gene Symbol`,rownames(exprSet))

exprSet.1 <- as.data.frame(t(exprSet[gene$`Gene Symbol`,]))
exprSet.1 <- exprSet.1  %>% 
  mutate(Group = substring(rownames(.),1,4)) %>% 
  select(Group,everything())
mydata<- exprSet.1 %>% 
  gather(key="gene",value="Expression",2:ncol(exprSet.1)) %>% 
  dplyr::select(gene,Expression,everything()) 
colnames(mydata)[3] <- "Group"

range(mydata$Expression)
####计算统计学意义
pvalues <- sapply(unique(mydata$gene), function(x) {
  res <- wilcox.test(Expression ~ Group, data = subset(mydata, gene == x))$p.value 
})
pv <- data.frame(gene = unique(mydata$gene), 
                 pvalue = pvalues)

pv$sigcode <- cut(pv$pvalue, c(0, 0.001, 0.01, 0.05, 1), 
                  labels=c('***', '**', '*', 'ns'))
mydata$Group <- factor(mydata$Group,
                       levels = c("GTEX","TCGA"),  
                       labels = c("GTEX_N","TCGA_C"),
                       ordered = F )
mydata$gene <- factor(mydata$gene,
                      levels = unique(gene$`Gene Symbol`),  
                      ordered = F )

table(gene$Category)

library(ggplot2)
p_bot <- ggplot(mydata, aes(gene, Expression, fill=Group)) + 
  geom_boxplot(aes(col = Group)) + 
  scale_color_manual(values = c("#C75DAA","#0FCFC0")) + # 设置透明色
  scale_fill_manual(values = c("#C75DAA","#0FCFC0")) +
  geom_text(aes(gene, y=max(mydata$Expression) * 1.1, 
                label=sigcode),
            data=pv, inherit.aes=F) +
  #annotate("text", x = 2, y= -0.35, label = "wilcox.test",size=4)+
  #labs(fill = "")+
  xlab(NULL)+
  ylab("log2(TPM+1)")+
  
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background=element_rect(fill='transparent',color='black'),
        axis.line = element_line(colour = "black"),
        legend.position = c(0.5,0.11),#"top",   #c(0.5,1)
        legend.direction = "horizontal",
        legend.background = element_blank(),
        legend.title = element_blank(),
        #axis.title.x = element_text(face = "bold",size = 12),
        axis.title.y=element_text(face = "bold",size = 11),     
        axis.text.y = element_text(face = "bold",size = 10),
        axis.text.x = element_text(face = "bold",hjust = 1,
                                   color=c(rep("#D53E4FFF",16),rep("#F46D43FF",9),rep("#FDAE61FF",5)
                                           ,rep("#9ECAE1",5),rep("#4292C6",13),rep("#08519C",23)),
                                   vjust = 0.5, angle = 90, size = 10)) 
# 用白色标记箱子的基本统计量
dat <- ggplot_build(p_bot)$data[[1]]
p_bot.1 <- p_bot + geom_segment(data=dat, aes(x=xmin, xend=xmax, y=middle, yend=middle), color="white", inherit.aes = F)
p_bot.1
###############################################################################
library(export)
graph2ppt(file="out/1.1.1 TCGA and GTEx箱式图.ppt")
#################################pca
library(plyr)
library(dplyr)
expr_df <- exprSet.1[,-1]
meta_df <- exprSet.1[,1,drop=F]
#用`prcomp`进行PCA分析
pca.results <- prcomp(expr_df, center = TRUE, scale. = FALSE)
###
pca.rotation <- pca.results$rotation
pca.rotation
pca.pv <- summary(pca.results)$importance[2,]
pca.pv
low_dim_df <- as.data.frame(pca.results$x[,c(1,2)])
low_dim_df$group <- meta_df$Group
#查看前3行
low_dim_df[1:3,]
add_ellipase <- function(p, x="PC1", y="PC2", group="group",
                         ellipase_pro = 0.95,
                         linetype="dashed",
                         colour = "black",
                         lwd = 2,...){
  obs <- p$data[,c(x, y, group)]
  colnames(obs) <- c("x", "y", "group")
  ellipse_pro <- ellipase_pro
  theta <- c(seq(-pi, pi, length = 50), seq(pi, -pi, length = 50))
  circle <- cbind(cos(theta), sin(theta))
  ell <- ddply(obs, 'group', function(x) {
    if(nrow(x) <= 2) {
      return(NULL)
    }
    sigma <- var(cbind(x$x, x$y))
    mu <- c(mean(x$x), mean(x$y))
    ed <- sqrt(qchisq(ellipse_pro, df = 2))
    data.frame(sweep(circle %*% chol(sigma) * ed, 2, mu, FUN = '+'))
  })
  names(ell)[2:3] <- c('x', 'y')
  
  ell <- ddply(ell, .(group) , function(x) x[chull(x$x, x$y), ])
  p <- p + geom_polygon(data = ell, aes(x=x,y=y,group = group), 
                        colour = colour,
                        alpha = 1,fill = NA,
                        linetype=linetype,
                        lwd =lwd)
  return(p)
}
#计算坐标轴标签
pc1.pv <- paste0(round(pca.pv['PC1'],digits = 3) * 100, "%")
pc2.pv <- paste0(round(pca.pv['PC2'],digits = 3) * 100, "%")

#画出各个样本在二维空间的点
p <- ggplot(low_dim_df) + 
  geom_point(aes(x=PC1, y=PC2, color=group), size=2.5, #点的大小
             shape=20,#点的形状
             alpha=0.5) +#设置点为半透明，出现叠加的效果
  #如果使用默认的颜色，就在下面这行前面加个#
  scale_color_manual(values = c("#C75DAA","#0FCFC0")) +
  #还能调整整体的颜色亮度
  #scale_colour_hue(l=45) + 
  theme_bw() + #去除背景色
  theme(panel.grid =element_blank()) + #去除网格线
  xlim(-30,30)+
  #添加标签，同样可以加到方法一的同一位置
  annotate("text",x=21,y=28,label = "GTEx_ovarian",color = "#C75DAA") +
  annotate("text",x=-21,y=14,label = "TCGA_OV",color = "#0FCFC0") +
  
  #图例
  guides(color=guide_legend(title = NULL)) +
  theme(legend.background = element_blank(), #移除整体边框
        #图例的左上角置于绘图区域的左上角
        legend.position = 'none',#c(0,1),
        #legend.justification = c(0,1),
        axis.title = element_text(face = "bold",size = 11),
        axis.text = element_text(size=10,face = "bold")) + #字体大小
  
  #调整坐标轴标签
  xlab(paste0("PC1 ( ", pc1.pv," variance )")) + 
  ylab(paste0("PC2 ( ", pc2.pv," variance )")) 
p
#画圈圈
p1 <- add_ellipase(p,ellipase_pro = 0.95,colour = "dimgrey",linetype=2,lwd=0.8)
p1
graph2ppt(file="out/1.1.2 TCGA and GTEx pca.ppt")
