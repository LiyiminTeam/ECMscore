rm(list = ls())
library(limma)
library(ggplot2)
library(tidyverse)
Sys.setenv(LANGUAGE = "en") #显示英文报错信息
options(stringsAsFactors = FALSE) #禁止chr转成factor
#load(file = "data/GEO exp/GSE26712.exp.Rdata") 
load(file = "data/GEO exp/ovarian exp.Rdata") 

load(file = "temp/ndegs.Rdata")
gene <- data.table::fread("resource/ECM.txt",data.table = F)
gene <- gene[gene$`Gene Symbol`%in% ndegs$Gene_ID,]
###########################
expr <- GSE18520.exprSet
group <- GSE18520.clin[,4,drop=F]

module.gene <- intersect(rownames(expr),gene$`Gene Symbol`)

expr <- as.data.frame(expr[module.gene,])
Group <- data.frame("id"=colnames(expr),
                    "Group" = ifelse(group$group=="Normal","Normal","Cancer"))

### 使用limma来做芯片的差异分析
### 1.创建分组
group <- Group$Group
### levels里面，把对照组放在前面
group <- factor(group,levels = c("Normal","Cancer"))
### 构建比较矩阵
design <- model.matrix(~group)
### 比较矩阵命名
colnames(design) <- levels(group)
design

### 2.线性模型拟合
fit <- lmFit(expr,design)
### 3.贝叶斯检验
fit2 <- eBayes(fit)
### 4.输出差异分析结果,其中coef的数目不能操过design的列数
### 此处的2代表的是design中第二列和第一列的比较
GSE18520.allDiff=topTable(fit2,adjust='fdr',coef=2,number=Inf)
###
#expr <- GSE26712.exprSet
#group <- GSE26712.clin[,1,drop=F]

#module.gene <- intersect(rownames(expr),gene$`Gene Symbol`)

#expr <- as.data.frame(expr[module.gene,])
#Group <- data.frame("id"=colnames(expr),
   #                 "Group" = ifelse(group$source_name_ch1=="Ovarian surface epithelial cells","Normal","Cancer"))

### 使用limma来做芯片的差异分析
### 1.创建分组
#group <- Group$Group
### levels里面，把对照组放在前面
#group <- factor(group,levels = c("Normal","Cancer"))
### 构建比较矩阵
#design <- model.matrix(~group)
### 比较矩阵命名
#colnames(design) <- levels(group)
#design

### 2.线性模型拟合
#fit <- lmFit(expr,design)
### 3.贝叶斯检验
#fit2 <- eBayes(fit)
### 4.输出差异分析结果,其中coef的数目不能操过design的列数
### 此处的2代表的是design中第二列和第一列的比较
GSE26712.allDiff=topTable(fit2,adjust='fdr',coef=2,number=Inf)

###
expr <- GSE40595.exprSet
group <- GSE40595.clin[,1,drop=F]

module.gene <- intersect(rownames(expr),gene$`Gene Symbol`)

expr <- as.data.frame(expr[module.gene,])
Group <- data.frame("id"=colnames(expr),
                    "Group" = ifelse(group$Group=="Human ovarian surface epthelium","Normal","Cancer"))

### 使用limma来做芯片的差异分析
### 1.创建分组
group <- Group$Group
### levels里面，把对照组放在前面
group <- factor(group,levels = c("Normal","Cancer"))
### 构建比较矩阵
design <- model.matrix(~group)
### 比较矩阵命名
colnames(design) <- levels(group)
design

### 2.线性模型拟合
fit <- lmFit(expr,design)
### 3.贝叶斯检验
fit2 <- eBayes(fit)
### 4.输出差异分析结果,其中coef的数目不能操过design的列数
### 此处的2代表的是design中第二列和第一列的比较
GSE40595.allDiff=topTable(fit2,adjust='fdr',coef=2,number=Inf)
##############################
GSE18520 <- GSE18520.allDiff %>%
            rownames_to_column("ID") %>%
            select(ID,logFC)
#GSE26712 <- GSE26712.allDiff %>%
#  rownames_to_column("ID") %>%
 # select(ID,logFC)
GSE40595 <- GSE40595.allDiff %>%
  rownames_to_column("ID") %>%
  select(ID,logFC)

logfc <- GSE18520 %>%
  # full_join(.,GSE26712,by ="ID") %>%
   full_join(.,GSE40595,by ="ID")%>%
   column_to_rownames("ID")
colnames(logfc) <- c("GSE18520",#"GSE26712",
                     "GSE40595")
#####################
GSE18520.p <- GSE18520.allDiff %>%
  rownames_to_column("ID") %>%
  select(ID,adj.P.Val)
#GSE26712.p <- GSE26712.allDiff %>%
#  rownames_to_column("ID") %>%
#  select(ID,adj.P.Val)
GSE40595.p <- GSE40595.allDiff %>%
  rownames_to_column("ID") %>%
  select(ID,adj.P.Val)

p.val <- GSE18520.p %>%
  #full_join(.,GSE26712.p,by ="ID") %>%
  full_join(.,GSE40595.p,by ="ID")%>%
  column_to_rownames("ID")
colnames(p.val) <- c("GSE18520",#"GSE26712",
                     "GSE40595")

genetype <- ndegs %>% column_to_rownames("Gene_ID")
genetype <- genetype[gene$`Gene Symbol`[71:1],]
genetype$type <- ifelse(genetype$Protective>genetype$Risky,"Protective","Risky")
########################按分类
logfcOrdered <- logfc[rownames(genetype),]
p.valOrdered <- p.val[rownames(genetype),]

## 对logfc进行分类
logfcCat <- apply(logfcOrdered, 2, function(x){
  cut(x, breaks = c(-Inf, -2, -1, 1, 2, Inf),
      labels = c("< -2", "-2 - -1", "-1 - 1", "1 - 2", "> 2"))
})
rownames(logfcCat) <- rownames(logfcOrdered)

## 确保两个数据集的列名和行名顺序一致
if(!identical(rownames(logfcCat), rownames(p.valOrdered))) p.valOrdered <- p.valOrdered[rownames(logfcCat),]
if(!identical(colnames(logfcCat), colnames(p.valOrdered))) p.valOrdered <- p.valOrdered[colnames(logfcCat)]

## 把P> 0.05的数据标记为 P > 0.05
logfcCat[p.valOrdered >= 0.05] <- "P >= 0.05"
#logfcCat[is.na(logfcCat)] <- "N/A"

## 查看有多少个分类
unique(matrix(logfcCat, ncol  = 1))

## 每个分类定义一个颜色
col_cat <- c("> 2" = "#A80C3A", "1 - 2" = "#ED5E57", "-1 - 1" = "#DDD3D2",
             "-2 - -1" = "#6B9AB7", "< -2" = "#2F5B89", "P >= 0.05" = "white")

cell_fun <- function(logfc, dataP, logfcCutoff = 1, PCutoff = 0.05, 
                     darkcol = "black", lightcol = "white", digit = 2, fontsize  = 6){
  function(j, i, x, y, width, height, fill){
    if(abs(logfc[i,j]) > logfcCutoff & dataP[i,j] < PCutoff){
      grid.text(round(logfc, digit)[i, j], x, y, 
                gp = gpar(fontsize = fontsize, col  = lightcol))
    }else{
      grid.text(round(logfc, digit)[i, j], x, y, 
                gp = gpar(fontsize = fontsize, col  = darkcol))
    }
  }
}

## 定义注释信息的颜色
an_col <- c("#0FCFC0","#C75DAA")
names(an_col) <- unique(genetype$type)
library(ComplexHeatmap)
## 定义注释信息
row_an <-  HeatmapAnnotation(type = genetype$type, ##注释信息的内容
                             show_annotation_name = F, ## 是否显示注释的标题
                             col = list(type = an_col), ## 注释信息的颜色
                             show_legend = T,  ## 是否显示注释信息的说明
                             annotation_legend_param = list(title = "HR"), ## 注释信息图例的标题
                             simple_anno_size = unit(0.25, "cm"),
                             which = "row") #对行或者列进行注释 


Heatmap(matrix = logfcCat, 
        name = "logFC", #主要图例的标题
        rect_gp = gpar(col = "NA", lwd = 1), #不画边框，或者用col = "grey"画灰色边框
        col = col_cat, #热图颜色
        row_names_side = "left", 
        cell_fun = cell_fun(logfcOrdered, p.valOrdered), 
        row_names_gp = gpar(fontsize = 8), #基因名字号
        column_names_gp = gpar(fontsize = 8), #肿瘤类型字号
        #row_names_rot = 90,
        column_names_rot = 90, #肿瘤类型呈45度
        left_annotation = row_an) #左侧基因分类，如果不画，就筛掉这个参数


library(export)
graph2ppt(file="out/GEO.exp.pptx")