
##---- This script was made for educational purposes ----##
##---- Data was taken after request from a published ----##
##---- Research Article N.Karagianni et al.(2019)    ----##
##---- https://doi.org/10.1371/journal.pcbi.1006933  ----##



##---- load libraries ----##
## quantiles normalization
library(preprocessCore)
## ploting
library(ggplot2)
## Dunnett post-hoc (optional package)
library(multcomp)
## heatmap.2
library(gplots)
## heatmap color clustering
library(RColorBrewer)
## kmeans clustering
library(factoextra)



##---- data preparation ----##
## set working directory
setwd("/home/touliopoulos/project/2pfsb/bioinformatics_workshop")

## load data
genes_data = read.delim("./Raw_common18704genes_antiTNF.tsv", header=T, row.names = 1, sep="\t")

## dimensions of dataframe
dim(genes_data)

## boxplot visualization
boxplot( genes_data, horizontal=T , las=1 , cex.axis=0.5 )

## keep gene and sample names
Gene = rownames(genes_data)
Sample = colnames(genes_data)

## normalization
genes_data = as.matrix(genes_data)
genes_data = normalize.quantiles(genes_data,copy=TRUE)
genes_data = data.frame(genes_data)

## apply gene and sample names to dataset after normalization
colnames(genes_data) = Sample
rownames(genes_data) = Gene

## remove variable
rm(Sample)

## boxplot visualization after normalization
boxplot( genes_data, horizontal=T , las=1 , cex.axis=0.5 )

## select number of rows for subsample
n = 18703
genes_data = head(genes_data , n)



##---- dimension reduction ----##
## keep only wt and tg samples
wt_tg_df = genes_data[, 1:23]

#after dataframe transposition columns must represent genes
wt_tg_df = t(wt_tg_df)

## create vector with groups
group = c(rep("A_Wt", 10), rep("B_Tg", 13) )

## PCA dimension reduction
wt_tg_df.pca = prcomp(wt_tg_df , scale. = FALSE)
summary(wt_tg_df.pca)
wt_tg_df.pca = data.frame("PC1" = wt_tg_df.pca$x[,1] , "PC2" = wt_tg_df.pca$x[,2] , "group" = group)

## plot PCA results
ggplot(wt_tg_df.pca , aes(x=PC1,y=PC2,color=group))+
  geom_point()+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank() )



##---- anova preparation ----##
## create matrix by excluding rownames and colnames
matrixdata = as.matrix(genes_data[1:n,])

## create groups
group = factor(c(rep("A_Wt", 10), rep("B_Tg", 13), rep("C_Proph_Ther_Rem", 3), rep("D_Ther_Rem", 10),
                rep("E_Ther_Hum", 10), rep("F_Ther_Enb", 10), rep("G_Ther_Cim", 10) ) )

## anova on first gene
gene1 = data.frame("gene_expression" = matrixdata[1,] , "group" = group)
geneaov = aov( gene_expression~group , data = gene1 )
summary(geneaov)

## calculate mean expression value for every group
group_mean_values = aggregate(gene1$gene_expression , list(gene1$group) , FUN=mean)
group_mean_values

## Tukeys post-hoc
TukeyHSD(geneaov , conf.level = 0.99)

## Access specific metrics from Tukeys post-hoc
tukey = TukeyHSD(geneaov)
tukey
tukey$group["B_Tg-A_Wt" ,  ]
tukey$group["B_Tg-A_Wt" , 1]
tukey$group["B_Tg-A_Wt" , 4]

tukey_data = c(tukey$group["B_Tg-A_Wt" , 1] , tukey$group["B_Tg-A_Wt" , 4] , 
               tukey$group["C_Proph_Ther_Rem-A_Wt" , 1] , tukey$group["C_Proph_Ther_Rem-A_Wt" , 4] ,
               tukey$group["D_Ther_Rem-A_Wt" , 1] , tukey$group["D_Ther_Rem-A_Wt" , 4] , 
               tukey$group["E_Ther_Hum-A_Wt" , 1] , tukey$group["E_Ther_Hum-A_Wt" , 4] , 
               tukey$group["F_Ther_Enb-A_Wt" , 1] , tukey$group["F_Ther_Enb-A_Wt" , 4] ,
               tukey$group["G_Ther_Cim-A_Wt" , 1] , tukey$group["G_Ther_Cim-A_Wt" , 4] )
tukey_data

## Dunnett's post-hoc
dunnett = glht(geneaov , linfct = mcp(group="Dunnett"))
summary(dunnett)
modgene = summary(dunnett)
modgene[[10]]$coefficients
modgene[[10]]$pvalues

## remove unwanted variables
rm(gene1,geneaov,group_mean_values,dunnett,modgene,tukey,tukey_data)

## Recursive anova on all genes
anova_table = data.frame()

for( i in 1:length(matrixdata[,1] ) ) 
{
  df = data.frame("gene_expression" = matrixdata[i,] , "group" = group)
  geneaov = aov( gene_expression~group , data = df )
  tukey = TukeyHSD(geneaov , conf.level = 0.99)
  
  tukey_data = c(tukey$group["B_Tg-A_Wt" , 1] , tukey$group["B_Tg-A_Wt" , 4] , 
                 tukey$group["C_Proph_Ther_Rem-A_Wt" , 1] , tukey$group["C_Proph_Ther_Rem-A_Wt" , 4] ,
                 tukey$group["D_Ther_Rem-A_Wt" , 1] , tukey$group["D_Ther_Rem-A_Wt" , 4] , 
                 tukey$group["E_Ther_Hum-A_Wt" , 1] , tukey$group["E_Ther_Hum-A_Wt" , 4] , 
                 tukey$group["F_Ther_Enb-A_Wt" , 1] , tukey$group["F_Ther_Enb-A_Wt" , 4] ,
                 tukey$group["G_Ther_Cim-A_Wt" , 1] , tukey$group["G_Ther_Cim-A_Wt" , 4] ,
                 
                 tukey$group["C_Proph_Ther_Rem-B_Tg" , 1] , tukey$group["C_Proph_Ther_Rem-B_Tg" , 4] , 
                 tukey$group["D_Ther_Rem-B_Tg" , 1] , tukey$group["D_Ther_Rem-B_Tg" , 4] ,
                 tukey$group["E_Ther_Hum-B_Tg" , 1] , tukey$group["E_Ther_Hum-B_Tg" , 4] , 
                 tukey$group["F_Ther_Enb-B_Tg" , 1] , tukey$group["F_Ther_Enb-B_Tg" , 4] , 
                 tukey$group["G_Ther_Cim-B_Tg" , 1] , tukey$group["G_Ther_Cim-B_Tg" , 4] )
  
  anova_table = rbind( anova_table , tukey_data )
}

## apply column names to anova table
colnames(anova_table) = c("Wt_Tg_diff" , "Wt_Tg_padj" ,
                          "Wt_Rem_P_diff" , "Wt_Rem_P_padj" , 
                          "Wt_Rem_diff" , "Wt_Rem_padj" , 
                          "Wt_Hum_diff" , "Wt_Hum_padj" , 
                          "Wt_Enb_diff" , "Wt_Enb_padj" , 
                          "Wt_Cim_diff" , "Wt_Cim_padj" ,
                          
                          "Tg_Rem_P_diff" , "Tg_Rem_P_padj" , 
                          "Tg_Rem_diff" , "Tg_Rem_padj" , 
                          "Tg_Hum_diff" , "Tg_Hum_padj" , 
                          "Tg_Enb_diff" , "Tg_Enb_padj" , 
                          "Tg_Cim_diff" , "Tg_Cim_padj")

## add column with gene names 
rownames(anova_table) = Gene

## remove unwanted variables
rm(df,geneaov,matrixdata,tukey,tukey_data,i)



##---- volcano plot dataframe preparation for DEGs from wt vs tg ----##
upWT = 0
downWT = 0
nochangeWT = 0

## filter genes based on mean diff and pvalue between wt and tg
upWT = which(anova_table[,1] < -1.0 & anova_table[,2] < 0.05)
downWT = which(anova_table[,1] > 1.0 & anova_table[,2] < 0.05)
nochangeWT = which(anova_table[,2] > 0.05 | (anova_table[,1] > -1.0 & anova_table[,1] < 1.0) )

## create vector to store states for each gene
state = vector(mode="character" , length=length(anova_table[,1]))
state[upWT] = "up_WT"
state[downWT] = "down_WT"
state[nochangeWT] = "nochange_WT"

## identify names of DEGs between wt and tg
genes_up_WT = c(rownames(anova_table)[upWT] )
genes_down_WT = c(rownames(anova_table)[downWT] )

## union of DEGs between wt and tg to acquire total Differential Expressed Genes
deg_wt_tg = c( genes_up_WT , genes_down_WT )

## subset dataframe based on specific degs
deg_wt_tg_df = subset( genes_data , Gene %in% deg_wt_tg )

## dataframe for volcano plot
volcano_data = data.frame( "padj" = anova_table[,2] , "DisWt" = anova_table[,1] , state=state )

## plot pvalues with mean diffs
ggplot( volcano_data , aes(x=DisWt , y=-log10(padj), colour=state))+
  geom_point()+
  geom_vline(xintercept = c(-1,1), color = 'red') + 
  geom_hline(yintercept = -log10(0.05), color = 'red')



##---- volcano plot dataframe preparation for DEGs from tg vs therapies ----##
upTHER = 0
downTHER = 0
nochangeTHER = 0

## filter genes based on mean diff and pvalue between th and therapies
upTHER = which( (anova_table[,13] < -1.0 & anova_table[,14] < 0.05) | 
                (anova_table[,15] < -1.0 & anova_table[,16] < 0.05) | 
                (anova_table[,17] < -1.0 & anova_table[,18] < 0.05) |
                (anova_table[,19] < -1.0 & anova_table[,20] < 0.05) |
                (anova_table[,21] < -1.0 & anova_table[,22] < 0.05) )

downTHER = which( (anova_table[,13] > 1.0 & anova_table[,14] < 0.05) | 
                  (anova_table[,15] > 1.0 & anova_table[,16] < 0.05) | 
                  (anova_table[,17] > 1.0 & anova_table[,18] < 0.05) |
                  (anova_table[,19] > 1.0 & anova_table[,20] < 0.05) |
                  (anova_table[,21] > 1.0 & anova_table[,22] < 0.05) )

nochangeTHER = which( ( (anova_table[,13] > -1.0 & anova_table[,13] < 1.0) | anova_table[,14] > 0.05) |
                  ( (anova_table[,15] > -1.0 & anova_table[,15] < 1.0) | anova_table[,16] > 0.05) |
                  ( (anova_table[,17] > -1.0 & anova_table[,17] < 1.0) | anova_table[,18] > 0.05) |
                  ( (anova_table[,19] > -1.0 & anova_table[,19] < 1.0) | anova_table[,20] > 0.05) |
                  ( (anova_table[,21] > -1.0 & anova_table[,21] < 1.0) | anova_table[,22] > 0.05) )

## create vector to store states for each gene
state = vector(mode="character" , length=length(anova_table[,1]))
state[upTHER] = "up_THER"
state[downTHER] = "down_THER"
state[nochangeTHER] = "nochange_THER"

## identify names of genes differentially expressed between tg and therapies
genes_up_THER = c(rownames(anova_table)[upTHER] )
genes_down_THER = c(rownames(anova_table)[downTHER] )

deg_tg_ther = c( genes_up_THER , genes_down_THER )

## subset dataframe based on these degs
deg_tg_ther_df = subset( genes_data , Gene %in% deg_tg_ther )

## remove unwanted variables
rm(downTHER,downWT,upTHER,upWT,state,nochangeTHER,nochangeWT,
   genes_down_THER,genes_up_THER,genes_down_WT,genes_up_WT)

## export subset dataframes to files
write.table( deg_wt_tg_df , file="./deg_wt_tg.tsv" , sep="\t" , quote=FALSE)
write.table( deg_tg_ther_df , file="./deg_tg_ther.tsv" , sep="\t" , quote=FALSE)



##---- PCA preparation after reducing dataset dimensions with DEGs ----##
## genes must be represented in columns
deg_wt_tg_df = t(deg_wt_tg_df)

## keep only observations from wt and tg
deg_wt_tg_df = deg_wt_tg_df[1:23,]

## group wt and tg as character and not factor
group = c(rep("A_Wt", 10), rep("B_Tg", 13) )

## dimension reduction with PCA for wt and tg dataframe
deg_wt_tg_df.pca = prcomp(deg_wt_tg_df , scale. = FALSE)
summary(deg_wt_tg_df.pca)
deg_wt_tg_df.pca = data.frame("PC1" = deg_wt_tg_df.pca$x[,1] , "PC2" = deg_wt_tg_df.pca$x[,2] , "group" = group)

## plot PCA results
ggplot(deg_wt_tg_df.pca , aes(x=PC1,y=PC2,color=group))+
  geom_point()+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank() )



## combine DEGs from tg and ther
DEGs = c( deg_tg_ther , deg_wt_tg )

## Data frame with all DEGs for clustering
DEGsFrame = anova_table[rownames(anova_table) %in% DEGs , ]
DEGsFrame = as.matrix(DEGsFrame)



##---- heatmap and hierarchical clustering ----##
hclustfunc = function(x) hclust(x, method="complete")
distfunc = function(x) dist(x, method="euclidean")

## perform clustering on rows and columns
cl.row = hclustfunc(distfunc(DEGsFrame[, c(1,3,5,7,9,11) ]))

## extract cluster assignments; i.e. k=8 (rows) k=5 (columns)
gr.row = cutree(cl.row, 5)

## apply a set of color pallete
col1 = brewer.pal(5, "Set3")

heatmap.2(DEGsFrame[, c(1,3,5,7,9,11) ], col = bluered(100), trace = "none", tracecol="black",
          density.info = "none", labCol = c("TG", "REM_P", "REM", "HUM", "ENB","CIM"),
          scale="none" , labRow="" , vline=0 , mar=c(6,2),
          RowSideColors = col1[gr.row])



##---- kmeans clustering ----##
kmeans = kmeans(DEGsFrame[ , c(1,3,5,7,9,11) ] , centers=6)
fviz_cluster(kmeans, data = (DEGsFrame[ , c(1,3,5,7,9,11) ] ), geom="point", show.clust.cent=TRUE)

## Extract genes from clusters
clusters = data.frame(kmeans$cluster)
colnames(clusters) = ("ClusterNo")

## extract specific cluster
cluster1 = data.frame(rownames(subset(clusters, ClusterNo==1)))
colnames(cluster1) = ("Gene")

