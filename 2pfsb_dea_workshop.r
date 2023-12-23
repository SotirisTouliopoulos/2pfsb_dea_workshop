
## -------
## Differential Expression Analysis Workshop
## Part 2
## -------

## set working directory
setwd("~")

## load normalized data
genes_data <- read.delim("~/Raw_common18704genes_antiTNF_normalized.tsv", header=T, sep="\t")

## boxplot visualization
boxplot( genes_data[,2:67] , horizontal=T , las=1 , cex.axis=0.5 )

## option to subsample for CPU/memory issues
n = 18703
genes_data = head(genes_data , n)
gene_names = genes_data[ , 1]

## create matrix by excluding rownames and colnames
matrixdata = as.matrix(genes_data[1:n , 2:67])

## create groups
group = c(rep("A_Wt", 10), rep("B_Tg", 13), rep("C_Proph_Ther_Rem", 3), rep("D_Ther_Rem", 10),
                rep("E_Ther_Hum", 10), rep("F_Ther_Enb", 10), rep("G_Ther_Cim", 10))

## anova on first gene
gene1 = data.frame("gene_expression" = matrixdata[1,] , "group" = group)
geneaov = aov( gene_expression~group , data = gene1 )
summary(geneaov)

## Tukeys post-hoc
TukeyHSD(geneaov , conf.level = 0.95)

## Access specific metrics from Tukeys post-hoc
tukey = TukeyHSD(geneaov)
tukey$group["B_Tg-A_Wt" ,  ]
tukey$group["B_Tg-A_Wt" , 1]
tukey$group["B_Tg-A_Wt" , 4]

## Access needed information all by once and store them in a vector
tukey_data = c(tukey$group["B_Tg-A_Wt" , 1] , tukey$group["B_Tg-A_Wt" , 4] , 
               tukey$group["C_Proph_Ther_Rem-A_Wt" , 1] , tukey$group["C_Proph_Ther_Rem-A_Wt" , 4] ,
               tukey$group["D_Ther_Rem-A_Wt" , 1] , tukey$group["D_Ther_Rem-A_Wt" , 4] , 
               tukey$group["E_Ther_Hum-A_Wt" , 1] , tukey$group["E_Ther_Hum-A_Wt" , 4] , 
               tukey$group["F_Ther_Enb-A_Wt" , 1] , tukey$group["F_Ther_Enb-A_Wt" , 4] ,
               tukey$group["G_Ther_Cim-A_Wt" , 1] , tukey$group["G_Ther_Cim-A_Wt" , 4] )
tukey_data


## Recursive anova on all genes

## create an empty dataframe
anova_table = data.frame()

## for every row (every gene)
for( i in 1:length(matrixdata[,1] ) ) 
{
  ## create dataframe from each row with two columns: expression and group 
  df = data.frame("gene_expression" = matrixdata[i,] , "group" = group)
  ## apply anova on the dataframe
  geneaov = aov( gene_expression~group , data = df )
  ## apply tukeys post-hoc test to anova results
  tukey = TukeyHSD(geneaov)
  ## select specific pairwise comparisons information (WT as control & diff,padj columns)
  tukey_data = c(tukey$group["B_Tg-A_Wt" , 1] , tukey$group["B_Tg-A_Wt" , 4] , 
                 tukey$group["C_Proph_Ther_Rem-A_Wt" , 1] , tukey$group["C_Proph_Ther_Rem-A_Wt" , 4] ,
                 tukey$group["D_Ther_Rem-A_Wt" , 1] , tukey$group["D_Ther_Rem-A_Wt" , 4] , 
                 tukey$group["E_Ther_Hum-A_Wt" , 1] , tukey$group["E_Ther_Hum-A_Wt" , 4] , 
                 tukey$group["F_Ther_Enb-A_Wt" , 1] , tukey$group["F_Ther_Enb-A_Wt" , 4] ,
                 tukey$group["G_Ther_Cim-A_Wt" , 1] , tukey$group["G_Ther_Cim-A_Wt" , 4] )
  ## append these 12 values as columns to dataframe
  anova_table = rbind( anova_table , tukey_data)
}

## change column names
colnames(anova_table) = c("Tg_diff" , "Tg_padj" ,
                          "Proph_Rem_diff" , "Proph_Rem_padj" , 
                          "Rem_diff" , "Rem_padj" , 
                          "Hum_diff" , "Hum_padj" , 
                          "Enb_diff" , "Enb_padj" , 
                          "Cim_diff" , "Cim_padj")




## keep only significant differences (p < 0.05) 
## categorize them in up / down regulated or no changed

## negative mean diff with WT as control ==> WT mean value higher
upWT = which(anova_table[,1] < 0 & anova_table[,2] < 0.05)
downWT = which(anova_table[,1] > 0 & anova_table[,2] < 0.05)
nochangeWT = which(anova_table[,2] > 0.05 )

## create vector to store changes information
state = vector(mode="character" , length=length(anova_table[,1]))
state[upWT] = "up"
state[downWT] = "down"
state[nochangeWT] = "nochange"

volcano_data = data.frame( anova_table[,1:2] , state=state )
colnames(volcano_data) = c("Tg_diff" , "Tg_padj" , "state")

library(ggplot2)
## volcano plot
ggplot( volcano_data , aes(x=Tg_diff , y=-log10(Tg_padj), colour=state))+
  geom_point()


## transpose dataframe if needed
anova_table = t(anova_table)

## add gene names to first column
colnames(anova_table) = gene_names

# Filtering DEGs with adjusted p-value <= 0.05 and log(FoldChange) >= 1
DegTG <- anova_table[1:2, which(abs(anova_table["Tg_diff",])>=0 &
                            anova_table["Tg_padj",] <= 0.05)]
DegRemP <- anova_table[3:4, which(abs(anova_table["Proph_Rem_diff",])>=0 &
                             anova_table["Proph_Rem_padj",]<=0.05)]
DegREM <- anova_table[5:6, which(abs(anova_table["Rem_diff",])>=0 &
                             anova_table["Rem_padj",]<=0.05)]
DegHUM <- anova_table[7:8, which(abs(anova_table["Hum_diff",])>=0 &
                             anova_table["Hum_padj",]<=0.05)]
DegENB <- anova_table[9:10, which(abs(anova_table["Enb_diff",])>=0 &
                              anova_table["Enb_padj",]<=0.05)]
DegCIM <- anova_table[9:10, which(abs(anova_table["Cim_diff",])>=0 &
                              anova_table["Cim_padj",]<=0.05)]

# Union of all the DEGs
DEGs <-data.frame(union(union(union(union(union(colnames(DegTG),
                                                colnames(DegRemP)),
                                                colnames(DegREM)),
                                                colnames(DegHUM)),
                                                colnames(DegENB)),
                                                colnames(DegCIM)))

colnames(DEGs) <- ("Genes")

# Data frame with all DEGs
DEGsFrame <- subset(anova_table, select=DEGs$Genes)


## heatmap and hierarchical clustering
library(gplots)
heatmap.2(t(DEGsFrame[c(1,3,5,7,9,11),]), col = bluered(100), trace = "none",
          density.info = "none", labCol = c("TG", "REM_P", "REM", "HUM", "ENB","CIM"),
          scale="none" , labRow="" , vline=0 , mar=c(6,2))


# kmeans clustering
library(factoextra)
kmeans <- kmeans(t(DEGsFrame[c(1,3,5,7,9,11),]), centers=6)
fviz_cluster(kmeans, data = t(DEGsFrame[c(1,3,5,7,9,11),]), geom="point", show.clust.cent=TRUE)

# Extract genes from clusters
clusters <- data.frame(kmeans$cluster)
colnames(clusters) <- ("ClusterNo")

cluster1 <- data.frame(rownames(subset(clusters, ClusterNo==1)))
colnames(cluster1) <- ("Gene")

cluster2 <- data.frame(rownames(subset(clusters, ClusterNo==2)))
colnames(cluster2) <- ("Gene")

cluster3 <- data.frame(rownames(subset(clusters, ClusterNo==3)))
colnames(cluster3) <- ("Gene")

cluster4 <- data.frame(rownames(subset(clusters, ClusterNo==4)))
colnames(cluster4) <- ("Gene")

cluster5 <- data.frame(rownames(subset(clusters, ClusterNo==5)))
colnames(cluster5) <- ("Gene")

cluster6 <- data.frame(rownames(subset(clusters, ClusterNo==6)))
colnames(cluster6) <- ("Gene")