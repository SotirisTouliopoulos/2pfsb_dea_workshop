## Data & Code for a Differential Expression Analysis workshop

Data was taken after request from a published Research Article: N.Karagianni et al.(2019). DOI: https://doi.org/10.1371/journal.pcbi.1006933

### First we have to load the required libraries for analysis
```
library(preprocessCore)
library(ggplot2)
library(multcomp)
library(gplots)
library(RColorBrewer)
library(factoextra)
```

### Loading od data with "read.delim" function
```
genes_data = read.delim("./Raw_common18704genes_antiTNF.tsv",
                        header=T,
                        row.names = 1,
                        sep="\t")
```

### Get the dimensions of the loaded dataframe
```
dim(genes_data)
```

### Visualize gene expression distributions with boxplots
```
boxplot(genes_data,
        horizontal=T,
        las=1,
        cex.axis=0.5 )
```

### Save gene and sample names in a vector
```
Gene = rownames(genes_data)
Sample = colnames(genes_data)
```

### Normalize data to make distributions comparable
```
genes_data = as.matrix(genes_data)
genes_data = normalize.quantiles(genes_data,copy=TRUE)
genes_data = data.frame(genes_data)
```

### Set the gene and sample names to the normalised dataframe
```
colnames(genes_data) = Sample
rownames(genes_data) = Gene
```

### Boxplot visualization after normalization
```
boxplot( genes_data, horizontal=T , las=1 , cex.axis=0.5 )
```


