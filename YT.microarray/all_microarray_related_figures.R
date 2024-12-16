library(tidyverse)
library(ggplot2)
library(limma)

# Data is deposited in ArrayExpress under accession number E-MTAB-14643
yt.data <- read.csv("./siYT_AT2.eset.csv")

yt.data <- yt.data  %>% column_to_rownames("X")
head(yt.data)

metadata <- yt.data %>%
    colnames() %>%
    as_tibble() %>%
    separate(value, into = c("disease", "gene", "experiment"), sep = "_",remove = FALSE, extra = "merge")  %>% mutate(condition=paste0(disease, "_", gene))


metadata$condition <- factor(metadata$condition)
yt.data <- yt.data[,metadata$value]

# Create design matrix
design <- model.matrix(~0 + condition, data = metadata)
colnames(design) <- levels(metadata$condition)
design

# Define contrasts
contrast.matrix <- makeContrasts(PBS_YT - PBS_SC, levels = design)
# Fit the linear model
fit <- lmFit(yt.data, design)
# Apply contrasts to the linear model
fit2 <- contrasts.fit(fit, contrast.matrix)
# Compute statistics
fit2 <- eBayes(fit2)
# Get the results
results <- topTable(fit2, adjust = "BH", number = Inf)

# View the results
results <- as.data.frame(results)  %>% arrange(P.Value) 
head(results, 20)
results  %>% write.csv(file="./pbs_siYT_vs_siSC_all.csv")




# Get the top 100 genes
top_genes <- results %>% top_n(-100, P.Value) %>% rownames_to_column("Gene") %>% pull(Gene)
# Subset the data for the top 100 genes
yt.data.top <- yt.data[top_genes, ]

yt.data.top %>% select(metadata[metadata$condition%in% c("PBS_YT","PBS_SC"),]$value) %>% as.data.frame() %>% rownames_to_column("Gene") %>% write.csv(file="./yt_pbs_top100.csv")

logFC_pbs100 <- yt.data.top %>% select(metadata[metadata$condition%in% c("PBS_YT","PBS_SC"),]$value) %>% 
                                as.data.frame()  %>% 
                                mutate(PBS_YT104 = log2(PBS_YT_S104/PBS_SC_S104), 
                                       PBS_YT105 = log2(PBS_YT_S105/PBS_SC_S105), 
                                       PBS_YT104_105 = log2(PBS_YT_S104_105/PBS_SC_S104_105)) %>% 
                                       select(c("PBS_YT104", "PBS_YT105", "PBS_YT104_105")) 
logFC_pbs100

## heatmap generation. Sample graph.
library(pheatmap)
pheatmap(logFC_pbs100, 
         color = colorRampPalette(c("purple", "grey", "yellow"))(100), 
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         cellwidth = 30,
         cellheight = 10,
         main = "Heatmap of logFC for top 100 genes")



# Perform Gene Ontology enrichment analysis
library(clusterProfiler)
library(org.Mm.eg.db)

top_genes1 <- results %>% filter(P.Value<0.001) %>% rownames_to_column("Gene") %>% pull(Gene)
ego <- enrichGO(gene = top_genes1,
                        OrgDb = org.Mm.eg.db,
                        keyType = "SYMBOL",
                        ont = "CC",
                        pAdjustMethod = "fdr",
                        qvalueCutoff = 0.05,
                        readable = TRUE)

# View the results
head(ego)

# Plot the results
#saveRDS(ego, file = "./ego.rds")
#ego <- readRDS("./ego.rds")
clusterProfiler::dotplot(ego, showCategory = 20)



# Perform the analysis for siYT vs siSC in Bleo AT2 cells.
contrast.matrix <- makeContrasts(Bleo_YT - Bleo_SC, levels = design)

# Fit the linear model
fit <- lmFit(yt.data, design)

# Apply contrasts to the linear model
fit2 <- contrasts.fit(fit, contrast.matrix)
# Compute statistics
fit2 <- eBayes(fit2)
# Get the results
results <- topTable(fit2, adjust = "fdr", number = Inf)

# View the results
results <- as.data.frame(results)  %>% arrange(P.Value) 
head(results, 20)

