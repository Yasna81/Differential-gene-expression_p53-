counts <- read.delim("https://www.ebi.ac.uk/gxa/experiments-content/E-MTAB-10041/resources/DifferentialSecondaryDataFiles.RnaSeq/raw-counts")
metadata <- read.delim("https://www.ebi.ac.uk/gxa/experiments-content/E-MTAB-10041/resources/ExperimentDesignFile.RnaSeq/experiment-design")
head(counts)
#counts
row.names(counts) <- counts$Gene.ID
genes <- counts[,c(1,2)]
counts <- counts[,-c(1,2)]
#metadata 
rownames(metadata) <- metadata$Run
metadata <- metadata [,c("Factor.Value.phenotype.","Factor.Value.compound."), drop= FALSE]
colnames(metadata) <- c("phenotype", "compound")
metadata$phenotype[metadata$phenotype=="p53R248Q heterozygous mutant"] <- "gainoffunction"
metadata$phenotype[metadata$phenotype=="p53 heterozygous mutant"] <- "lossoffunction"
metadata$phenotype[metadata$phenotype=="wild type phenotype"] <- "normal"
metadata$compound[metadata$compound == "none"] <- "control"
metadata$compound[metadata$compound == "Nutlin-3a 150 milligram per kilogram per day"] <- "Treated"
#turning to factor
head(metadata)
metadata$phenotype <- factor(metadata$phenotype,levels = c("normal","lossoffunction","gainoffunction"))
metadata$phenotype
metadata$compound <- factor(metadata$compound,levels = c("control","Treated"))
metadata$compound
#Run DESQ2
library(DESeq2)
dim(metadata)
dim(counts)
#metadata had 16 rows while counts had 12 rows so we edit that : 
samples_to_keep <- colnames(counts)
metadata_filtered <- metadata[rownames(metadata) %in% samples_to_keep, ]
dim(metadata_filtered)

dds <- DESeqDataSetFromMatrix(countData=counts, colData=metadata_filtered, design=~phenotype + compound)
dds <- dds[rowSums(counts(dds)) > 10, ]
dds <- DESeq(dds)
res = results(dds, contrast=c("compound", "Treated", "control"), alpha=1e-5)
levels(dds$phenotype)
table(dds$phenotype)
dds$phenotype <- relevel(dds$phenotype, ref ="lossoffunction")
res <-results(dds)
res = results(dds,contrast = c("phenotype","gainoffunction","lossoffunction"), alpha = 1e-5)
#normal should be first and refrence solve it. 
resultsNames(dds)
plotMA(res)
plotDispEsts(dds)
#merging genes
res_df = as.data.frame(res)
head(res_df)
head(genes)
res_df = merge(res_df, genes, by='row.names')
head(res_df)

library(EnhancedVolcano)
EnhancedVolcano(res, lab=rownames(res), x='log2FoldChange', y='pvalue')
genes_to_check = c("Gm3373", "Gm6937","Erdr1")
res_df[res_df$Gene.Name %in% genes_to_check, ]
#heatmap 
# get normalized counts
rld <- rlog(dds, blind = FALSE)

# select top genes 
topVarGenes <- head(order(rowVars(assay(rld)), decreasing = TRUE), 30)

library(pheatmap)
colnames(genes) <- c("Gene.ID", "Gene.Name")


vsd <- vst(dds, blind = TRUE)

topVarGenes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 30)
top_mat <- assay(vsd)[topVarGenes, ]

# Replace Ensembl IDs with gene names
gene_names <- genes$Gene.Name[match(rownames(top_mat), genes$Gene.ID)]

# keep Ensembl ID if no gene name is found
gene_names[is.na(gene_names)] <- rownames(top_mat)[is.na(gene_names)]
rownames(top_mat) <- gene_names


pheatmap(top_mat,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         annotation_col = as.data.frame(colData(dds)[, c("phenotype", "compound")]))
