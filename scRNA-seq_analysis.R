## load libraries
library(dplyr)
library(Seurat)
library(patchwork)
library(biomaRt)
library(calibrate)
library(cowplot)
library(ggplot2)
library(uwot)
library(Matrix)
library(DropletUtils)
### Note: Performed with Serurat4.4.0

########## scRNA-seq_create_UMAP ############################################
#############################################################################
## Load 10xdata and metadata
r1.data <- Read10X(data.dir = "filtered_feature_bc_matrix")
r1.metadata <- read.csv(file = "meta.csv", row.names = 1)

## Create Seurat object
r1 <- CreateSeuratObject(r1.data, project = "r1", meta.data = r1.metadata, min.features =  100)
r1
View(r1@meta.data)

## Add number of genes per UMI for each cell to metadata
r1$log10GenesPerUMI <- log10(r1$nFeature_RNA) / log10(r1$nCount_RNA)

## Compute percent mito ratio
r1$mitoRatio <- PercentageFeatureSet(r1, pattern = "^mt-")
r1$mitoRatio <- r1@meta.data$mitoRatio / 100

## Create metadata dataframe
metadata <- r1@meta.data
metadata <- metadata %>% dplyr::rename(seq.folder = orig.ident, nUMI = nCount_RNA, nGene = nFeature_RNA)
r1@meta.data <- metadata

## Visualize the number of cell counts per sample
metadata %>% 
  ggplot(aes(x = stage, fill = stage)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  ggtitle("nCells_r1")

## Visualize the number UMIs counts per cell via histogram
metadata %>% 
  ggplot(aes(color = stage, x = nUMI, fill = stage)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 7300)+
  ggtitle("nUMIs_per_nCells")

## Visualize the distribution of genes detected per cell via histogram
metadata %>% 
  ggplot(aes(color = stage, x = nGene, fill = stage)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() + 
  geom_vline(xintercept = 2300) +
  ggtitle("nCells_vs_nGenes")

## Visualize the correlation between genes detected and number of UMIs and determine whether strong presence of cells with low numbers of genes/UMIs
metadata %>% 
  ggplot(aes(x = nUMI, y = nGene, color = mitoRatio)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method = lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 7300) +
  geom_hline(yintercept = 2300) +
  facet_wrap(~stage) +
  ggtitle("correlation_of_nGene_and_nUMI")

## Visualize the distribution of mitochondrial gene expression detected per cell via histogram
metadata %>% 
  ggplot(aes(color = stage, x = mitoRatio, fill = stage)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  geom_vline(xintercept = 0.1) +
  ggtitle("mitoRatio")

## Visualize the overall complexity of the gene expression by visualizing the genes detected per UMI via histogram
metadata %>%
  ggplot(aes(x = log10GenesPerUMI, color = stage, fill = stage)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 0.8) +
  ggtitle("genes_detected_per_UMI")


## Filter out low quality reads using selected thresholds - these will change with experiment
r1 <- subset(x = r1, 
             subset = (nUMI >= 7300) & 
               (nGene >= 2300) & 
               (log10GenesPerUMI > 0.8) & 
               (mitoRatio < 0.1))
r1

## Visualize the number of cell counts per sample after filtration
r1@meta.data %>% 
  ggplot(aes(x = stage, fill = stage)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  ggtitle("nCells_filtered")

VlnPlot(r1, features = c("nUMI", "nGene", "mitoRatio"), ncol = 3, pt.size = 0.2)

## Normalizing data (Log-normalization)
r1 <- NormalizeData(r1, normalization.method = "LogNormalize", scale.factor = 10000)
r1 <- NormalizeData(r1)

## Scaling data
r1 <- ScaleData(r1, vars.to.regress = c("mitoRatio", "CC.Difference"), features = rownames(r1), verbose = TRUE)

## PCA
r1 <- RunPCA(r1, features = rownames(r1), npcs =100, verbose = FALSE)

## Examining and visualizing PCA results a few different ways
print(x = r1[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(object = r1, dims = 1:2, reduction = "pca")
DimPlot(object = r1, reduction = "pca")

## Determining the 'dimensionality' of the dataset
ElbowPlot(object = r1, ndims = 40)

## Cluster the cells
r1.dim40 <- FindNeighbors(object = r1, reduction = "pca", dims = 1:40, nn.eps = 0)

## Determine resolution
r1.dim40.res0.01<- FindClusters(r1.dim40, resolution = 0.01, n.start = 10, algorithm = 3)

## Visualize UMAP
r1.dim40.res0.01.UMAP <- RunUMAP(object = r1.dim40.res0.01, dims = 1:40)
DimPlot(object = r1.dim40.res0.01.UMAP, reduction = "umap", label = TRUE, pt.size = 0.2, group.by = "seurat_clusters", label.size = 5) + ggtitle("UMAP")

## Export coordination data of UMAP
write.csv(r1.dim40.res0.01.UMAP@reductions$umap@cell.embeddings, file = paste(" ", ".csv", sep=""))
write.csv(r1.dim40.res0.01.UMAP@active.ident, file = paste(" ", ".csv", sep=""))
#### These files were imported into the loupe browser to view gene expression.



########## Bmpr1a+ AEC vs Bmpr1a-AEC #########################################
#############################################################################

## Extraction of AEC cluster
subset.r1 <- subset(r1.dim40.res0.01.UMAP, idents = c(0, 1, 2, 4), invert = TRUE)
subset.r1_Dll4<- subset(x = subset.r1, subset = Dll4 > 0.1)
subset.r1_Dll4

## Bmpr1a+ AEC vs Bmpr1a- AEC
Idents(subset.r1_Dll4, WhichCells(object = subset.r1_Dll4, expression = Bmpr1a > 0.23, slot = 'data')) <- 'Bmpr1a.pos'
Idents(subset.r1_Dll4, WhichCells(object = subset.r1_Dll4, expression = Bmpr1a <= 0.23, slot = 'data')) <- 'Bmpr1a.neg'
Bmpr1a.genes <- FindMarkers(subset.r1_Dll4, ident.1 = 'Bmpr1a.pos', ident.2 = 'Bmpr1a.neg',min.pct = 0.25,logfc.threshold = 0)
Bmpr1a.genes
write.csv(x = head(x = Bmpr1a.genes), file = ".csv")

## Plotting volcano
with(Bmpr1aCdh5Dll4.Pos.Neg, plot(avg_log2FC, -log10(p_val_adj), pch = 1, cex = 1.12, main = "Volcano_plot_Bmpr1a.Pos_Neg",xlim = c(-1.5, 1.5), ylim = c(0, 3.5) ))
abline(v = 0.58496, col= 'black', lty=2) 
abline(v = -0.58496, col= 'black', lty=2)
abline(h = 1.301, col ='black', lty=2 )
with(subset(Bmpr1aCdh5Dll4.Pos.Neg, (avg_log2FC) < -0.58496 & (-log10(p_val_adj)) > 1.301), points(avg_log2FC, -log10(p_val_adj), pch = 20, cex = 1.5, col = "blue"))
with(subset(Bmpr1aCdh5Dll4.Pos.Neg, (avg_log2FC) > 0.58496 & (-log10(p_val_adj)) > 1.301) ,points(avg_log2FC, -log10(p_val_adj), pch = 20, cex = 1.5, col = "Red"))
with(subset(Bmpr1aCdh5Dll4.Pos.Neg, (avg_log2FC) > 0.58496 & (-log10(p_val_adj)) > 1.301) , textxy(avg_log2FC, -log10(p_val_adj), labs = genes , cex = 0.5))



########## Ligand_Exp_in Bmpr1a-EC and Receptor_Exp_in Bmpr1a+ HEC ##########
#############################################################################

## Extraction of Cdh5+ EC cluster
subset.r1 <- subset(r1.dim40.res0.01.UMAP, idents = c(0, 1, 2, 4), invert = TRUE)
subset.r1_Cdh5<- subset(x = subset.r1, subset = Cdh5 > 0.29)
subset.r1_Cdh5

## Definition of Bmpr1a- EC and Bmpr1a+Runx1+ HEC
Idents(r1.dim40.res0.01.UMAP, WhichCells(object = subset.r1_Cdh5, expression = Bmpr1a <= 0.23, slot = 'data')) <- 'Bmpr1a.neg'
Idents(r1.dim40.res0.01.UMAP, WhichCells(object = subset.r1_Cdh5, expression = Bmpr1a > 0.23 & Runx1 > 0.12, slot = 'data')) <- 'Bmpr1a.posRunx1.pos'

## Geneset import : Get the receptor-ligand list from Dimitrov Nature communication 2022
## Extraction of average expression levels, expression ratios, and normalized expression levels of ligand genes
features1<- (file= "Ligand_list.csv")
g1 <- DotPlot(object = r1.dim40.res0.01.UMAP, features = features1, assay="RNA")
View(g1$data)
write.csv(g1$data, file = "Ligand.csv")
#### This csv file was used to perform screening of ligand-receptor sets based on the ave. Exp. and pct. Exp. value of the ligand expression in Bmpr1a- EC.

## Extraction of average expression levels, expression ratios, and normalized expression levels of receptor genes
features2<- (file= "Screened_receptor_list.csv")
g2 <- DotPlot(object = r1.dim40.res0.01.UMAP, features = features2, assay="RNA")
View(g1$data)
write.csv(g2$data, file = "Receptor.csv")
#### This csv file was used to perform screening of ligand-receptor sets based on the ave. Exp. and pct. Exp. value of the receptor expression in Bmpr1a+Runx1+ HEC.
















