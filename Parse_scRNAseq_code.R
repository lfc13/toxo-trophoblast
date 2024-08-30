#Parse single cell analysis with Seurat

#loading necessary packages
library(igraph)
library(dplyr)
library(patchwork)
library(DESeq2)
library(EnhancedVolcano)
library(Seurat)
library(Matrix)
library(ape)
library(ggplot2)
library(monocle3)
library(SeuratWrappers)
library(ggplot2)
library(patchwork)
library(magrittr)
library(reticulate)
library(ggforce)
library(ggpubr)
library(BiocManager)
library(gridExtra)
library(cluster)
library(factoextra)

rm(list = ls()) #clear environment

#convenience functions for saving objects and plots
SaveFigure <- function(plots, name, type = "png", width, height, res){
  if(type == "png") {png(paste0(fig_path, name, ".", type),
                         width = width, height = height, units = "in", res = 300)
  } else {
    pdf(paste0(fig_path, name, ".", type),
        width = width, height = height)
  }
  print(plots)
  dev.off()
}

SaveObject <- function(object, name){
  saveRDS(object, paste0(data_path, name, ".RDS"))
}

ReadObject <- function(name){
  readRDS(paste0(data_path, name, ".RDS"))
}

#where objects and figures are saved
data_path <- "path"
fig_path <- "path"

#color palette
preEVT <- "#F564E3"
CTB <- "#39BA38"
preSTB1 <- "#619CFF"
preSTB2 <- "#B79F01"
EVT <- "#3FBFC4"
STB <- "#F6766D"

#_______________________________________________________________________________

#Reading in data (requires Seurat v > 4.0.8)
mat_path <- "path_to_DGE_filtered"
mat <- ReadParseBio(mat_path)
mat

#Check to see if empty gene names are present, add name if so.
table(rownames(mat) == "")
rownames(mat)[rownames(mat) == ""] <- "unknown"

#read in cell meta data
cell_meta <- read.csv(paste0(mat_path, "/cell_metadata.csv"), row.names = 1)

#create object
trophall <- CreateSeuratObject(mat, min_genes = 100, min_cells = 5, names.field = 1, meta.data = cell_meta)

#setting our initial cell class to a single type, this will change after clustering
trophall@meta.data$orig.ident <- factor(rep("trophall", nrow(trophall@meta.data)))
Idents(trophall) <- trophall@meta.data$orig.ident
Idents(trophall)

SaveObject(trophall, "trophall_before_QC")
trophall <- ReadObject("trophall_before_QC")

#_______________________________________________________________________________

#cell quality control to prevent outliers from influence downstream analyses
trophall[["percent.mt"]] <- PercentageFeatureSet(trophall, pattern = "^MT-")
plot <- VlnPlot(trophall, pt.size = 0.10, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), group.by = "sample", ncol = 3)
SaveFigure(plot, "vln_QC_beforerfilter_bysample", width = 12, height = 6)
plot

#FeatureScatter used to visualize feature-feature relationships
plot1 <- FeatureScatter(trophall, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = "sample")
plot2 <- FeatureScatter(trophall, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "sample")
SaveFigure((plot1 + plot2), "scatter_QC_beforefilter_bysample", width = 12, height = 6, res = 200)

#Reassigning object based off of QC filtering cutoffs
trophall #number of cells before filtering
plot <- VlnPlot(trophall, pt.size = 0.10, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), group.by = "sample", ncol = 3)

trophall <- subset(trophall, subset = nFeature_RNA < 3500 & nCount_RNA < 15000 & percent.mt < 25)
trophall <- subset(trophall, subset = nFeature_RNA > 100 & nCount_RNA > 1000)
trophall #number of cells after filtering
plot <- VlnPlot(trophall, pt.size = 0.10, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), group.by = "sample", ncol = 3)

SaveObject(trophall, "trophall_afterfiltering")
trophall <- ReadObject("trophall_afterfiltering")

#_______________________________________________________________________________

#Normalizing the data
#then multiplies this by a scale factor (10,000 by default), and log-transforms the result. Normalized values are stored in pbmc[["RNA"]]@data
trophall <- NormalizeData(trophall, normalization.method = "LogNormalize", scale.factor = 10000)

#Identification of highly variable features
trophall <- FindVariableFeatures(trophall, selection.method = "vst", nfeatures = 2000)
trophall[["RNA"]]@meta.features %>% arrange(desc(vst.mean)) #list of transcripts ordered from highest to lowest total expression across all cells, can change sort

#top 10 most variable genes
top10 <- head(VariableFeatures(trophall), 10)
top20 <- head(VariableFeatures(trophall), 20)
top100 <- head(VariableFeatures(trophall), 50)

#plot variable features with and without labels
plot1 <- VariableFeaturePlot(trophall) #without
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE) #with_10
plot3 <- LabelPoints(plot = plot1, points = top20, repel = TRUE, max.overlaps = 100) #with_20
plot4 <- LabelPoints(plot = plot1, repel = TRUE, points = top100 , max.overlaps = 100)
plot1 + plot2 #1 x 2 layout
plot1 + plot3 #1 x 2 layout

#scale the data (transform) prior to dimensionality reduction (scales mean to 0 and variance across cells to 1). This is an essential step for PCA.
trophall <- ScaleData(trophall)

#perform linear dimensional reduction (PCA)
trophall <- RunPCA(trophall)

SaveObject(trophall, "seurat_obj_after_PCA")
trophall <- ReadObject("seurat_obj_after_PCA")

#_______________________________________________________________________________

#examine and visualize PCA results a few different ways
print(trophall[["pca"]], dims = 1:5, nfeatures = 5)

plot <- VizDimLoadings(trophall, dims = 1:5, reduction = "pca", ncol = 5)
SaveFigure(plot, "viz_PCA_loadings", width = 10, height = 8)

plot <- DimPlot(trophall, reduction = "pca", group.by = "sample")
SaveFigure(plot, "pc1_2_scatter", width = 8, height = 6)

plot <- DimHeatmap(trophall, dims = 1, cells = 500, balanced = TRUE, fast = FALSE) #this heatmap allows for easy exploration of the heterogeneity in a dataset. 
SaveFigure(plot, "dim_heatmap1", width = 8, height = 6)

plot <- DimHeatmap(trophall, dims = 1:15, cells = 500, balanced = TRUE) #this shows dimensions 1 through 15, each of these is a PC
SaveFigure(plot, "dim_heatmap1_15", width = 8, height = 6)

#Determine dimensionality of a dataset. Assigns a signficiance to each principle component. 
#"Significance" by this re-sampling randomly permutated test is defined as PCs that have strong enrichment of low-pvalue features (genes)
trophall <- JackStraw(trophall, num.replicate = 100)
trophall <- ScoreJackStraw(trophall, dims = 1:20)
plot <- JackStrawPlot(trophall, dims = 1:15) #displays pvalues of PCs
SaveFigure(plot, "PC_jackstraw_plot", width = 8, height = 10)

plot <- ElbowPlot(trophall) #alternatively, ElbowPlot makes a ranking of PCs based on percentage variance explained by each one. 
SaveFigure(plot, "PC_elbow_plot", width = 8, height = 10)

#Clustering
#See Seurat tutorial for more about clustering strategy 
trophall <- FindNeighbors(trophall, dims = 1:20)
trophall <- FindClusters(trophall, resolution = 0.25) 
head(Idents(trophall), 5) #Look at cluster IDs for first 5 cells
trophall <- BuildClusterTree(trophall, reorder = TRUE, reorder.numeric = TRUE)

#Run non-linear dimensional reduction (UMAP/tSNE)
trophall <- RunUMAP(trophall, dims = 1:20) #use the same dims here as you used in PCA
plot <- DimPlot(trophall, reduction = "umap", label = TRUE, label.size = 4, pt.size = 0.5)

SaveObject(trophall, "trophall_afterUMAP")
trophall <- ReadObject("trophall_afterUMAP")

#_______________________________________________________________________________
#finding differentially expressed features
trophall@active.ident #current identities

#various ways to switch active idents
trophall <- SetIdent(trophall, value = trophall$collapsed) #set to cluster identities
new_ids <- c("TS-STB", "proto TS-STB 1", "TS-EVT", "TS-EVT predisposed", "proto TS-STB 2", "TS-CTB")
trophall <- RenameIdents(trophall, new_ids) #proper ident syntax/grammar

trophall <- SetIdent(trophall, value = trophall$TgGRA1) #set to TgGRA1+/- identity

trophall <- SetIdent(trophall, value = trophall$sample) #set to sample identity
new_ids <- c("TS-STB", "TS-EVT", "TS-CTB", "infectedTS-CTB")
trophall <- RenameIdents(trophall, new_ids) #proper ident syntax/grammar

trophall@active.ident

#finding cluster markers
markers <- FindAllMarkers(trophall, min.pct = 0.25, logfc.threshold = 0.25)
markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
markers

#RidgePlot
plot <- RidgePlot(trophall, features = c("CGA", "EGFR", "TEAD1"), cols = c(STB, preSTB2, preSTB1, CTB, preEVT, EVT))
SaveFigure(plot, "RidgePlot_CGA_EGFR_TEAD1", height = 3, width = 11)

#dot plot of markers, currently in order of STB - CTB - EVT
plot <- DotPlot(trophall, features = c("CGA", "PAPPA", "CGB1", "CGB3", "CGB5", "CGB7", "KISS1", "SDC1", "CYP19A1",       
                                       "ERVW-1", "ERVFRD-1","PEG10","EGFR", "TFAP2C", "SP6", "OVOL1", "LYN", "TEAD4",
                                       "YAP1", "LRP2", "ITGA6", "TP63", "VGLL1", "CDH1", "CDH2", "NOTCH1",  "ITGB4",  
                                       "rna_TgGRA1", "ITGA2", "BACH2", "BCAT1", "LMCD1", "MKI67","TOP2A", "TEAD1", "TCF7L2", "GCM1", "FN1", "MMP2", 
                                       "HLA-G", "HTRA4", "DIO2", "ITGA5", "NOTUM")) + 
  #coord_flip() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position = "left")
plot
SaveFigure(plot, "dotplot_horizontal_6-24-24", width = 13, height = 3, res = 300)

#reorder clusters
trophall@active.ident <- factor(trophall@active.ident, levels = c("TS-STB", "proto TS-STB 2", "proto TS-STB 1", "TS-CTB", "TS-EVT predisposed", "TS-EVT"))

#subsetting by sample
infCTB_subset <- subset(trophall, subset = sample == "infCTB")
CTB_subset <- subset(trophall, subset = sample == "CTB")
STB_subset <- subset(trophall, subset = sample == "STB")
EVT_subset <- subset(trophall, subset = sample == "EVT")

markers <- FindMarkers(trophall, ident.1 = "TS-STB", ident.2 = "TS-EVT", min.pct = 0.25) #find markers of one sample versus another 
head(markers, n = 100) %>% arrange(avg_log2FC) #view in window
write.csv(markers, "filename.csv") #will save to working directory

#_______________________________________________________________________________

#filtering doublets and reassigning to trophall
trophall <- subset(trophall, TgGRA1 > 0.0001 & sample == "SYN", invert = TRUE) #removing cells that are "SYN" and TgGRA1(+)
trophall <- subset(trophall, TgGRA1 > 0.0001 & sample == "EVT", invert = TRUE) #above, but "EVT"
trophall <- subset(trophall, TgGRA1 > 0.0001 & sample == "CYT", invert = TRUE) #above, but "CYT"

#subset TgGRA1+ cells
gra1pos <- subset(trophall, TgGRA1 > 1.5) #filtering trophall then subsetting TgGRA1(+)
AverageExpression(gra1pos, features = c("rna_TgGRA1")) #calculate average expression of TgGRA1 across idents
plot <- FeaturePlot(gra1pos, features = c("rna_TgGRA1")) #check expression across clusters

plot <- VlnPlot(gra1pos, features = c("rna_TgGRA1"), cols = c(preSTB2, preSTB1, CTB, preEVT), pt.size = 1) + ggtitle("TgGRA1 expression among infected cells")
SaveFigure(plot, "UMAP_gra1pos_across_clusters", height = 4, width = 4, res = 300)

SaveObject(gra1pos, "gra1pos_object")
gra1pos <- ReadObject("gra1pos_object")

#designate new metadata classification (ident) based on TgGRA1 expression = positive or negative
#this allows us to compare all TgGRA1(+) to all others and find markers
trophall.gra1pos <- WhichCells(trophall, expression = rna_TgGRA1 > 1.5)
trophall$TgGRA1 <- ifelse(colnames(trophall) %in% trophall.gra1pos, "TgGRA1-pos", "TgGRA1-neg")
trophall@active.ident 

#show expression of TgGRA1(+) across samples --> proxy for rate of infection/number of parasites inside cells
plot <- VlnPlot(gra1pos, features = c("rna_TgGRA1"), cols = c(preSTB2, preSTB1, CTB, preEVT), pt.size = 1) + ggtitle("TgGRA1 expression among infected cells")
SaveFigure(plot, "VlnPlot_TgGRA1amonginfcells_3-21-24", height = 6, width = 7)

#find markers of TgGRA1(+) cells among all
markers <- FindMarkers(trophall, ident.1 = "TgGRA1-pos", min.pct = 0.25)
head(markers, n = 100) %>% arrange(desc(avg_log2FC)) #view in window
write.csv(markers, "filename.csv") 

#find markers of TgGRA1(+) cells among infected CTB cluster
markers <- FindMarkers(infCYT_subset, ident.1 = "TgGRA1-pos", min.pct = 0.10) #find all the markers of one cluster among all cells
head(markers, n = 100) %>% arrange(avg_log2FC) #view in window
write.csv(markers, "filename.csv") 

#Create downsampled dataset
infCYT_downsample <- subset(infCYT_subset, downsample = 273)
DimPlot(infCYT_downsample) #visualize size of TgGRA1(+) sample versus TgGRA1(-) sample
table(infCYT_subset$TgGRA1) #compare to below
table(infCYT_downsample$TgGRA1) #compare to above
Idents(infCYT_downsample) <- infCYT_downsample$TgGRA1 #change idents according to what metadata you want visualized
plot <- DimPlot(infCYT_downsample, reduction = 'umap', cols = c("blue", "darkgrey"))
SaveFigure(plot, "filename", height = 5, width = 5.5)

#_______________________________________________________________________________

#regression line of SYN markers in infCYT between TgGRA1(+) vs TgGRA1(-) cells
avgmarkers <- AverageExpression(infCYT_downsample, slot = "counts", features = c("ERVW-1", "ERVFRD-1", "GCM1", "MSX2", 
                                                                                 "SDC1", "CGB3", "CGB5", "CGB7", "CGB1", "GDF15", 
                                                                                 "LEP", "TFPI", "COBLL1", "TREM1", "PSG3", "PSG5",
                                                                                 "HOPX", "TFPI2", "CYP19A1", "CGA")) 
df <- data.frame(avgmarkers)
plot <- ggplot(df, aes(x = RNA.TgGRA1.neg, y = RNA.TgGRA1.pos)) + 
  geom_point(size = 3) + 
  xlim(0,10) + ylim(0,10) +
  stat_smooth(method = "lm",
              formula = y ~ x,
              geom = "smooth", col = "red") +
  geom_label_repel(aes(label = rownames(df)),
                   segment.color = "grey50", max.overlaps = 25, min.segment.length = 1, box.padding = 0.2,
                   nudge_x = 0, nudge_y = 0.25, direction = "both", size = 2) +
  theme_classic() +
  stat_regline_equation(label.x = 5, label.y = 2, col = "red", size = 5) +
  ylab("Average expression TgGRA1(+) cells (counts)") + xlab("Average expression TgGRA1(-) cells (counts)") +
  ggtitle("SYN markers among infected TS-CYT") + theme(plot.title = element_text(hjust = 0.5, face = "bold"))
p <- p + geom_abline(intercept = 0, linetype = "dashed", col = "grey", size = 0.75)
p
p2 <- p + scale_x_break(c(2.5,5), scales = "free") 

stats <- lm(formula = df$RNA.TgGRA1.pos ~ df$RNA.TgGRA1.neg, data = df)
summary(stats) #F-statistic, R2, and p-value

SaveFigure(plot, "filename", width = 4, height = 4)

#_______________________________________________________________________________

#visualizations
par(mfrow = c(1, 1))
plotA <- FeaturePlot(trophall, features = c("ITGA2", "rna_TgGRA1"), blend = TRUE, blend.threshold = .1, cols = c("red", "blue"))
plotB <- FeaturePlot(trophall, features = c("MKI67", "rna_TgGRA1"), blend = TRUE, blend.threshold = .1, cols = c("red", "blue"))
plotC <- FeaturePlot(trophall, features = c("TEAD1", "rna_TgGRA1"), blend = TRUE, blend.threshold = .1, cols = c("red", "blue"))
plot <- ggarrange(plotA, plotB, plotC, nrow = 3)
SaveFigure(plot, "FeaturePlot_BLEND_ITGA2_MKI67_TEAD1_TgGRA1", height = 15, width = 17)


#Name clusters and change idents
new_ids <- c("TS-STB", "proto TS-STB 1", "TS-EVT", "TS-EVT predisposed", "proto TS-STB 2", "TS-CTB")
new_id_list <- list(TS.SYN = 1, 
                    proto.TS.SYN.1 = 2, 
                    TS.EVT = 3, 
                    TS.EVT.predisposed = 4, 
                    proto.TS.SYN.2 = 5,
                    TS.CYT = 6)
for (i in 1:length(new_id_list)) {
  ind <- which(trophall@meta.data$tree.ident %in% new_id_list[[i]])
  trophall@meta.data$collapsed[ind] <- names(new_id_list)[i]
}
trophall@meta.data$collapsed <- factor(
  trophall@meta.data$collapsed, levels = names(new_id_list), ordered = TRUE)
Idents(trophall) <- trophall@meta.data$collapsed
names(new_ids) <- levels(trophall)
trophall <- RenameIdents(trophall, new_ids) #switch back to other idents
head(trophall)
trophall@active.ident


plot <- DimPlot(trophall, reduction = "umap", cols = c("grey", "blue"), pt.size = 1)
SaveFigure(plot, "DimPlot_allcells_gra1posneg_4-10-24", height = 8, width = 8.5)



#WORKING AREA
gene <- "MKI67"
AverageExpression(trophall, features = (gene))
b <- FeaturePlot(trophall, features = c(gene))
AverageExpression(infCYT_subset, features = c(gene))
FeaturePlot(infCYT_subset, features = (gene))
AverageExpression(gra1pos, features = c(gene))
FeaturePlot(gra1pos, features = c(gene)) 
a+b
Idents(trophall)


markers <- FindMarkers(trophall, ident.1 = "TS-SYN", min.pct = 0.25) 
head(markers, n = 100) %>% arrange(desc(avg_log2FC))
markers <- FindMarkers(trophall, ident.1 = "proto TS-SYN 2", min.pct = 0.25)
head(markers, n = 100) %>% arrange(desc(avg_log2FC))
markers <- FindMarkers(trophall, ident.1 = "proto TS-SYN 1", min.pct = 0.25)
head(markers, n = 100) %>% arrange(desc(avg_log2FC))
markers <- FindMarkers(trophall, ident.1 = "TS-CYT", min.pct = 0.25)
head(markers, n = 100) %>% arrange(desc(avg_log2FC))
markers <- FindMarkers(trophall, ident.1 = "TS-EVT predisposed", min.pct = 0.25)
head(markers, n = 100) %>% arrange(desc(avg_log2FC))
markers <- FindMarkers(trophall, ident.1 = "TS-EVT", min.pct = 0.25)
head(markers, n = 100) %>% arrange(desc(avg_log2FC))

SaveObject(trophall, "trophall_7-28-23")
trophall <- ReadObject("trophall_7-28-23")
SaveObject(trophall, "trophall_final")
trophall <- ReadObject("trophall_final")

#EVT progenitor markers (in smooth chorion)
FeaturePlot(trophall, features = c("SP6")) #SUPYN
plot <- VlnPlot(trophall, features = c("ERVW-1", "ERVFRD-1",  
                                       #"MFSD2A", "ERVV-1", "ERVV-2",
                                       #"ERVK13-1", "SLC1A5", 
                                       "ERVH48-1"), ncol = 3) 
SaveFigure(plot, "Feature_ERVs_forqPCR", height = 6, width = 16)
FeaturePlot(trophall, features = c("ADAM12"))
FeaturePlot(trophall, features = c("TOP2A", "MKI67", "CTNNB1"))
genes <- "SLC2A3"
VlnPlot(trophall, features = c(genes)) | FeaturePlot(trophall, features = c(genes))

SaveObject(gra1pos, "gra1pos_7-28-23")
gra1pos <- ReadObject("gra1pos_7-28-23")


my_levels <- c("TS-STB", "proto TS-STB 2", "proto TS-STB 1", "TS-CTB", "TS-EVT predisposed", "TS-EVT") #reorder clusters
trophall@active.ident <- factor(x = trophall@active.ident, levels = my_levels)
gra1pos@active.ident <- factor(x = gra1pos@active.ident, levels = my_levels)
plot <- RidgePlot(trophall, features = c("CGA", "EGFR", "TEAD1"), cols = c(STB, preSTB2, preSTB1, CTB, preEVT, EVT))
SaveFigure(plot, "RidgePlot_CGA_EGFR_TEAD1_3-20-24", height = 3, width = 11)





#_______________________________________________________________________________

#starting Monocle3 integration
cds <- as.cell_data_set(trophall)
colData(cds)
fData(cds)
rownames(fData(cds))[1:10]
fData(cds)$gene_short_name <- rownames(fData(cds))
counts(cds)
recreate.partition <- rep(1, length(cds@colData@rownames))

#assign partitions
names(recreate.partition) <- cds@colData@rownames
recreate.partition <- as.factor(recreate.partition)
recreate.partition
cds@clusters$UMAP$partitions <- recreate.partition

#assign the cluster info
list_cluster <- trophall@active.ident
cds@clusters$UMAP$clusters <- list_cluster

#assign UMAP coordinate - cell embeddings
cds@int_colData@listData$reducedDims$UMAP <- trophall@reductions$umap@cell.embeddings

#plot
cluster.before.trajectory <- plot_cells(cds,
                                        color_cells_by = "cluster",
                                        label_cell_groups = FALSE,
                                        group_label_size = 5) +
  theme(legend.position = "right")

cluster.names <- plot_cells(cds,
                            color_cells_by = "seurat_clusters",
                            label_cell_groups = FALSE,
                            group_label_size = 5) +
  theme(legend.position = "right")

#this will compare the seurat clustering to the clustering that you forced monocle to recreate (from Seurat).
#No need to use second string if you haven't consolidated clusters and renamed them from Seurat

#learn trajectory graph
cds <- learn_graph(cds, use_partition = FALSE)

umap <-plot_cells(cds,
                  color_cells_by = "cluster",
                  label_branch_points = FALSE,
                  label_roots = TRUE,
                  label_leaves = TRUE,
                  group_label_size = 5,
                  label_cell_groups = TRUE,
                  trajectory_graph_color = "black",
                  trajectory_graph_segment_size = 2,
                  cell_size = 0.5)
umap <- umap + theme(legend.position = "right")
umap
SaveFigure(umap, "umap_pseudotime_trophall_test", height = 8, width = 8, res = 300)

#order the cells in pseudotime
cds <- order_cells(cds, reduction_method = "UMAP", root_cells = colnames(cds[,clusters(cds) == "TS-CYT"])) # == "" needs to be set to your root
SaveObject(cds, "cds_monocle3_trophall_test")

plot <- plot_cells(cds,
                   color_cells_by = "pseudotime", #change to "sample" or other col in colData(cds)
                   label_branch_points = FALSE,
                   label_roots = FALSE,
                   label_leaves = FALSE,
                   group_label_size = 10,
                   label_cell_groups = TRUE,
                   cell_size = 0.5,
                   show_trajectory_graph = FALSE,
                   trajectory_graph_segment_size = 2,
                   trajectory_graph_color = "black")
plot <- plot + theme(legend.position = "right")
plot
SaveFigure(plot, "umap_pseudotime_trophall_test_6-13", height = 5, width = 6, res = 300)

#cells ordered by monocle3 pseudotime
pseudotime(cds)
cds$monocle3_pseudotime <- pseudotime(cds)
data.pseudo <- as.data.frame(colData(cds))
data.pseudo
colData(cds)

ordered.pseudo <- ggplot(data.pseudo, aes(monocle3_pseudotime, reorder(ident, monocle3_pseudotime, median), fill = ident)) + geom_boxplot()
ordered.pseudo <- ordered.pseudo + xlab("pseudotime") + ylab("ident")
ordered.pseudo

plot
SaveFigure(plot, "pseudotime_bypseudotime_umap_root3", width = 8, height = 8, res = 300)

ordered.pseudo 
SaveFigure(ordered.pseudo, "pseudotime_bycluster_boxplot_root3", width = 5, height = 8, res = 300)

#finding genes that change as a function of pseudotime
deg_troph <- graph_test(cds, neighbor_graph = "principal_graph", cores = 6)
deg_troph

deg_troph_ordered <- deg_troph %>%
  arrange(q_value) %>%
  filter(status == "OK")

deg_troph_ordered
write.csv(deg_troph_ordered, file = "pseudotime_deg_troph_test.csv")

FeaturePlot(trophall, features = c("DDAH1", "JPH2", "LMCD1", "TgGRA1")) #features from deg_troph_ordered

trophall$pseudotime <- pseudotime(cds)

pseudotime_seurat_plot <- FeaturePlot(trophall, features = "pseudotime", pt.size = 1) + ggtitle("pseudotime_seurat")
SaveFigure(pseudotime_seurat_plot, "pseudotime_seurat_umap_trophall_test", width = 8, height = 8, res = 300)

#3d trajectory
cds_3d <- reduce_dimension(cds, preprocess_method = c("PCA"), reduction_method = c("UMAP"), max_components = 3)
cds_3d <- cluster_cells(cds_3d)
cds_3d <- learn_graph(cds_3d)
cds_3d <- order_cells(cds_3d, root_pr_nodes=get_earliest_principal_node(cds))
cds_3d_plot_obj <- plot_cells_3d(cds_3d, color_cells_by="monocle3_pseudotime", show_trajectory_graph = FALSE)
plot <- plot_cells_3d(cds_3d, genes = c("TFPI"), min_expr = 0.0001) #equivalent of feature plot
SaveObject(cds, "cds_trophall_test")
SaveObject(cds_3d, "cds_3d_trophall_test")

htmlwidgets::saveWidget(plot, "plotly_TgGRA1.html") #save as interactive html file. file will only be usable from its origin save location
htmlwidgets::saveWidget(cds_3d_plot_obj, "plotly_pseudotime_trophall_test.html")

#_______________________________________________________________________________


#I've avoided naming the clusters up until now for the sake of not obscuring nuance. 
#Here, I name them to use in Scanpy for dendrogram analysis

new_ids <- c("TS-STB", "proto TS-STB 1", "TS-EVT", "TS-EVT predisposed", "proto TS-STB 2", "TS-CTB")

new_id_list <- list(TS.STB = 1, 
                    proto.TS.STB.1 = 2, 
                    TS.EVT = 3, 
                    TS.EVT.predisposed = 4, 
                    proto.TS.STB.2 = 5,
                    TS.CTB = 6)

for (i in 1:length(new_id_list)) {
  ind <- which(trophall@meta.data$tree.ident %in% new_id_list[[i]])
  trophall@meta.data$collapsed[ind] <- names(new_id_list)[i]
}

trophall@meta.data$collapsed <- factor(
  trophall@meta.data$collapsed, levels = names(new_id_list), ordered = TRUE)
Idents(trophall) <- trophall@meta.data$collapsed
names(new_ids) <- levels(trophall)
trophall <- RenameIdents(trophall, new_ids)
trophall@active.ident

plot <- DimPlot(trophall, reduction = "umap", label = TRUE, cols = c(STB, preSTB2, preSTB1, CTB, preEVT, EVT))
SaveFigure(plot, "DimPlot_correctlabeling_3-20-24", height = 6, width = 7)
SaveObject(trophall, "trophall_final")

#converting seurat object to anndata format via tutorial at: https://mojaveazure.github.io/seurat-disk/articles/convert-anndata.html
trophall <- ReadObject("trophall_final")
SaveH5Seurat(trophall, filename = "trophall2.h5Seurat")
seurat_ad <- Convert("trophall2.h5Seurat", dest = "h5ad")


#slingshot pseudotime analysis
if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("kstreet13/slingshot")
library(slingshot)
library(RColorBrewer)
BiocManager::install("zellkonverter")
library(zellkonverter)

sce <- readH5AD("trophall2.h5ad")
sds <- slingshot(reducedDim(sce, "umap"), clusterLabels = trophall$seurat_clusters, start.clus = 3, spar = 1.1)

a <- FeaturePlot(trophall, features = c("ITGA2", "SDC1"), ncol = 2, pt.size = 0.75)
SaveFigure(a, "Feature_ITGA2_SDC1", height = 8, width = 15, res = 300)

library(Seurat)
library(SeuratData)
library(SeuratDisk)
library(reticulate)

SaveH5Seurat(trophall, filename = "trophall_test_h5Seurat")
trophall_h5ad <- Convert("trophall_test_h5Seurat.h5seurat", dest = "h5ad")
ad <- import("anndata", convert = FALSE)

head(trophall)



