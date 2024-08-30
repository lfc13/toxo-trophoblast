#10X genomics multiome simultaneous scATAC and scRNA data analysis using Signac and Seurat

#convenience functions for saving objects and plots
SaveFigure <- function(plots, name, type = "png", width, height, res){
  if(type == "png") {png(paste0(fig_path, name, ".", type),
                         width = width, height = height, units = "in", bg = "transparent", res = 300) # add bg = "transparent", for transparent background
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
data_path <- "OneDrive - University of Pittsburgh/Leah Cabo/Leah PhD/Data/single cell/10x multiome/Robjects/"
fig_path <- "OneDrive - University of Pittsburgh/Leah Cabo/Leah PhD/Data/single cell/10x multiome/figures/"
macs2.path <- "miniconda3/envs/PeakCalling_analysis/bin/macs2"

#installation of Signac according to stuartlab.org
setRepositories(ind=1:3)
install.packages("Signac")
install.packages("SeuratObject") #had to update this package to load Signac
#successful

#installing other packages
install.packages("biovizBase")
BiocManager::install("Herper")
install_CondaTools(tools="macs2", env="PeakCalling_analysis", pathToMiniConda="miniconda3/") #miniconda is already installed so use this to install tools into miniconda through R without re-installing miniconda or routing through python
install.packages("Matrix", type = "source")
install.packages("irlba", type = "source")
install.packages("RSpectra")
install.packages("JASPAR2020")
BiocManager::install("JASPAR2020")
BiocManager::install("TFBSTools")
BiocManager::install("motifmatchr")
BiocManager::install("ggplot2")
BiocManager::install("ggseqlogo", force = TRUE)
install.packages("scCustomize")
BiocManager::install("GreenleafLab/chromVAR")
devtools::install_github("rpolicastro/scProportionTest")
remotes::install_github("bcgov/elucidate")
install.packages("scCustomize")

#update Seurat
remotes::install_github("satijalab/seurat", "seurat5", quiet = TRUE)
#igraph had non-zero exit status during dependency updates but Seurat loaded fine

#install genome assembly and annotation packages
BiocManager::install(c('BSgenome.Hsapiens.UCSC.hg19', 'EnsDb.Hsapiens.v75')) #human hg19
BiocManager::install(c('BSgenome.Hsapiens.UCSC.hg38', 'EnsDb.Hsapiens.v86')) #human hg38
remotes::install_github("stuart-lab/signac", ref = "develop")

#download filtered feature bc matrix, atac fragments, and atac fragments index from 10x cloud
#done via CyberDuck into OneDrive

#load packages
library(Signac)
library(Seurat)
library(SeuratObject)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(biovizBase)
library(Herper)
library(Matrix)
library(irlba)
library(dplyr)
library(ggplot2)
library(JASPAR2020)
library(TFBSTools)
library(patchwork)
library(monocle3)
library(SeuratWrappers)
library(ggseqlogo)
library(chromVAR)
library("scProportionTest")
library(elucidate)
library(scCustomize)

#load RNA and ATAC data
counts <- Read10X_h5("OneDrive - University of Pittsburgh/Leah Cabo/Leah PhD/Data/single cell/10x multiome/files_for_Signac/filtered_feature_bc_matrix.h5")
fragpath <- "OneDrive - University of Pittsburgh/Leah Cabo/Leah PhD/Data/single cell/10x multiome/files_for_Signac/new_atac_fragments.tsv.gz"

head(counts)

#get annotations for hg38
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
#seqlevels(annotation) <- paste0("chr", seqlevels(annotation)) #this step is present in the vignette and adds the prefix "chr" to the chromosome category but puts the annotations in the chromosome categories out of agreement

seqlevelsStyle(multiome)
seqlevelsStyle(annotation) <- "UCSC"
#seqlevelsStyle(annotation) <- "NCBI"
seqlevels(annotation) # <- check chromosome names

#create a Seurat object containing the RNA data
multiome <- CreateSeuratObject(counts = counts$`Gene Expression`, assay = "RNA", project = "multiome", names.delim = "-", names.field = 2) #names.delim and names.field allow you to split into sample types. aggr assigns each barcode a "-#" for which sample it came from
Idents(multiome) <- 'orig.ident'
multiome@active.ident
multiome <- RenameIdents(multiome, 
                         `1` = "Tg-infected TS-CTB", 
                         `2` = "mock-infected TS-CTB", 
                         `3` = "Day4 TS-EVT", 
                         `4` = "Day8 TS-EVT") #rename idents based on sample name
#this only temporarily changes the IDs

#mismatch in dataframe size
new_ids <- c("Day8 TS-EVT", "Day4 TS-EVT", "mock-infected TS-CTB", "Tg-infected TS-CTB")
new_id_list <- list(D8.TS.EVT = 1, 
                    D4.TS.EVT = 2, 
                    MI.TS.CTB = 3, 
                    TGI.TS.CTB = 4)
for (i in 1:length(new_id_list)) {
  ind <- which(multiome@meta.data$orig.ident %in% new_id_list[[i]])
  multiome@meta.data$sample[ind] <- names(new_id_list)[i]
}
multiomel@meta.data$sample <- factor(
  multiome@meta.data$sample, levels = names(new_id_list), ordered = TRUE)
Idents(multiome) <- multiome@meta.data$sample
names(new_ids) <- levels(multiome)


#create ATAC assay and add it to the object
multiome[["ATAC"]] <- CreateChromatinAssay(counts = counts$Peaks, sep = c(":", "-"), 
                                           fragments = fragpath, 
                                           annotation = annotation)
                                           #tolerance = 0.95) <- this would redudce the number of cells that have to match but was not used during this analysis

#ran into errors at this point: not all cells request found in fragment file is the main error
#attempt 1: changed "fragments" to exact file path instead of fragpath <- did not fix
#attempt 2: changed "fragments" to .tbi folder <- error = fragment file is not indexed <- did not fix
#attempt 3: changed fragpath to be the outer folder with both .tbi and .tsv files <- did not fix
#attempt 4: changed fragpath back to be .tsv file and reduced tolerance slot to 0.95 <- did not fix
#attempt 5: unzip (gzip), rezip (bgzip), re-index (tabix) fragment file in case of damage between server transfers. Did this on toxohunter and retransferred files using Fetch
  #see suggestion source at https://github.com/stuart-lab/signac/issues/242#issuecomment-693666737
  #note that I cleared all previous file versions out of the /files_for_Signac/ directory and had both the fragment file and index in the same directory (suggestion on github)
  #THIS WORKED HALLELUJAH

DefaultAssay(multiome) <- "ATAC"

multiome <- NucleosomeSignal(multiome)
multiome <- TSSEnrichment(multiome) #error: ecdf(x = object$TSS.enrichment) 'x' must have 1 or more non-missing values, error avoided after not adding prefix to seqlevels of annotation (line 56)

head(Fragments(multiome)[[1]])
head(Annotation(multiome))
#chromosome annotations need to match in each of these (NCBI style versus UCSC style)

SaveObject(multiome, "multiome_afterATAC")
multiome <- ReadObject("multiome_afterATAC")

#checking if chromosome annotation is the same to troubleshoot TSSEnrichment error
head(Annotation(multiome))
chroms <- unique(Annotation(multiome))
chroms@seqnames #shows chromosome names

str(multiome)
head(Fragments(multiome)[[1]], n = 10000)

VlnPlot(object = multiome,
        features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal"),
        ncol = 4,
        pt.size = 0)

# filter out low quality cells, lost 11.5K cells here, consider adjusting
multiome <- subset(
  x = multiome, subset = 
    nCount_ATAC < 100000 &
    nCount_RNA < 50000 &
    nCount_ATAC > 1000 &
    nCount_RNA > 1000 &
    nucleosome_signal < 2 &
    TSS.enrichment > 1)

VlnPlot(object = multiome,
        features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal"),
        ncol = 4,
        pt.size = 0)

SaveObject(multiome, "multiome_afterfilter")
multiome <- ReadObject("multiome_afterfilter")

#peak calling

#call using MACS2, this requires installation of macs2 with python using the herper package in R (kind of like reticulate)
peaks <- CallPeaks(multiome, macs2.path = macs2.path) #macs2 function path called by going within /bin/ inside new peak calling env see line 26
#fyi this step takes >20 minutes

#remove peaks on nonstandard chromosomes and in genomic blacklist regions
peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_hg38_unified, invert = TRUE)

#quantify coutns in each peak
macs2_counts <- FeatureMatrix(fragments = Fragments(multiome),
                              features = peaks,
                              cells = colnames(multiome))

#create new assay using the MACS2 peak set and add it to the Seurat obejct
multiome[["peaks"]] <- CreateChromatinAssay(counts = macs2_counts,
                                            fragments = fragpath,
                                            annotation = annotation)

SaveObject(multiome, "multiome_afterMACS")
multiome <- ReadObject("multiome_afterMACS")

#gene expression and data processing
DefaultAssay(multiome) <- "RNA"
multiome <- SCTransform(multiome) #vector memory exhausted. Followed suggestion at: https://stackoverflow.com/questions/51295402/r-on-macos-error-vector-memory-exhausted-limit-reached
multiome <- RunPCA(multiome)

#dna accessibility data processing (performing latent semantic indexing)
DefaultAssay(multiome) <- "peaks"
multiome <- FindTopFeatures(multiome, min.cutoff = 5)
multiome <- RunTFIDF(multiome)
multiome <- RunSVD(multiome)
#Error in irlba(A = t(x = object), nv = n, work = irlba.work, tol = tol) :function 'as_cholmod_sparse' not provided by package 'Matrix'
#attempt 1: manually install Matrix and irlba from source and restart R <- This worked! did not need to do below
#attempt 2: install and run Rspectra in preference to irlba for spectral initialization: https://github.com/jlmelville/uwot/issues/115 

#build a joint neighbor graph using both assays
multiome <- FindMultiModalNeighbors(object = multiome,
                                reduction.list = list("pca", "lsi"), 
                                dims.list = list(1:50, 2:40),
                                modality.weight.name = "RNA.weight",
                                verbose = TRUE)

SaveObject(multiome, "multiome_WNN_5-7-24")
multiome <- ReadObject("multiome_WNN_5-7-24")

multiome <- NormalizeData(multiome, normalization.method = "LogNormalize", scale.factor = 10000) #have default assay on "RNA" to generate data layer
cluster1 <- NormalizeData(cluster1)

gene.activities <- GeneActivity(multiome)
multiome[['RNA']] <- CreateAssayObject(counts = gene.activities)

multiome@assays
DefaultAssay(multiome) <- "peaks"
DefaultAssay(multiome) <-"SCT"
DefaultAssay(multiome) <- "RNA"
DefaultAssay(multiome) <- "ATAC"

#build a joint UMAP visualization <- at this point idents(SeuratProject) so will be all cells
multiome <- RunUMAP(object =  multiome,
                nn.name = "weighted.nn",
                assay = "RNA",
                verbose = TRUE)

multiome <- FindClusters(multiome, graph.name = "wsnn", algorithm = 2, verbose = FALSE, resolution = 0.06)
#algorithm options: 1 = original Louvain, 2 = Louvain with multilevel refinement, 3 = SLM, 4 = Leiden (requires leidenalg python)
plot <- DimPlot(multiome, label = TRUE, repel = TRUE)
SaveFigure(plot, "dimplot_bycluster_5-16-24", height = 8, width = 9)

#remove clustering at other resolutions
multiome@meta.data
multiome$wsnn_res.0.08 <- NULL
multiome$wsnn_res.0.06 <- NULL
multiome$wsnn_res.0.1 <- NULL
multiome$wsnn_res.0.05 <- NULL

plot <- FeaturePlot(multiome, features = c("DKK1"))
SaveFigure(plot, "featureplot_GRA1_2-29-24", height = 8, width = 8)

gene <- "CDH2"
a <- FeaturePlot(multiome, features = c(gene))
b <- FeaturePlot(multiome, features = c("rna_TgGRA1"))
b <- VlnPlot(multiome, features = c(gene))
c <- RidgePlot(multiome, features = c(gene))
plot <- a + b + c
plot
SaveFigure(plot, "featureplot_ridgeplot_ITGA2_2-29-24", height = 15, width = 8)

SaveObject(multiome, "multiome_runUMAP_findClusters")
multiome <- ReadObject("multiome_runUMAP_findClusters")

#initial plotting
DefaultAssay(multiome) <-"RNA"
DimPlot(multiome, label = TRUE, repel = TRUE)
FeaturePlot(multiome, features = c("TgGRA1"), pt.size = 1)

#optional Seurat specific processing
multiome <- FindVariableFeatures(multiome, selection.method = "vst", nfeatures = 2000)
top20 <- head(VariableFeatures(multiome), 20)
plot <- VariableFeaturePlot(multiome)
plot1 <- LabelPoints(plot = plot, points = top20, repel = TRUE)
all.genes <- rownames(multiome)
multiome <- ScaleData(multiome, features = all.genes)

#find cluster markers
head(Idents(multiome))
cluster0markers <- FindMarkers(multiome, ident.1 = 0)
head(cluster0markers, n = 100) %>% arrange(avg_log2FC)
cluster1markers <- FindMarkers(multiome, ident.1 = 1)
head(cluster1markers, n = 100) %>% arrange(avg_log2FC)
cluster2markers <- FindMarkers(multiome, ident.1 = 2)
head(cluster2markers, n = 100) %>% arrange(avg_log2FC)
cluster3markers <- FindMarkers(multiome, ident.1 = 3)
head(cluster3markers, n = 100) %>% arrange(avg_log2FC)
cluster4markers <- FindMarkers(multiome, ident.1 = 4)
head(cluster4markers, n =100) %>% arrange(avg_log2FC)
cluster5markers <- FindMarkers(multiome, ident.1 = 5)
head(cluster5markers, n =100) %>% arrange(avg_log2FC)
cluster7markers <- FindMarkers(multiome, ident.1 = 7)
head(cluster7markers, n =100) %>% arrange(avg_log2FC)
cluster8markers <- FindMarkers(multiome, ident.1 = 8)
head(cluster8markers, n =100) %>% arrange(avg_log2FC)
cluster9markers <- FindMarkers(multiome, ident.1 = 9)
head(cluster9markers, n =100) %>% arrange(avg_log2FC)

cluster1_tgGRA1posvneg <- FindMarkers(cluster1, ident.1 = "TgGRA1-pos", ident.2 = "TgGRA1-neg")
head(cluster1_tgGRA1posvneg, n = 100) %>% arrange(avg_log2FC)
write.csv(cluster1_tgGRA1posvneg, "OneDrive - University of Pittsburgh/Leah Cabo/Leah PhD/Data/single cell/10x multiome/markerscluster1_TgGRA1posvneg_5-8-24.csv")

allcells_tgGRA1posvneg <- FindMarkers(multiome, ident.1 = "TgGRA1-pos", ident.2 = "TgGRA1-neg")
head(allcells_tgGRA1posvneg, n = 100) %>% arrange(avg_log2FC)
write.csv(allcells_tgGRA1posvneg, "OneDrive - University of Pittsburgh/Leah Cabo/Leah PhD/Data/single cell/10x multiome/markersallcells_TgGRA1posvneg_5-8-24.csv")

plot <- VlnPlot(gra1pos, features = c("rna_TgGRA1"), cols = c("#F6766D", "#A4A601", "#3CC07C", "#E76BF3"))
SaveFigure(plot, "VlnPlot_TgGRA1acrosscluster_5-16-24", height = 7, width = 8)

#linking peaks to genes
DefaultAssay(multiome) <- "peaks"
DefaultAssay(multiome_c1) <- "RNA"

#first compute the GC content for each peak
seqlevelsStyle(BSgenome.Hsapiens.UCSC.hg38) <- "NCBI" #need to run this for motif analysis
seqnames(BSgenome.Hsapiens.UCSC.hg38)

#Identity setting

#make a new metadata classification of TgGRA1(+) versus TgGRA1(-)
multiome.gra1pos <- WhichCells(multiome, expression = rna_TgGRA1 > 2)
multiome$TgGRA1 <- ifelse(colnames(multiome) %in% multiome.gra1pos, "TgGRA1-pos", "TgGRA1-neg")
multiome <- SetIdent(multiome, value = multiome$TgGRA1) #set to gra1+/- identity
Idents(multiome)

SaveObject(multiome, "multiome_afterGRA1")
multiome <- ReadObject("multiome_afterGRA1")

#subset cluster 1 and set idents to TgGRA1+/-
plot <- DimPlot(multiome, reduction = "umap", label = F, cols = c("lightgrey", "blue")) #check cluster numbers
multiome_c1 <- subset(multiome, idents = 1) #subset cluster 1
multiome_c1 <- SetIdent(multiome_c1, value = multiome_c1$TgGRA1) #set to gra1+/- identity
SaveFigure(plot, "DimPlot_TgGRA1pos_5-15-24", width = 8.5, height = 8)

#subset infected CTB sample
TGI_CTB <- subset(multiome, orig.ident == "1")
DimPlot(TGI_CTB, reduction = "umap")
TGI_CTB <- SetIdent(TGI_CTB, value = TGI_CTB$TgGRA1)
TGI_CTB <- SetIdent(TGI_CTB, value = TGI_CTB$wsnn_res.0.06)

MI_CTB <- subset(multiome, orig.ident == "2")
D4_EVT <- subset(multiome, orig.ident == "3")
D8_EVT <- subset(multiome, orig.ident == "4")

current <- D8_EVT
new_ids <- c("spontaneous TS-STB", "TS-CTB", "TS-EVT","abberant TS-EVT", "proto TS-STB")
names(new_ids) <- levels(current)
current <- RenameIdents(object = current, new_ids)
current@active.ident
plot <- DimPlot(current, reduction = "umap")
SaveFigure(plot, "DimPlot_D8_EVT_5-21-24", height = 5, width = 5.5)

table(multiome@meta.data$wsnn_res.0.06)
table(multiome@meta.data$wsnn_res.0.06, multiome@meta.data$TgGRA1)

#remove cluster 6 because it's literally two cells and is misleading
multiome <- subset(multiome, wsnn_res.0.07 == "6", invert = TRUE)
DimPlot(multiome, reduction = 'umap')

#switching between identities
multiome <- SetIdent(multiome, value = multiome$wsnn_res.0.06) #set to cluster identity
multiome <- SetIdent(multiome, value = multiome$TgGRA1) #set to gra1+/- identity
multiome <- SetIdent(multiome, value = multiome$orig.ident) #set to sample identity
multiome <- RenameIdents(multiome, 
                         `1` = "Tg-infected TS-CTB", 
                         `2` = "mock-infected TS-CTB", 
                         `3` = "Day4 TS-EVT", 
                         `4` = "Day8 TS-EVT") #rename orig.idents based on sample name, only temporary
multiome@active.ident

#PEAK LINKING
DefaultAssay(multiome) <- "peaks"
DefaultAssay(multiome) <- "RNA"
DefaultAssay(cluster1) <- "peaks"
DefaultAssay(cluster1) <- "RNA"

VlnPlot(cluster1, features = c("CHD1", "CHD2", "SMARCA1", "SMARCA2", "SMARCA4", "SMARCA5", "EP400", "INO80"), pt.size = 0, sort = "increasing")

multiome <- RegionStats(multiome, genome = BSgenome.Hsapiens.UCSC.hg38)

#link peaks to genes
cluster1 <- LinkPeaks(object = cluster1, peak.assay = "peaks", expression.assay = "SCT", genes.use = c("DKK1"))

SaveObject(multiome, "multiome_afterLinkPeaks")
multiome <- ReadObject("multiome_afterLinkPeaks")

gra1pos <- subset(multiome, rna_TgGRA1 > 2)
gra1pos <- SetIdent(gra1pos, value = gra1pos$wsnn_res.0.06)
plot <- DimPlot(gra1pos, cols = c("#F6766D", "#A4A601", "#3CC07C", "#E76BF3"))
SaveFigure(plot, "DimPlot_subsetgra1pos_5-21-24", height = 8, width = 8.5)
VlnPlot(gra1pos, features = c("rna_TgGRA1"), sort = "increasing")
SaveFigure

new_ids <- c("spontaneous TS-STB", "TS-CTB", "TS-EVT", "proto TS-STB")
names(new_ids) <- levels(gra1pos)
gra1pos <- RenameIdents(object = gra1pos, new_ids)
gra1pos@active.ident

idents.plot <- c("TgGRA1-neg", "TgGRA1-pos")
idents.plot <- c("0", "1", "2", "3", "4")
multiome@active.ident

VlnPlot(multiome, features = c("DKK1"))
FeaturePlot(multiome, features = c("DKK1"))
FeaturePlot(multiome, features = c("DKK1"))
plot <- CoveragePlot(object = multiome,
                     region = "DKK1",
                     features = "DKK1",
                     expression.assay = "SCT",
                     idents = idents.plot,
                     links = TRUE,
                     peaks = TRUE)
                     #extend.upstream = 0,
                    #extend.downstream = 0)


plot <- CoverageBrowser(object = multiome,
                     region = "FOS",
                     features = "FOS",
                     expression.assay = "SCT",
                     idents = idents.plot,
                     extend.upstream = 10000,
                     extend.downstream = 10000)

SaveFigure(plot, "CoveragePlot_NRG1_cluster1subset_5-8-24", height = 4, width = 10)

#TF MOTIF ENRICHMENT ANALYSIS
#adding motif analysis information to Seurat object
#get a list of motif position frequecy matrices from the JASPAR database
pfm <- getMatrixSet(x = JASPAR2020, opts = list(collection = "CORE", tax_group = "vertebrates", all_versions = FALSE))

#add motif information
seqlevelsStyle(BSgenome.Hsapiens.UCSC.hg38) <- "NCBI" #this will match the reference genome seqlevels to the multiome object seqlevels
multiome <- AddMotifs(object = multiome, genome = BSgenome.Hsapiens.UCSC.hg38, pfm = pfm) #needed to install motifmatchr

#finding overrepresented motifs between identities
Idents(multiome)
Idents(multiome_c1)
DimPlot(multiome, reduction = "umap")
da_peaks <- FindMarkers(object = multiome,
                        ident.1 = 3,
                        #ident.2 = 1,
                        only.pos = "TRUE",
                        test.use = "LR",
                        min.pct = 0.05,
                        latent.vars = "nCount_peaks")

SaveObject(da_peaks, "da_peaks_multiomecluster3")
da_peaks <- ReadObject("da_peaks_TgGRA1-pos_cluster1subset")

da_peaks %>% arrange(desc(avg_log2FC))

#get top differentially accessible peaks
top.da.peak <- rownames(da_peaks[da_peaks$p_val<0.0005,])

FeaturePlot(multiome, features = c("ERVV-1", "ERVV-2"))

#test enrichment
enriched.motifs <- FindMotifs(object = multiome, features = top.da.peak)

#plot the position weight matrices for motifs

plot <- MotifPlot(object = multiome, motifs = head(rownames(enriched.motifs), n = 25))
SaveFigure(plot, "motifplot_n=25_multiomecluster3_EVTs_3-19-24", height = 7, width = 12)

SaveObject(multiome, "multiome_afterAddMotif_5-8-24")
SaveObject(cluster1, "cluster1_beforeMotifFootprinting_5-8-24")
cluster1 <- ReadObject("cluster1_beforeMotifFootprinting_5-8-24")
multiome <- ReadObject("multiome_afterAddMotif")

#Motif footprinting
multiome@active.ident
cluster1 <- Footprint(object = cluster1, motif.name = c("BACH2"), genome = BSgenome.Hsapiens.UCSC.hg38, in.peaks = TRUE)
plot <- PlotFootprint(cluster1, features = c("BACH2"))
SaveFigure(plot, "Footprint_ZNF384+BACH2_cluster1_gra1posneg", height = 5, width = 11, res = 300)
data <- GetFootprintData(cluster1, features = c("ZNF384"))
write.csv(data, "OneDrive - University of Pittsburgh/Leah Cabo/Leah PhD/Data/single cell/10x multiome/Footprintdata_cluster1_ZNF384.csv")

#motif activity

multiome <- RunChromVAR(object = multiome, genome = BSgenome.Hsapiens.UCSC.hg38)
DefaultAssay(multiome) <- "chromvar"
plot <- FeaturePlot(multiome, features = c("MA0497.1"), min.cutoff = "q10", max.cutoff = "q90")


#Working room
#switching between identities
multiome <- SetIdent(multiome, value = multiome$wsnn_res.0.07) #set to cluster identity
multiome <- SetIdent(multiome, value = multiome$TgGRA1) #set to gra1+/- identity
multiome@active.ident

plot <- FeaturePlot(multiome, features = c("VGLL3", "rna_TgGRA1"), blend = TRUE, blend.threshold = .1, cols = c("red", "blue"))
plot
SaveFigure(plot, "FeaturePlot_BLEND_BACH2_TgGRA1", height = 5, width = 17)

#PEAK LINKING
DefaultAssay(multiome) <- "peaks"
DefaultAssay(multiome) <- "RNA"
DefaultAssay(multiome) <- "ATAC"
DefaultAssay(multiome_c1) <- "peaks"
DefaultAssay(multiome_c1) <- "RNA"
DefaultAssay(multiome_c1) <- "ATAC"
DefaultAssay(cluster1) <- "RNA"
DefaultAssay(cluster1) <- "peaks"
DefaultAssay(cluster1) <- "ATAC"

FeaturePlot(multiome, feature = c("NPPB"))
VlnPlot(multiome, features = c("BACH2","DKK1","JDP2", "TGFB2", "NRG1", "IL6", "CXCL8"), pt.size = 0)
  #coord_flip() + 
  #theme(axis.text.x = element_text(angle = 45, hjust = 1))

markers <- FindMarkers(multiome, ident.1 = "3", min.pct = 0.25) 
write.csv(markers, "markers_TgGRA1+_multiome.csv")

head(markers, n = 100) %>% arrange(avg_log2FC)

cluster1$
cluster1 <- subset(multiome, ident = "1")
cluster1 <- SetIdent(cluster1, value = cluster1$TgGRA1)


plot <- VlnPlot(cluster1, features = c("CGA", "PAPPA", "CGB1", "CGB3", "CGB5", "CGB7", "KISS1",  "OVOL1", "SDC1", "CYP19A1",       #currently in order of SYN - CYT - EVT
                                       "ERVW-1", "ERVFRD-1"), pt.size = 0)

plot <- VlnPlot(cluster1, features = c("PEG10","EGFR", "TFAP2C", "LYN", "TEAD4", "YAP1", "LRP2", "ITGA6", "TP63", "VGLL1", "CDH1"), pt.size = 0)

plot <- VlnPlot(cluster1, features = c("CDH2", "NOTCH1",  "ITGB4", "rna_TgGRA1", "ITGA2", "BACH2", "BCAT1", "LMCD1"), pt.size = 0) 

plot <- VlnPlot(cluster1, features = c("MKI67","TOP2A", "TEAD1", "TCF7L2", "GCM1", "FN1", "MMP2", "HLA-G", "HTRA4", "DIO2", "ITGA5","NOTUM"), pt.size = 0)

CTBmarkers <- c("PEG10","EGFR", "TFAP2C", "LYN", "TEAD4", "YAP1", "LRP2", "ITGA6", "TP63", "VGLL1", "CDH1", "ITGB4")
STBmarkers <- c("CGA", "PAPPA", "CGB1", "CGB3", "CGB5", "CGB7", "KISS1", "OVOL1", "SDC1", "CYP19A1", "ERVW-1", "ERVFRD-1")
iCTBmarkers <- c("CDH2", "NOTCH1", "ITGA2", "BACH2", "MKI67","TOP2A", "NRG1", "JUN", "FOS", "SMAD2", "SMAD3", "JDP2", "TEAD1", "TCF7L2", "GATA3") 
EVTmarkers <- c("GCM1", "ASCL2", "FN1", "MMP2", "MMP9", "HLA-G", "HTRA4", "DIO2", "ITGA5", "NOTUM", "ITGA1", "CSH1", "EPAS1")

features <- list("pre-/STB markers" = STBmarkers, 
              "CTB markers" = CTBmarkers, 
              "pre-EVT markers" = iCTBmarkers,
              "pre-/EVT markers" = EVTmarkers)

features <- c(STBmarkers, CTBmarkers, iCTBmarkers, EVTmarkers)

plot <- DotPlot(multiome, features = features, cluster.idents = TRUE) + 
  coord_flip() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position = "right")
plot
SaveFigure(plot, "DotPlot_vertical_samefeaturesastrophall_bysample_withclassifyers_5-16-24", height = 15, width = 4)
plot <- DimPlot(multiome, reduction = "umap")
SaveFigure(plot, "DotPlot_clusters5-15-24", height = 3, width = 15)

DefaultAssay(multiome) <- "RNA"
multiome.tree <- BuildClusterTree(object = multiome, features = features)
PlotClusterTree(multiome.tree, 
                type = "tidy", 
                show.node.label = FALSE, 
                cex = 1.5, 
                direction = "rightwards", 
                edge.width = 3,
                tip.color = c("#F6766D", "#A4A601", "#3CC07C", "#4EADF0", "#E76BF3"))


plot <- RidgePlot(multiome, features = c("CGA", "EGFR", "TEAD1"))
SaveFigure(plot, "RidgePlot_samefeaturesastrophall_5-14-24", height = 3, width = 8)

plot <- VlnPlot(gra1pos, features = c("rna_TgGRA1"), pt.size = 1)
SaveFigure(plot, "violin_TgGRA1_5-14-24", height = 6, width = 6, res = 300)


plotA <- FeaturePlot(multiome, features = c("BACH2", "rna_TgGRA1"), blend = TRUE, blend.threshold = .1, cols = c("lightgrey", "red", "blue"))
plotB <- FeaturePlot(multiome, features = c("NRG1", "rna_TgGRA1"), blend = TRUE, blend.threshold = .1, cols = c("lightgrey", "red", "blue"))
plotC <- FeaturePlot(multiome, features = c("TEAD1", "rna_TgGRA1"), blend = TRUE, blend.threshold = .5, cols = c("lightgrey", "red", "blue"))

SaveFigure(plotA, "FeaturePlot_blend_BACH2TgGRA1_5-15-24", height = 5, width = 16)

#name clusters

new_ids <- c("spontaneous TS-STB", "TS-CTB", "TS-EVT", "abberant TS-EVT", "proto TS-STB")
names(new_ids) <- levels(multiome)
multiome <- RenameIdents(object = multiome, new_ids)
DimPlot(multiome)


FeaturePlot(multiome, features = c("SP6"), reduction = "umap")

SaveObject(multiome, "multiome_working")
SaveObject(gra1pos, "gra1pos_working")
SaveObject(TGI_CTB, "TGI_CTB_working")
SaveObject(MI_CTB, "MI_CTB_working")
SaveObject(D4_EVT, "D4_EVT_working")
SaveObject(D8_EVT, "D8_EVT_working")

multiome <- ReadObject("multiome_working")

proptest <- sc_utils(multiome)
proptest <- permutation_test(proptest, cluster_identity = "named.clusters", sample_1 = "1", sample_2 = "2", sample_identity = "orig.ident", n_permutations = 10000)
plot <- permutation_plot(proptest, log2FD_threshold = log2(1.4142), order_clusters = FALSE) + 
  ggtitle("Mock-infected TS-CTB vs. Tg-infected TS-CTB") + 
  theme(plot.title = element_text(hjust = 0.5, size = 11, face = "bold"), 
        panel.grid.minor = element_line(size = 0.5), 
        panel.grid.major = element_line(size = 0.5),
        panel.border = element_rect(linewidth = 1, color = "black"),
        legend.position = "bottom",
        axis.text.x = element_text(face = "bold"),
        axis.text.y = element_text(face = "bold"),
        legend.box.background = element_rect(color = "black", linewidth = 1)) +
  scale_color_manual(values = c("red", "darkgrey"))
SaveFigure(plot, "permutation_test_MITSCTBvTGTSCTB", height = 3, width = 6)

#PSEUDOTIME ANALYSIS
#Monocle3 using script from trophall
cds <- as.cell_data_set(multiome)
colData(cds)
fData(cds)
rownames(fData(cds))[1:10]
fData(cds)$gene_short_name <- rownames(fData(cds))
counts(cds)
recreate.partition <- rep(1, length(cds@colData@rownames))

SaveObject(cds, "multiome_cds")
cds <- ReadObject("multiome_cds")

names(recreate.partition) <- cds@colData@rownames
recreate.partition <- as.factor(recreate.partition)
recreate.partition
cds@clusters$UMAP$partitions <- recreate.partition

list_cluster <- multiome@active.ident
cds@clusters$UMAP$clusters <- list_cluster

cds@int_colData@listData$reducedDims$UMAP <- multiome@reductions$umap@cell.embeddings

cluster.before.trajectory <- plot_cells(cds,
                                        color_cells_by = "cluster",
                                        label_cell_groups = FALSE,
                                        group_label_size = 5) +
  theme(legend.position = "right")

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

cds <- order_cells(cds, reduction_method = "UMAP", root_cells = colnames(cds[,clusters(cds) == 1])) # == "" needs to be set to your root

plot <- plot_cells(cds,
                   color_cells_by = "pseudotime", #change to "sample" or other col in colData(cds)
                   label_branch_points = FALSE,
                   label_roots = FALSE,
                   label_leaves = FALSE,
                   group_label_size = 10,
                   label_cell_groups = FALSE,
                   cell_size = 0.5,
                   show_trajectory_graph = FALSE,
                   trajectory_graph_segment_size = 2,
                   trajectory_graph_color = "black")
plot <- plot + theme(legend.position = "right")
SaveFigure(plot, "pseudotime_2-29-24", height = 8, width = 8)

pseudotime(cds)
cds$monocle3_pseudotime <- pseudotime(cds)
data.pseudo <- as.data.frame(colData(cds))
data.pseudo

ordered.pseudo <- ggplot(data.pseudo, aes(monocle3_pseudotime, reorder(ident, monocle3_pseudotime, median), fill = ident)) + geom_boxplot()
ordered.pseudo <- ordered.pseudo + xlab("pseudotime") + ylab("ident")
ordered.pseudo

#list of decidual related genes: NRG1, DDK1, BACH2, JDP2, IL6, IL8, BMP2, INHBA, TGFB2, LIF, INHBB

#trying to permanently rename identities

multiome <- RenameIdents(multiome, `0` = "spontaneous TS-STB", `1` = "TS-CTB", `2` = "TS-EVT", `3` = "abberant TS-EVT", `4` = "proto TS-STB")
multiome@meta.data$wsnn_res.0.06
multiome@meta.data$seurat_clusters

table(multiome$seurat_clusters)
table(multiome$wsnn_res.0.06)

#this worked thank god
multiome <- Rename_Clusters(multiome, new_idents = new_ids, meta_col_name = "named.clusters")
multiome$named.clusters
