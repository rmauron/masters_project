library(ArchR)
library(tidyverse)
library(ggpubr)

#addArchRGenome("hg38") #for human
addArchRGenome("mm10") #for mouse
library(BSgenome.Mmusculus.UCSC.mm10)
addArchRThreads(32) #rocommended 1/2 of total available cores
set.seed(1)

## In case you have to reload parts of the project
# projMulti1 <- loadArchRProject(path = "./Save-ProjMulti1", force = FALSE, showLogo = F)
# projMulti2 <- loadArchRProject(path = "./Save-ProjMulti2", force = FALSE, showLogo = F)
# projMulti3 <- loadArchRProject(path = "./Save-ProjMulti3", force = FALSE, showLogo = F)
# projMulti4 <- loadArchRProject(path = "./Save-ProjMulti4", force = FALSE, showLogo = F)


#when working on server, need to set the working directory aka the path to your input files:
#setwd("/Users/raphael.mauron/archr/archr_multiome/multiome_data")

#Get Input Fragment Files. Here, input the ATAC-seq data.
#insert as a vector c("sample1", "sample2", "sample3")
inputFiles <- getInputFiles(c("/Users/raphaelmauron/archr_rr3/data-raw/ATAC/CF1",
                              "/Users/raphaelmauron/archr_rr3/data-raw/ATAC/CF2",
                              "/Users/raphaelmauron/archr_rr3/data-raw/ATAC/CF7",
                              "/Users/raphaelmauron/archr_rr3/data-raw/ATAC/CG8",
                              "/Users/raphaelmauron/archr_rr3/data-raw/ATAC/CG9"))


#insert corresponding name c("name1", "name2", "name3")
names(inputFiles) <- c("CF1",
                       "CF2",
                       "CF7",
                       "CG8",
                       "CG9")


#check the inputs
inputFiles


ArrowFiles <- createArrowFiles(
    inputFiles = inputFiles, #contains only the ATAC data for now
    sampleNames = names(inputFiles), #names(atacFiles) #if you use mock dataset
    minTSS = 4,
    minFrags = 1000,
    addTileMat = TRUE,
    addGeneScoreMat = TRUE
)


projMulti1 <- ArchRProject(ArrowFiles = ArrowFiles)


#create an ArchR object containing the RNA data

#setwd("/Users/raphael.mauron/archr/archr_multiome/multiome_data") #path to data when working on server

seRNA <- import10xFeatureMatrix(
    input = c("/Users/raphaelmauron/archr_rr3/data-raw/ATAC/CF1/filtered_feature_bc_matrix.h5",
              "/Users/raphaelmauron/archr_rr3/data-raw/ATAC/CF2/filtered_feature_bc_matrix.h5",
              "/Users/raphaelmauron/archr_rr3/data-raw/ATAC/CF7/filtered_feature_bc_matrix.h5",
              "/Users/raphaelmauron/archr_rr3/data-raw/ATAC/CG8/filtered_feature_bc_matrix.h5",
              "/Users/raphaelmauron/archr_rr3/data-raw/ATAC/CG9/filtered_feature_bc_matrix.h5"),

    names = c("CF1",
              "CF2",
              "CF7",
              "CG8",
              "CG9"),

    strictMatch = TRUE
)


#Count per condition
condition_sample <- table(projMulti1@cellColData$Sample)
condition_sample <- as.data.frame(condition_sample)
condition_sample <- cbind(condition_sample, before_filter = c("before", "before", "before", "before", "before"))
condition_sample <- setnames(condition_sample, c("Sample", "nCells", "filter"))
condition_sample


p_sample <- ggplot(condition_sample,
                   aes(x = Sample, y = nCells)) +
    geom_bar(stat = "identity")+
    ggtitle("Number of cells per sample")
#p_sample


cat("Total number of cells among the samples before doublets removal is:", nCells(projMulti1))


p_frags <- plotGroups(
    projMulti1,
    groupBy = "Sample",
    colorBy = "cellColData",
    name = "nFrags",
    plotAs = "violin",
    alpha = 0.4,
    addBoxPlot = T,
    log = F
)
#p_frags

p_TSS <- plotGroups(
    projMulti1,
    groupBy = "Sample",
    colorBy = "cellColData",
    name = "TSSEnrichment",
    plotAs = "violin",
    alpha = 0.4,
    addBoxPlot = T
)
#p_TSS

p_multi_frags <- plotGroups(
    projMulti1,
    groupBy = "Sample",
    colorBy = "cellColData",
    name = "nMultiFrags",
    plotAs = "violin",
    alpha = 0.4,
    addBoxPlot = T
)
#p_multi_frags


plotPDF(p_frags, p_TSS, p_multi_frags, name = "QC-sample", addDOC = FALSE)
rm(p_frags, p_multi_frags, p_TSS)


#adding conditions Flight vs Ground as a CellColData
projMulti1 <- addCellColData(ArchRProj = projMulti1, data = paste0(gsub("0|1|2|3|4|5|6|7|8|9","", projMulti1@cellColData$Sample)), name = "Conditions", cells = getCellNames(projMulti1), force = TRUE)


p_frags <- plotGroups(
    projMulti1,
    groupBy = "Conditions",
    colorBy = "cellColData",
    name = "nFrags",
    plotAs = "violin",
    alpha = 0.4,
    addBoxPlot = T,
    log = F
)
#p_frags


#Count per condition
condition_table <- table(projMulti1@cellColData$Conditions)
condition_table <- as.data.frame(condition_table)
condition_table <- setnames(condition_table, c("Condition", "nCells"))
condition_table


p_cond <- ggplot(condition_table,
                 aes(x = Condition,
                     y = nCells)) +
    geom_bar(stat = "identity")+
    ggtitle("Number of cells per condition")
#p_cond
rm(condition_table)


plotPDF(p_sample, p_frags, p_cond, name = "QC-condition-sample", addDOC = FALSE)
rm(p_sample, p_frags, p_cond)


#save project since we will change project now
saveArchRProject(ArchRProj = projMulti1, outputDirectory = "Save-ProjMulti1", load = F)


#load the project from last saved project (if necessary)
#library(ArchR)
#projMulti1 <- loadArchRProject(path = "./Save-ProjMulti1", force = FALSE, showLogo = F)


length(which(getCellNames(projMulti1) %ni% colnames(seRNA)))


cellsToKeep <- which(getCellNames(projMulti1) %in% colnames(seRNA))
projMulti2 <- subsetArchRProject(
    ArchRProj = projMulti1,
    cells = getCellNames(projMulti1)[cellsToKeep],
    outputDirectory = "Save-ProjMulti2",
    force = TRUE)


#adding the RNA data to the project together with the ATAC data
projMulti2 <- addGeneExpressionMatrix(input = projMulti2, seRNA = seRNA, strictMatch = TRUE, force = TRUE)


projMulti2 <- addDoubletScores(projMulti2, force = TRUE)
projMulti2 <- filterDoublets(projMulti2)


cat("Total number of cells among the samples after doublets removal is:", nCells(projMulti2),". \n")

cat("A total of", nCells(projMulti1)-nCells(projMulti2), "doublets were removed.")


#Count per condition
condition_after_filter <- table(projMulti2@cellColData$Sample)
condition_after_filter <- as.data.frame(condition_after_filter)
condition_after_filter <- cbind(condition_after_filter, after_filter = c("after", "after", "after", "after", "after"))
condition_after_filter <- setnames(condition_after_filter, c("Sample", "nCells", "filter"))

filter_diff <- condition_sample[,"nCells"] - condition_after_filter[,"nCells"]
filter_df <- as.data.frame(cbind(Sample = c("CF1", "CF2", "CF7", "CG8", "CG9"),
                                 nCells = filter_diff,
                                 filter = c("diff", "diff", "diff", "diff", "diff")))

super_frame <- rbind(condition_sample, condition_after_filter)
super_frame$filter = factor(super_frame$filter, levels = c("before", "after"), ordered = TRUE)

p_before_after <- ggplot(super_frame,
                         aes(x=Sample,
                             y=nCells,
                             fill=filter)) +
    geom_bar(position = "dodge", stat="identity", width=0.9) +
    ggtitle("Sample before-after Quality Control") +
    theme(legend.position="right") +
    theme_minimal() +
    scale_fill_manual(values=c("#5ec962", "#3b528b"))

# Export dataframe as csv (for report)
write.csv(super_frame, "/Users/raphaelmauron/masters_project/doc/archr_multi/df_nCells_before_after.csv", row.names=T)

rm(condition_after_filter, filter_diff, filter_df, super_frame, condition_sample)


plotPDF(p_before_after, name = "Sample-before-after-QC", ArchRProj = projMulti2, addDOC = FALSE)
rm(p_before_after)


#LSI ATAC
projMulti2 <- addIterativeLSI(
    ArchRProj = projMulti2,
    clusterParams = list(
        resolution = 0.2,
        sampleCells = 10000,
        n.start = 10
    ),
    saveIterations = FALSE,
    useMatrix = "TileMatrix",
    depthCol = "nFrags",
    name = "LSI_ATAC",
    seed = 5,
    force = T
)


#LSI RNA
projMulti2 <- addIterativeLSI(
    ArchRProj = projMulti2,
    clusterParams = list(
        resolution = 0.2,
        sampleCells = 10000,
        n.start = 10
    ),
    saveIterations = FALSE,
    useMatrix = "GeneExpressionMatrix",
    depthCol = "Gex_nUMI",
    varFeatures = 2500,
    firstSelection = "variable",
    binarize = FALSE,
    scaleDims = FALSE,
    name = "LSI_RNA",
    seed = 3,
    force = T
)


#LSI from both ATAC and RNA
projMulti2 <- addCombinedDims(projMulti2, reducedDims = c("LSI_ATAC", "LSI_RNA"), name =  "LSI_Combined")


#generate UMAP from the reduced dimensions performed with LSI
projMulti2 <- addUMAP(projMulti2, reducedDims = "LSI_ATAC", name = "UMAP_ATAC", nNeighbors = 40, minDist = 0.1, force = TRUE)
projMulti2 <- addUMAP(projMulti2, reducedDims = "LSI_RNA", name = "UMAP_RNA", nNeighbors = 40, minDist = 0.1, force = TRUE)
projMulti2 <- addUMAP(projMulti2, reducedDims = "LSI_Combined", name = "UMAP_Combined", nNeighbors = 40, minDist = 0.1, force = TRUE)


projMulti2 <- addClusters(projMulti2, reducedDims = "LSI_ATAC", name = "Clusters_ATAC", resolution = 0.3, maxClusters = 20, force = TRUE)
projMulti2 <- addClusters(projMulti2, reducedDims = "LSI_RNA", name = "Clusters_RNA", resolution = 0.3, maxClusters = 20,force = TRUE)
projMulti2 <- addClusters(projMulti2, reducedDims = "LSI_Combined", name = "Clusters_Combined", resolution = 0.3, maxClusters = 20, force = TRUE)


p1 <- plotEmbedding(projMulti2, name = "Clusters_ATAC", embedding = "UMAP_ATAC", size = 1, labelAsFactors=F, labelMeans=F)
p2 <- plotEmbedding(projMulti2, name = "Clusters_RNA", embedding = "UMAP_RNA", size = 1, labelAsFactors=F, labelMeans=F)
p3 <- plotEmbedding(projMulti2, name = "Clusters_Combined", embedding = "UMAP_Combined", size = 1, labelAsFactors=F, labelMeans=F)

p <- lapply(list(p1,p2,p3), function(x){
    x + guides(color = "none", fill = "none") +
        theme_ArchR(baseSize = 6.5) +
        theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm")) +
        theme(
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank()
        )
})

#do.call(cowplot::plot_grid, c(list(ncol = 3),p))
rm(p)


plotPDF(p1, p2, p3, name = "UMAP-scATAC-scRNA-Combined", addDOC = FALSE)
rm(p1, p2, p3)


cM_atac_rna <- confusionMatrix(paste0(projMulti2$Clusters_ATAC), paste0(projMulti2$Clusters_RNA))
cM_atac_rna <- cM_atac_rna / Matrix::rowSums(cM_atac_rna)

library(pheatmap)
p_atac_rna <- pheatmap::pheatmap(
    mat = as.matrix(cM_atac_rna),
    color = paletteContinuous("whiteBlue"),
    border_color = "black"
)


plotPDF(p_atac_rna, name = "Confusion-Matrix-scATAC-scRNA", addDOC = FALSE)
rm(p_atac_rna)


#by sample and by clusters
p1 <- plotEmbedding(ArchRProj = projMulti2, colorBy = "cellColData", name = "Sample", embedding = "UMAP_Combined")

p2 <- plotEmbedding(ArchRProj = projMulti2, colorBy = "cellColData", name = "Clusters_Combined", embedding = "UMAP_Combined")

p3 <- plotEmbedding(ArchRProj = projMulti2, colorBy = "cellColData", name = "Conditions", embedding = "UMAP_Combined", alpha = 0.5, pal = c("CF" = "white","CG" = "green"))

p4 <- plotEmbedding(ArchRProj = projMulti2, colorBy = "cellColData", name = "Conditions", embedding = "UMAP_Combined", alpha = 0.5, pal = c("CF" = "blue","CG" = "white"))

p5 <- plotEmbedding(ArchRProj = projMulti2, colorBy = "cellColData", name = "Conditions", embedding = "UMAP_Combined", alpha = 0.35, pal = c("CF" = "blue","CG" = "green"))

# p1
# p2
# p3
# p4
# p5


plotPDF(p1, p2, p3, p4, p5, name = "UMAP-Combined-sample-clusters", addDOC = FALSE)
rm(p1, p2, p3, p4)


#by sample and by clusters
p1 <- plotEmbedding(ArchRProj = projMulti2, colorBy = "cellColData", name = "Sample", embedding = "UMAP_ATAC")

p2 <- plotEmbedding(ArchRProj = projMulti2, colorBy = "cellColData", name = "Clusters_ATAC", embedding = "UMAP_ATAC")

p3 <- plotEmbedding(ArchRProj = projMulti2, colorBy = "cellColData", name = "Conditions", embedding = "UMAP_ATAC", alpha = 0.5, pal = c("CF" = "white","CG" = "green"))

p4 <- plotEmbedding(ArchRProj = projMulti2, colorBy = "cellColData", name = "Conditions", embedding = "UMAP_ATAC", alpha = 0.5, pal = c("CF" = "blue","CG" = "white"))

p5 <- plotEmbedding(ArchRProj = projMulti2, colorBy = "cellColData", name = "Conditions", embedding = "UMAP_ATAC", alpha = 0.35, pal = c("CF" = "blue","CG" = "green"))

# p1
# p2
# p3
# p4
# p5


plotPDF(p1, p2, p3, p4, p5, name = "UMAP-ATAC-sample-clusters", addDOC = FALSE)
rm(p1, p2, p3, p4)


#by sample and by clusters
p1 <- plotEmbedding(ArchRProj = projMulti2, colorBy = "cellColData", name = "Sample", embedding = "UMAP_RNA")

p2 <- plotEmbedding(ArchRProj = projMulti2, colorBy = "cellColData", name = "Clusters_RNA", embedding = "UMAP_RNA")

p3 <- plotEmbedding(ArchRProj = projMulti2, colorBy = "cellColData", name = "Conditions", embedding = "UMAP_RNA", alpha = 0.5, pal = c("CF" = "white","CG" = "green"))

p4 <- plotEmbedding(ArchRProj = projMulti2, colorBy = "cellColData", name = "Conditions", embedding = "UMAP_RNA", alpha = 0.5, pal = c("CF" = "blue","CG" = "white"))

p5 <- plotEmbedding(ArchRProj = projMulti2, colorBy = "cellColData", name = "Conditions", embedding = "UMAP_RNA", alpha = 0.35, pal = c("CF" = "blue","CG" = "green"))

# p1
# p2
# p3
# p4
# p5


plotPDF(p1, p2, p3, p4, p5, name = "UMAP-RNA-sample-clusters", addDOC = FALSE)
rm(p1, p2, p3, p4)


#Harmony Combined
projMulti2 <- addHarmony(
    ArchRProj = projMulti2,
    reducedDims = "LSI_Combined",
    name = "Harmony",
    groupBy = "Sample",
    theta = 0.001,
    lambda = 1,
    force = T
)

#Harmony ATAC
projMulti2 <- addHarmony(
    ArchRProj = projMulti2,
    reducedDims = "LSI_ATAC",
    name = "Harmony_ATAC",
    groupBy = "Sample",
    theta = 0.001,
    lambda = 1,
    force = T
)

#Harmony RNA
projMulti2 <- addHarmony(
    ArchRProj = projMulti2,
    reducedDims = "LSI_RNA",
    name = "Harmony_RNA",
    groupBy = "Sample",
    force = T
)


#dimensionality reduction (Harmony, not LSI) COMBINED
projMulti2 <- addUMAP(
    ArchRProj = projMulti2,
    reducedDims = "Harmony",
    name = "UMAPHarmony",
    nNeighbors = 40,
    minDist = 0.1,
    metric = "cosine",
    force = T
)

#dimensionality reduction (Harmony, not LSI) ATAC
projMulti2 <- addUMAP(
    ArchRProj = projMulti2,
    reducedDims = "Harmony_ATAC",
    name = "UMAPHarmony_ATAC",
    nNeighbors = 40,
    minDist = 0.1,
    metric = "cosine",
    force = T
)

#dimensionality reduction (Harmony, not LSI) RNA
projMulti2 <- addUMAP(
    ArchRProj = projMulti2,
    reducedDims = "Harmony_RNA",
    name = "UMAPHarmony_RNA",
    nNeighbors = 40,
    minDist = 0.1,
    metric = "cosine",
    force = T
)


p1 <- plotEmbedding(ArchRProj = projMulti2,
                    colorBy = "cellColData",
                    name = "Sample",
                    embedding = "UMAPHarmony")

p2 <- plotEmbedding(ArchRProj = projMulti2,
                    colorBy = "cellColData",
                    name = "Clusters_Combined",
                    embedding = "UMAPHarmony")

p3 <- plotEmbedding(ArchRProj = projMulti2,
                    colorBy = "cellColData",
                    name = "Conditions",
                    embedding = "UMAPHarmony",
                    alpha = 0.35, pal = c("CF" = "blue","CG" = "green")
)

p4 <- plotEmbedding(ArchRProj = projMulti2,
                    colorBy = "cellColData",
                    name = "Conditions",
                    embedding = "UMAPHarmony",
                    alpha = 0.35, pal = c("CF" = "white","CG" = "green")
)

p5 <- plotEmbedding(ArchRProj = projMulti2,
                    colorBy = "cellColData",
                    name = "Conditions",
                    embedding = "UMAPHarmony",
                    alpha = 0.35, pal = c("CF" = "blue","CG" = "white")
)
# p1
# p2
# p3
# p4
# p5


plotPDF(p1, p2, p3, p4, p5, name = "UMAP-Harmony-Combined-sample-clusters", addDOC = FALSE)
rm(p1, p2, p3, p4, p5)



p1 <- plotEmbedding(ArchRProj = projMulti2,
                    colorBy = "cellColData",
                    name = "Sample",
                    embedding = "UMAPHarmony_ATAC")

p2 <- plotEmbedding(ArchRProj = projMulti2,
                    colorBy = "cellColData",
                    name = "Clusters_ATAC",
                    embedding = "UMAPHarmony_ATAC")

p3 <- plotEmbedding(ArchRProj = projMulti2,
                    colorBy = "cellColData",
                    name = "Conditions",
                    embedding = "UMAPHarmony_ATAC",
                    alpha = 0.35, pal = c("CF" = "blue","CG" = "green")
)
p4 <- plotEmbedding(ArchRProj = projMulti2,
                    colorBy = "cellColData",
                    name = "Conditions",
                    embedding = "UMAPHarmony_ATAC",
                    alpha = 0.35, pal = c("CF" = "white","CG" = "green")
)
p5 <- plotEmbedding(ArchRProj = projMulti2,
                    colorBy = "cellColData",
                    name = "Conditions",
                    embedding = "UMAPHarmony_ATAC",
                    alpha = 0.35, pal = c("CF" = "blue","CG" = "white")
)
# p1
# p2
# p3
# p4
# p5


plotPDF(p1, p2, p3, p4, p5, name = "UMAP-Harmony-ATAC-sample-clusters", addDOC = FALSE)
rm(p1, p2, p3, p4, p5)


p1 <- plotEmbedding(ArchRProj = projMulti2,
                    colorBy = "cellColData",
                    name = "Sample",
                    embedding = "UMAPHarmony_RNA")

p2 <- plotEmbedding(ArchRProj = projMulti2,
                    colorBy = "cellColData",
                    name = "Clusters_RNA",
                    embedding = "UMAPHarmony_RNA")

p3 <- plotEmbedding(ArchRProj = projMulti2,
                    colorBy = "cellColData",
                    name = "Conditions",
                    embedding = "UMAPHarmony_RNA",
                    alpha = 0.35, pal = c("CF" = "blue","CG" = "green")
)

p4 <- plotEmbedding(ArchRProj = projMulti2,
                    colorBy = "cellColData",
                    name = "Conditions",
                    embedding = "UMAPHarmony_RNA",
                    alpha = 0.35, pal = c("CF" = "white","CG" = "green")
)

p5 <- plotEmbedding(ArchRProj = projMulti2,
                    colorBy = "cellColData",
                    name = "Conditions",
                    embedding = "UMAPHarmony_RNA",
                    alpha = 0.35, pal = c("CF" = "blue","CG" = "white")
)
# p1
# p2
# p3
# p4
# p5


plotPDF(p1, p2, p3, p4, p5, name = "UMAP-Harmony-RNA-sample-clusters", addDOC = FALSE)
rm(p1, p2, p3, p4, p5)


#include your marker gene list here
my_markers <- c(
    "Wfs1", # RR3 markers
    "Dkk3",
    "Prox1",
    "Egfr",
    "Ptn",
    "Efna5",
    "Vegfa",
    "Zic1",
    "Zic2",
    "Atoh1",
    "Sox2",
    "Nr4a2",
    "Pparg",
    "Rxra",
    "Nr2f6",
    "Hic1",
    "Sox6",
    "GBA", #Parkinson markers
    "LRRK2",
    "PINK1",
    "UCHL1",
    "VPS35",
    "ApoE"
)
my_markers <- unique(my_markers)


#bias correction
se <- getMarkerFeatures(ArchRProj = projMulti2,
                        groupBy = "Clusters_ATAC", #"Clusters_Combined",
                        bias = c("TSSEnrichment", #snATAC-seq data quality
                                 "log10(nFrags)", #read depth for snATAC
                                 "log10(Gex_nUMI)") #read depth for snRNA
)

#heatmap of all the marker genes
heatmap_gex <- plotMarkerHeatmap(
    seMarker = se,
    cutOff = "FDR <= 0.01 & Log2FC >= 0.5",
    transpose = T,
    labelMarkers = my_markers #specify the markers you want a label on
)

#draw(heatmap_gex, heatmap_legend_side = "bot", annotation_legend_side = "bot")


plotPDF(heatmap_gex,
        name = "GeneExpression-Marker-Heatmap",
        width = 12,
        height = 7,
        addDOC = FALSE)
rm(heatmap_gex)


# list of marker genes per cluster (ranked in comparison to the markers from the other markers)
markerList <- getMarkers(se, cutOff = "FDR <= 0.01 & Log2FC >= 1.25")

i <- 1
top <- 20 #select the top X markers per cluster
for (cluster in 1:length(markerList@listData)) {
    if (nrow(markerList@listData[[cluster]]["name"]) >= top) {
        if (i == 1) {
            my_row <- head(markerList@listData[[cluster]]["name"],top)
            df_top <- as.data.frame(my_row)
            colnames(df_top)[1] <- gsub(" ", "",paste("Cluster", cluster))
        } else {
            df_top <- cbind(df_top, head(markerList@listData[[cluster]]["name"],top))
            colnames(df_top)[i] <- gsub(" ", "",paste("Cluster", cluster))
        }
    } else {
        next
    }
    i <- i + 1
}
rm(cluster, i)

# ROWNAMES
row_list <- c()
for (row in 1:top) {
    row_list <- append(row_list, gsub(" ", "",paste("MarkerGene", row)))
}

rownames(df_top) <- row_list
df_top
rm(row, row_list, top)


# list of marker genes per cluster (ranked in comparison to the markers from the other markers)
markerList <- getMarkers(se, cutOff = "FDR <= 0.01 & Log2FC >= 1.25")

i <- 1
top <- 3
for (cluster in 1:length(markerList@listData)) {
    if (nrow(markerList@listData[[cluster]]["name"]) >= top) {
        if (i == 1) {
            my_row <- head(markerList@listData[[cluster]]["name"],top)
            df_top3 <- as.data.frame(my_row)
            colnames(df_top3)[1] <- gsub(" ", "",paste("Cluster", cluster))
        } else {
            df_top3 <- cbind(df_top3, head(markerList@listData[[cluster]]["name"],top))
            colnames(df_top3)[i] <- gsub(" ", "",paste("Cluster", cluster))
        }
    } else {
        next
    }
    i <- i + 1
}
rm(cluster, i)

# ROWNAMES
row_list <- c()
for (row in 1:top) {
    row_list <- append(row_list, gsub(" ", "",paste("MarkerGene", row)))
}
rownames(df_top3) <- row_list
df_top3


# take every element of the small list and create 1 list
list_top3_df <- as.list(df_top3)
list_top3 <- c()

for (c in 1:length(list_top3_df)) {
    for (i in 1:top) {
        list_top3 <- append(list_top3, list_top3_df[[c]][i])
        list_top3 <- unique(list_top3)
    }
}
list_top3

# remove variable
rm(row, row_list, top, list_top3_df, c, i)


# Export dataframe as csv (for report)
write.csv(t(df_top3), "/Users/raphaelmauron/masters_project/doc/archr_multi/df_top3_markergenes.csv", row.names=T)


p <- plotEmbedding(
    ArchRProj = projMulti2,
    colorBy = "GeneExpressionMatrix",
    name = c(my_markers, list_top3),
    embedding = "UMAPHarmony_ATAC", #can select the dimensionality reduction here
    quantCut = c(0.01, 0.95),
    imputeWeights = NULL
)


# one embedding per selected genes
p <- lapply(p, function(x){
    x + guides(color = FALSE, fill = FALSE) +
        theme_ArchR(baseSize = 6.5) +
        theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
        theme(
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank()
        )
})

#uncomment to plot here, but might overload RStudio
#do.call(cowplot::plot_grid, c(list(ncol = 1),p))


plotPDF(plotList = p,
        name = "Plot-UMAP-Marker-Genes-WO-Imputation.pdf",
        ArchRProj = projMulti2,
        addDOC = FALSE,
        width = 5,
        height = 5)
rm(p)


projMulti2 <- addImputeWeights(ArchRProj = projMulti2,
                               reducedDims = "LSI_ATAC" #"LSI_Combined"
)


p <- plotEmbedding(
    ArchRProj = projMulti2,
    colorBy = "GeneExpressionMatrix",
    name = c(my_markers, list_top3),
    embedding = "UMAPHarmony_ATAC", #can select the dimensionality reduction here
    imputeWeights = getImputeWeights(projMulti2)
)


#Rearrange for grid plotting
p <- lapply(p, function(x){
    x + guides(color = FALSE, fill = FALSE) +
        theme_ArchR(baseSize = 6.5) +
        theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
        theme(
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank()
        )
})

#uncomment to plot, might overload the session
#do.call(cowplot::plot_grid, c(list(ncol = 3),p))


plotPDF(plotList = p,
        name = "Plot-UMAP-Marker-Genes-W-Imputation.pdf",
        ArchRProj = projMulti2,
        addDOC = FALSE,
        width = 5,
        height = 5)
rm(p)


p <- plotBrowserTrack(
    ArchRProj = projMulti2,
    groupBy = "Clusters_ATAC", #Choose the clustering
    geneSymbol = c(my_markers, list_top3),
    upstream = 50000,
    downstream = 50000
)


plotPDF(plotList = p,
        name = "Plot-Tracks-Marker-Genes.pdf",
        ArchRProj = projMulti2,
        addDOC = FALSE,
        width = 5,
        height = 5)
rm(p)


#save project since we will change project now
saveArchRProject(ArchRProj = projMulti2, outputDirectory = "Save-ProjMulti2", load = F)


#load the project from last saved project (if necessary)
#library(ArchR)
#projMulti2 <- loadArchRProject(path = "./Save-ProjMulti2", force = FALSE, showLogo = F)


#loading the database file in seurat format (.Rds)

#setwd("/Users/raphaelmauron/archr_rr3/data-raw/ref_data/mousebrain_tax4_subsampled.Rds")

se_RNA_mousebrain_data <- readRDS(file = "/Users/raphaelmauron/archr_rr3/data-raw/ref_data/mousebrain_tax4_subsampled.Rds")
#se_RNA_mousebrain_data


#check requirement for seurat object before integration
#1
"RNA" %in% names(se_RNA_mousebrain_data@assays)

#2 all the cell types found our sample:
unique(se_RNA_mousebrain_data@meta.data$TaxonomyRank4)


#unconstrained method (8.1.1 in manual)
projMulti3 <- addGeneIntegrationMatrix(
    ArchRProj = projMulti2,
    useMatrix = "GeneExpressionMatrix",
    matrixName = "GeneIntegrationMatrix_db_1",
    reducedDims = "Harmony", #"LSI_Combined",
    seRNA = se_RNA_mousebrain_data,
    addToArrow = FALSE,
    groupRNA = "TaxonomyRank1",
    nameCell = "predictedCell_Un_1", #store cell ID from matched RNA expr
    nameGroup = "predictedGroup_Un_1", #store the group ID
    nameScore = "predictedScore_Un_1" #store the cross-platform integration score
)

p1 <- plotEmbedding(
    projMulti3,
    colorBy = "cellColData",
    name = "predictedGroup_Un_1",
    embedding = "UMAPHarmony",
    pal = paletteDiscrete(values = unique(projMulti3@cellColData@listData[["predictedGroup_Un_1"]]))
)

p1_1 <- plotEmbedding(
    projMulti3,
    colorBy = "cellColData",
    name = "predictedGroup_Un_1",
    embedding = "UMAPHarmony_ATAC",
    pal = paletteDiscrete(values = unique(projMulti3@cellColData@listData[["predictedGroup_Un_1"]]), set = "stallion2")
)
#p1
#p1_1


projMulti3 <- addGeneIntegrationMatrix(
    ArchRProj = projMulti2,
    useMatrix = "GeneExpressionMatrix",
    matrixName = "GeneIntegrationMatrix_db_2",
    reducedDims = "Harmony", #"LSI_Combined"
    seRNA = se_RNA_mousebrain_data,
    addToArrow = FALSE,
    groupRNA = "TaxonomyRank2",
    nameCell = "predictedCell_Un_2",
    nameGroup = "predictedGroup_Un_2",
    nameScore = "predictedScore_Un_2"
)

p2 <- plotEmbedding(
    projMulti3,
    colorBy = "cellColData",
    name = "predictedGroup_Un_2",
    embedding = "UMAPHarmony",
    pal = paletteDiscrete(values = unique(projMulti3@cellColData@listData[["predictedGroup_Un_2"]]))
)

p2_1 <- plotEmbedding(
    projMulti3,
    colorBy = "cellColData",
    name = "predictedGroup_Un_2",
    embedding = "UMAPHarmony_ATAC",
    pal = paletteDiscrete(values = unique(projMulti3@cellColData@listData[["predictedGroup_Un_2"]]))
)
# p2
# p2_1


projMulti3 <- addGeneIntegrationMatrix(
    ArchRProj = projMulti2,
    useMatrix = "GeneExpressionMatrix",
    matrixName = "GeneIntegrationMatrix_db_4",
    reducedDims = "Harmony", #"LSI_Combined"
    seRNA = se_RNA_mousebrain_data,
    addToArrow = FALSE, #set TRUE to the Taxonomy you want to add to the ArrowFiles
    groupRNA = "TaxonomyRank4",
    nameCell = "predictedCell_Un_4",
    nameGroup = "predictedGroup_Un_4",
    nameScore = "predictedScore_Un_4"
)


p4 <- plotEmbedding(
    projMulti3,
    colorBy = "cellColData",
    name = "predictedGroup_Un_4",
    embedding = "UMAPHarmony",
    pal = paletteDiscrete(values = unique(projMulti3@cellColData@listData[["predictedGroup_Un_4"]]), set = "stallion")

)

p4_1 <- plotEmbedding(
    projMulti3,
    colorBy = "cellColData",
    name = "predictedGroup_Un_4",
    embedding = "UMAPHarmony_ATAC",
    pal = paletteDiscrete(values = unique(projMulti3@cellColData@listData[["predictedGroup_Un_4"]]), set = "stallion")

)
# p4
# p4_1


projMulti3 <- addGeneIntegrationMatrix(
    ArchRProj = projMulti2,
    useMatrix = "GeneExpressionMatrix",
    matrixName = "GeneIntegrationMatrix_db_3",
    reducedDims = "Harmony_ATAC", #"LSI_Combined"
    seRNA = se_RNA_mousebrain_data,
    addToArrow = TRUE,
    force = TRUE, #to allow rerun
    groupRNA = "TaxonomyRank3",
    nameCell = "predictedCell_Un_3",
    nameGroup = "predictedGroup_Un_3",
    nameScore = "predictedScore_Un_3"
)

p3 <- plotEmbedding(
    projMulti3,
    colorBy = "cellColData",
    name = "predictedGroup_Un_3",
    embedding = "UMAPHarmony",
    pal = paletteDiscrete(values = unique(projMulti3@cellColData@listData[["predictedGroup_Un_3"]]))
)

p3_1 <- plotEmbedding(
    projMulti3,
    colorBy = "cellColData",
    name = "predictedGroup_Un_3",
    embedding = "UMAPHarmony_ATAC",
    pal = paletteDiscrete(values = unique(projMulti3@cellColData@listData[["predictedGroup_Un_3"]])),
    #legendSize = 50,
    labelSize = 5,
    colorTitle = ""

) +
    theme(legend.text = element_text(size = 14)) +
    guides(color = guide_legend(ncol = 2, override.aes = list(size = 4)) )

# p3_1
ggsave("tax3.png", path = "./Save-ProjMulti3/Plots", width = 10, height = 12)


plotPDF(p1, p1_1, p2, p2_1, p3, p3_1, p4, p4_1,
        name = "Plot-HarmonyUMAP-Cell-Annotation3.pdf",
        ArchRProj = projMulti3,
        addDOC = FALSE,
        width = 5,
        height = 5)
rm(p1, p1_1, p2, p2_1, p3, p3_1, p4, p4_1)


#confusion matrix between our scATAC-seq clusters and the predictedGroup
# cM <- confusionMatrix(projMulti3$Clusters_ATAC, projMulti3$predictedGroup_Un_4)
# labelOld <- rownames(cM)
# labelOld


#identify the cell type from predictedGroup which best defines that cluster
# labelNew <- colnames(cM)[apply(cM, 1, which.max)]
# labelNew


#rename the clusters from taxonomy level 4 with the celltypes printed above
# remapClust <- c(
#     "Astrocytes" = "Astrocytes",
#     "Cholinergic and monoaminergic neurons" = "Cholinergic and monoaminergic neurons",
#     "Choroid epithelial cells" = "Immature neural",
#     "Dentate gyrus granule neurons" = "Neurons",
#     "Dentate gyrus radial glia−like cells" = "Radial glia−like cells",
#     "Di− and mesencephalon excitatory neurons" = "Di- and mesencephalon neurons",
#     "Di− and mesencephalon inhibitory neurons" = "Di- and mesencephalon neurons",
#     "Ependymal cells" = "Astrocytes",
#     "Glutamatergic neuroblasts" = "Immature neural",
#     "Microglia" = "Microglia",
#     "Non−glutamatergic neuroblasts" = "Immature neural",
#     "Oligodendrocyte precursor cells" = "Oligodendrocyte precursor cells"            ,
#     "Oligodendrocytes" = "Oligodendrocytes",
#     "Peptidergic neurons" = "Neurons",
#     "Perivascular macrophages" = "Microglia",
#     "Subventricular zone radial glia−like cells" = "Radial glia−like cells",
#     "Telencephalon inhibitory interneurons" = "Telencephalon interneurons",
#     "Telencephalon projecting excitatory neurons" = "Telencephalon projecting neurons",
#     "Telencephalon projecting inhibitory neurons" = "Telencephalon projecting neurons",
#     "Vascular and leptomeningeal cells" = "Vascular cells",
#     "Vascular smooth muscle cells" = "Vascular cells"
# )
#
# remapClust <- remapClust[names(remapClust) %in% labelNew]



# labelNew2 <- mapLabels(labelNew, oldLabels = names(remapClust), newLabels = remapClust)
# labelNew2



#can use the mapLabels() function again to create new cluster labels in cellColData
# projMulti3$Clusters_ATAC_2 <- mapLabels(projMulti3$Clusters_ATAC, newLabels = labelNew2, oldLabels = labelOld)



# p1 <- plotEmbedding(projMulti3, colorBy = "cellColData", name = "Clusters_ATAC_2", embedding = "UMAPHarmony_ATAC")
# p1




# plotPDF(p1,
#         name = "Annotated_UMAP_HarmonyATAC.pdf",
#         ArchRProj = projMulti3,
#         addDOC = FALSE)
#rm(p1)



#new category which specifies the condition and the cell type
# Condition_CellType
projMulti3 <- addCellColData(
    ArchRProj = projMulti3,
    data = gsub(" ", "_", paste0(projMulti3@cellColData@listData[["predictedGroup_Un_3"]], "_x_", projMulti3@cellColData$Sample)),
    name = "CellType_Sample",
    cells = getCellNames(projMulti3),
    force = TRUE)

projMulti3 <- addCellColData(
    ArchRProj = projMulti3,
    data = gsub(" ", "_", paste0(projMulti3@cellColData@listData[["predictedGroup_Un_3"]], "_x_", projMulti3@cellColData$Conditions)),
    name = "CellType_Condition",
    cells = getCellNames(projMulti3),
    force = TRUE)

#return the unique tags from this new category
unique(projMulti3@cellColData@listData[["CellType_Condition"]])


df_cond_cell <- as.data.frame(table(projMulti3@cellColData@listData[["CellType_Condition"]]))
df_cond_cell <- setnames(df_cond_cell, c("CellType_Condition", "nCells"))
df_cond_cell$CellType_Condition <- as.character(df_cond_cell$CellType_Condition)

df_cond_cell <- df_cond_cell %>%
    mutate(Conditions = case_when(
        endsWith(CellType_Condition, "CF") ~ "CF",
        endsWith(CellType_Condition, "CG") ~ "CG"))

df_cond_cell


p_cond <- ggplot(df_cond_cell,
                 aes(x = CellType_Condition,
                     y = nCells,
                     fill = Conditions)) + #new Sample
    geom_bar(stat = "identity")+
    ggtitle("Number of cells per condition") +
    theme(axis.text.x = element_text(angle = 60, hjust=1))+
    theme(legend.position="bottom")
rm(df_cond_cell)


df <- as.data.frame(table(projMulti3@cellColData@listData[["CellType_Sample"]]))
df <- setnames(df, c("CellType_Sample", "nCells"))
df$CellType_Sample <- as.character(df$CellType_Sample)

sample_list <- c()
for (sample in df$CellType_Sample) {
    sample_list <- append(sample_list, sub("[0-9].*", "", sample))
}

df <- df %>% add_column(Samples = sample_list)

df <- df %>%
    mutate(Condition = case_when(
        endsWith(CellType_Sample, "CF1") ~ "CF1",
        endsWith(CellType_Sample, "CF2") ~ "CF2",
        endsWith(CellType_Sample, "CF7") ~ "CF7",
        endsWith(CellType_Sample, "CG8") ~ "CG8",
        endsWith(CellType_Sample, "CG9") ~ "CG9"
    ))


rm(sample, sample_list)


p_sample <- ggplot(df,
                   aes(x = Samples,
                       y = nCells,
                       fill = Condition)) + #new Sample
    geom_bar(position = "stack", stat = "identity")+
    ggtitle("Number of cells per condition and per sample") +
    theme(axis.text.x = element_text(angle = 60, hjust=1)) +
    theme(legend.position="bottom") +
    theme_minimal() +
    scale_fill_viridis_d()
rm(df)


p_comb <- ggarrange(p_cond, p_sample, ncol = 2,nrow = 1)


plotPDF(p_cond, p_sample, p_comb,
        name = "Plot-nCells_Condition_CellType.pdf",
        ArchRProj = projMulti3,
        addDOC = FALSE,
        width = 25,
        height = 10)
rm(p_cond, p_sample, p_comb)


#here want to perform metrics of each clusters
metrics <- as.data.frame(table(projMulti3@cellColData@listData[["predictedGroup_Un_3"]]))
metrics <- as.data.frame(t(metrics))

rownames(metrics)[1] <- "Cell type"
rownames(metrics)[2] <- "Total cells"

# renames columns
for (i in 1:length(metrics)) {
    colnames(metrics)[i] <- paste0("Cluster", i)
}

# Add cells from CF
my_table <- table(projMulti3@cellColData@listData[["CellType_Condition"]])
matches <- my_table[names(my_table) %in% grep("CF", sort(unique(projMulti3@cellColData@listData[["CellType_Condition"]])), value = T)]
metrics <- rbind(metrics, t(as.data.frame(matches)[2])[1,])
rownames(metrics)[3] <- "CF cells"

# Add cells from CG
my_table <- table(projMulti3@cellColData@listData[["CellType_Condition"]])
matches <- my_table[names(my_table) %in% grep("CG", sort(unique(projMulti3@cellColData@listData[["CellType_Condition"]])), value = T)]
metrics <- rbind(metrics, t(as.data.frame(matches)[2])[1,])
rownames(metrics)[4] <- "CG cells"

# Add ratio CF/CG
metrics <- rbind(metrics,  formatC(as.integer(metrics[3,])/as.numeric(metrics[4,]), digits = 3, format = "f"))
rownames(metrics)[5] <- "Ratio CF/CG"

# Add cells from CF1
my_table <- table(projMulti3@cellColData@listData[["CellType_Sample"]])
matches <- my_table[names(my_table) %in% grep("CF1", sort(unique(projMulti3@cellColData@listData[["CellType_Sample"]])), value = T)]
metrics <- rbind(metrics, t(as.data.frame(matches)[2])[1,])
rownames(metrics)[6] <- "CF1 cells"

# Add cells from CF2
my_table <- table(projMulti3@cellColData@listData[["CellType_Sample"]])
matches <- my_table[names(my_table) %in% grep("CF2", sort(unique(projMulti3@cellColData@listData[["CellType_Sample"]])), value = T)]
metrics <- rbind(metrics, t(as.data.frame(matches)[2])[1,])
rownames(metrics)[7] <- "CF2 cells"

# Add cells from CF7
my_table <- table(projMulti3@cellColData@listData[["CellType_Sample"]])
matches <- my_table[names(my_table) %in% grep("CF7", sort(unique(projMulti3@cellColData@listData[["CellType_Sample"]])), value = T)]
metrics <- rbind(metrics, t(as.data.frame(matches)[2])[1,])
rownames(metrics)[8] <- "CF7 cells"

# Add cells from CG8
my_table <- table(projMulti3@cellColData@listData[["CellType_Sample"]])
matches <- my_table[names(my_table) %in% grep("CG8", sort(unique(projMulti3@cellColData@listData[["CellType_Sample"]])), value = T)]
metrics <- rbind(metrics, t(as.data.frame(matches)[2])[1,])
rownames(metrics)[9] <- "CG8 cells"

# Add cells from CG9
my_table <- table(projMulti3@cellColData@listData[["CellType_Sample"]])
matches <- my_table[names(my_table) %in% grep("CG9", sort(unique(projMulti3@cellColData@listData[["CellType_Sample"]])), value = T)]
metrics <- rbind(metrics, t(as.data.frame(matches)[2])[1,])
rownames(metrics)[10] <- "CG9 cells"

# Add ratio CF1/total
metrics <- rbind(metrics,  formatC(as.integer(metrics[6,])/as.numeric(metrics[2,]), digits = 3, format = "f"))
rownames(metrics)[11] <- "Ratio CF1/total"

# Add ratio CF2/total
metrics <- rbind(metrics,  formatC(as.integer(metrics[7,])/as.numeric(metrics[2,]), digits = 3, format = "f"))
rownames(metrics)[12] <- "Ratio CF2/total"

# Add ratio CF7/total
metrics <- rbind(metrics,  formatC(as.integer(metrics[8,])/as.numeric(metrics[2,]), digits = 3, format = "f"))
rownames(metrics)[13] <- "Ratio CF7/total"

# Add ratio CG8/total
metrics <- rbind(metrics,  formatC(as.integer(metrics[9,])/as.numeric(metrics[2,]), digits = 3, format = "f"))
rownames(metrics)[14] <- "Ratio CG8/total"

# Add ratio CG9/total
metrics <- rbind(metrics, formatC(as.integer(metrics[10,])/as.numeric(metrics[2,]), digits = 3, format = "f"))
rownames(metrics)[15] <- "Ratio CG9/total"

# print dataframe
metrics

# Export dataframe as csv (for report)
write.csv(metrics, "/Users/raphaelmauron/masters_project/doc/archr_multi/metrics.csv", row.names=T)

rm(i, metrics, my_table, matches)


saveArchRProject(ArchRProj = projMulti3, outputDirectory = "Save-ProjMulti3", load = F)


#load the project from last saved project (if necessary)
#library(ArchR)
#projMulti3 <- loadArchRProject(path = "./Save-ProjMulti3", force = FALSE, showLogo = TRUE)


projMulti4 <- projMulti3

#add coverage (point 9 in manual)
projMulti4 <- addGroupCoverages(
    ArchRProj = projMulti4,
    groupBy = "CellType_Condition", #"Clusters_ATAC" #in CellColData,
)

#get path to the python tool installed through miniconda (MACS2)
projMulti4 <- addReproduciblePeakSet(
    ArchRProj = projMulti4,
    groupBy = "CellType_Condition", #"Clusters_ATAC", #in CellColData,
    pathToMacs2 = "/usr/local/bin/macs2"#"/Users/raphaelmauron/Library/Python/3.9/bin/macs2"
)

#add peaks just found
projMulti4 <- addPeakMatrix(ArchRProj = projMulti4)

#peak2genes (point 15.3 in manual)
projMulti4 <- addPeak2GeneLinks(
    ArchRProj = projMulti4,
    reducedDims = "Harmony_ATAC", #"LSI_Combined",
    useMatrix = "GeneExpressionMatrix")

p2g <- getPeak2GeneLinks(ArchRProj = projMulti4)

p2g[[1]]


#check the available matrices in the project. Should get "PeakMatrix"
getAvailableMatrices(projMulti4)


#Our scRNA labels to remind ourselves of the cell types that we are working with in projMulti4
table(projMulti4@cellColData@listData[["CellType_Condition"]])


# Create list used for 2 parameters needed below.
useGroups_list <- grep("CF", sort(unique(projMulti4@cellColData@listData[["CellType_Condition"]])), value = T)
bgdGroups_list <- grep("CG", sort(unique(projMulti4@cellColData@listData[["CellType_Condition"]])), value = T)


#pairwise testing for each cluster between conditions
i <- 1
marker_significant <- list() #list of dataframes of marker_pairwise_significant
markerTest <- list() #list of dataframes of marker_pairwise used for Volcano

for (i in 1:length(useGroups_list)) {
    #print names of group that is running
    print(useGroups_list[i])
    print(bgdGroups_list[i])

    #marker_pairwise for each iteration with cell names
    nam <- gsub(" ", "", paste0("marker_pairwise", "_", useGroups_list[i]))
    nam <- getMarkerFeatures(
        ArchRProj = projMulti4,
        useMatrix = "PeakMatrix",
        groupBy = "CellType_Condition",
        testMethod = "wilcoxon",
        bias = c("TSSEnrichment",
                 "log10(nFrags)",
                 "log10(Gex_nUMI)"),
        useGroups = useGroups_list[i],
        bgdGroups = bgdGroups_list[i]
    )
    markerTest[[i]] <- nam

    #Volcano for each clusters and save them
    pv <- plotMarkers(seMarker = nam,
                      name = useGroups_list[i],
                      cutOff = "FDR <= 0.1 & abs(Log2FC) >= 1",
                      plotAs = "Volcano")
    #print(pv) #print the plot
    plot_type <- "volcano-by-cluster"
    plotPDF(pv, name = gsub(" ", "-",paste(plot_type, useGroups_list[i])), width = 5, height = 5, ArchRProj = projMulti4, addDOC = FALSE)


    #marker_pairwise_significant for each iteration with cell names
    nam_significant <- gsub(" ", "", paste0("marker_pairwise_significant", "_", useGroups_list[i]))
    nam_significant <- as.data.frame(nam@assays@data@listData)
    nam_significant <- setnames(nam_significant, c("Log2FC", "Mean", "FDR", "Pval", "MeanDiff", "AUC", "MeanBGD"))
    nam_significant <- nam_significant[((nam_significant$FDR <= 0.1) & (nam_significant$Log2FC >= 1)),]
    nam_significant <- nam_significant[((nam_significant$Log2FC >= 1)),]
    nam_significant <- nam_significant %>% drop_na()
    marker_significant[[i]] <- nam_significant

    i <- i+1
}
rm(plot_type, nam_significant, nam, i, pv, useGroups_list)


#indentify marker peaks by calling addMarkerFeatures()
markersPeaks <- getMarkerFeatures(
    ArchRProj = projMulti4,
    useMatrix = "PeakMatrix",
    groupBy = "CellType_Condition", #in CellColData
    bias = c("TSSEnrichment",
             "log10(nFrags)",
             "log10(Gex_nUMI)"),
    testMethod = "wilcoxon"
)


heatmapPeaks <- plotMarkerHeatmap(
    seMarker = markersPeaks,
    cutOff = "FDR <= 0.1 & Log2FC >= 0.5",
    transpose = TRUE
)


#draw(heatmapPeaks, heatmap_legend_side = "bot", annotation_legend_side = "bot")


plotPDF(heatmapPeaks, name = "Peak-Marker-Heatmap", width = 8, height = 6, ArchRProj = projMulti4, addDOC = FALSE)
rm(heatmapPeaks)


markersPeaksCondition <- getMarkerFeatures(
    ArchRProj = projMulti4,
    useMatrix = "PeakMatrix",
    groupBy = "Conditions", #in CellColData
    bias = c("TSSEnrichment",
             "log10(nFrags)",
             "log10(Gex_nUMI)"),
    testMethod = "wilcoxon",
    useGroups = "CF", #cell type 1
    bgdGroups = "CG" #cell type 2
)


#Volcano
pv <- plotMarkers(
    seMarker = markersPeaksCondition,
    name = "CF", #that is in the parameter "groupBy" above. Up/down regulation according to this name.
    cutOff = "(FDR <= 0.1) & (Log2FC >= 1 | Log2FC <= 1)",#"FDR <= 0.1 & abs(Log2FC) >= 0.5",
    plotAs = "Volcano")


plotPDF(pv, name = "Overall-Conditions-PeakMarkers-Volcano", width = 5, height = 5, ArchRProj = projMulti4, addDOC = FALSE)
rm(pv)


markersGenesCondition <- getMarkerFeatures(
    ArchRProj = projMulti4,
    useMatrix = "GeneExpressionMatrix", #normally: "PeakMatrix", but can also generate differences in gene expression
    groupBy = "Conditions", #in CellColData
    bias = c("TSSEnrichment",
             "log10(nFrags)",
             "log10(Gex_nUMI)"),
    testMethod = "wilcoxon",
    useGroups = "CF", #cell type 1
    bgdGroups = "CG" #cell type 2
)


#Volcano
pv <- plotMarkers(
    seMarker = markersGenesCondition,
    name = "CF", #change cell type
    cutOff = "(FDR <= 0.001) & (Log2FC >= 1 | Log2FC <= 1)",#"FDR <= 0.1 & abs(Log2FC) >= 0.5",
    plotAs = "Volcano")


plotPDF(pv, name = "Overall-Conditions-GeneExpression-Volcano", width = 5, height = 5, ArchRProj = projMulti4, addDOC = FALSE)
rm(pv)


#add information about which peaks contain motifs
projMulti4 <- addMotifAnnotations(
    ArchRProj = projMulti4,
    motifSet = "cisbp",
    name = "Motif",
    force = T)


# NOTE: if only very few marker on the graph, run "options(ggrepel.max.overlaps = Inf)" in the console before hand.
options(ggrepel.max.overlaps = Inf)


motifsUp <- list() #list of motifsUp after iteration in each markerTest
motifsDo <- list() #list of motifsDo after iteration in each markerTest
df_motifsUp <- list() #list of dataframe to make the plots
df_motifsDo <- list() #list of dataframe to make the plots
plot_type <- "Markers-Motifs-Enriched-UP"


for (i in seq_along(markerTest)) {
    ################# UP:
    ## iterate for the Up-regulated motifs
    up <- gsub(" ", "", paste0("motifsUp", "_", useGroups_list[i]))
    up <- peakAnnoEnrichment(
        seMarker = markerTest[[i]],
        ArchRProj = projMulti4,
        peakAnnotation = "Motif",
        cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
    )
    motifsUp[[i]] <- up

    ## dataframe to be plotted. Also store it in list of dataframe
    df_Up <- data.frame(TF = rownames(up), mlog10Padj = assay(up)[,1])
    df_Up <- df_Up[order(df_Up$mlog10Padj, decreasing = TRUE),]
    df_Up$rank <- seq_len(nrow(df_Up))
    df_motifsUp[[i]] <- df_Up

    ## plot
    ggUp <- ggplot(df_Up, aes(rank, mlog10Padj, color = mlog10Padj)) +
        geom_point(size = 1) +
        ggrepel::geom_label_repel(
            data = df_Up[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF),
            size = 1.5,
            nudge_x = 2,
            color = "black"
        ) + theme_ArchR() +
        ylab("-log10(P-adj) Motif Enrichment") +
        xlab("Rank Sorted TFs Enriched") +
        scale_color_gradientn(colors = paletteContinuous(set = "comet"))

    plotPDF(ggUp, name = gsub(" ", "-",paste(plot_type, useGroups_list[i])), width = 5, height = 5, ArchRProj = projMulti4, addDOC = FALSE)


    ################# DOWN:
    ## iterate for the Down-regulated motifs
    down <- gsub(" ", "", paste0("motifsDo", "_", useGroups_list[i]))
    down <- peakAnnoEnrichment(
        seMarker = markerTest[[i]],
        ArchRProj = projMulti4,
        peakAnnotation = "Motif",
        cutOff = "FDR <= 0.1 & Log2FC <= -0.5"
    )
    motifsDo[[i]] <- down

    ## dataframe to be plotted. Also store it in list of dataframe
    df_Do <- data.frame(TF = rownames(down), mlog10Padj = assay(down)[,1])
    df_Do <- df_Do[order(df_Do$mlog10Padj, decreasing = TRUE),]
    df_Do$rank <- seq_len(nrow(df_Do))
    df_motifsDo[[i]] <- df_Do

    ## plot
    ggDo <- ggplot(df_Do, aes(rank, mlog10Padj, color = mlog10Padj)) +
        geom_point(size = 1) +
        ggrepel::geom_label_repel(
            data = df_Do[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF),
            size = 1.5,
            nudge_x = 2,
            color = "black"
        ) + theme_ArchR() +
        ylab("-log10(P-adj) Motif Enrichment") +
        xlab("Rank Sorted TFs Enriched") +
        scale_color_gradientn(colors = paletteContinuous(set = "comet"))

    plot_type <- "Markers-Motifs-Enriched-DOWN"
    plotPDF(ggDo, name = gsub(" ", "-",paste(plot_type, useGroups_list[i])), width = 5, height = 5, ArchRProj = projMulti4, addDOC = FALSE)

}
rm(i, plot_type, up, down, df_Up, df_Do, ggUp, ggDo, motifsDo, motifsUp)


marker_motifs_top5 <- c()
for (i in seq_along(df_motifsUp)) {
    top5 <- df_motifsUp[[i]][["TF"]][1:5]
    top5 <- gsub("_.*", "", top5) #remove the characters after the "_" in order to make them compatible for the rest of the code,
    marker_motifs_top5 <- append(marker_motifs_top5, top5)
}

for (i in seq_along(df_motifsDo)) {
    top5 <- df_motifsDo[[i]][["TF"]][1:5]
    top5 <- gsub("_.*", "", top5) #remove the characters after the "_" in order to make them compatible for the rest of the code,
    marker_motifs_top5 <- append(marker_motifs_top5, top5)
}


marker_motifs_top5 <- unique(marker_motifs_top5)
print(marker_motifs_top5)
rm(i, top5)


marker_up_enrichr <- data.frame()
marker_down_enrichr <- data.frame()
row_names <- c()
i <- 1
for (i in seq_along(df_motifsUp)) {
    topUP30 <- df_motifsUp[[i]][["TF"]][1:30]
    topUP30 <- gsub("_.*", "", topUP30) #remove the characters after the "_" in order to make them compatible for the rest of the code,
    marker_up_enrichr <- rbind(marker_up_enrichr, topUP30)
    row_names <- append(row_names,  paste("Cluster", as.character(i)))
}

i <- 1
for (i in seq_along(df_motifsUp)) {
    topDO30 <- df_motifsDo[[i]][["TF"]][1:30]
    topDO30 <- gsub("_.*", "", topDO30) #remove the characters after the "_" in order to make them compatible for the rest of the code,
    marker_down_enrichr <- rbind(marker_down_enrichr, topDO30)
}

i <- 1
my_col_names <- c()
for (i in 1:length(topUP30)) {
    my_col_names <- append(my_col_names, paste("Motif", as.character(i)))
}

rownames(marker_up_enrichr) <- row_names
colnames(marker_up_enrichr) <- my_col_names

rownames(marker_down_enrichr) <- row_names
colnames(marker_down_enrichr) <- my_col_names

# Export dataframe as csv (for report)
write.csv(marker_up_enrichr, "/Users/raphaelmauron/masters_project/doc/archr_multi/marker_up_enrichr.csv", row.names=T)

write.csv(marker_down_enrichr, "/Users/raphaelmauron/masters_project/doc/archr_multi/marker_down_enrichr.csv", row.names=T)

rm(i, topUP30, topDO30, marker_up_enrichr, marker_down_enrichr, row_names, my_col_names)


#motifs enrichment from differential testing among conditions UP and DOWN
motifsUpClusters <- peakAnnoEnrichment(
    seMarker = markersPeaksCondition, #markerTest,
    ArchRProj = projMulti4,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
)

motifsDoClusters <- peakAnnoEnrichment(
    seMarker = markersPeaksCondition, #markerTest,
    ArchRProj = projMulti4,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.1 & Log2FC <= -0.5"
)


#To prepare this data for plotting with ggplot we can create a simplified data.frame object containing the motif names, the corrected p-values, and the significance rank
df_motifsUpClusters <- data.frame(TF = rownames(motifsUpClusters), mlog10Padj = assay(motifsUpClusters)[,1])
df_motifsUpClusters <- df_motifsUpClusters[order(df_motifsUpClusters$mlog10Padj, decreasing = TRUE),]
df_motifsUpClusters$rank <- seq_len(nrow(df_motifsUpClusters))

df_motifsDoClusters <- data.frame(TF = rownames(motifsDoClusters), mlog10Padj = assay(motifsDoClusters)[,1])
df_motifsDoClusters <- df_motifsDoClusters[order(df_motifsDoClusters$mlog10Padj, decreasing = TRUE),]
df_motifsDoClusters$rank <- seq_len(nrow(df_motifsDoClusters))


#to test the most enriched motifs in the peaks
# head(df_motifsUpClusters)
# head(df_motifsDoClusters)


#Using ggplot we can plot the rank-sorted TF motifs and color them by the significance of their enrichment. Here we use ggrepel to label each TF motif.
ggUp_CF <- ggplot(df_motifsUpClusters, aes(rank, mlog10Padj, color = mlog10Padj)) +
    geom_point(size = 1) +
    ggrepel::geom_label_repel(
        data = df_motifsUpClusters[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF),
        size = 1.5,
        nudge_x = 2,
        color = "black",

    ) + theme_ArchR() +
    ylab("-log10(P-adj) Motif Enrichment") +
    xlab("Rank Sorted TFs Enriched") +
    scale_color_gradientn(colors = paletteContinuous(set = "comet")) +
    ggtitle("Motif enrichment among every cells in CF")


#Using ggplot we can plot the rank-sorted TF motifs and color them by the significance of their enrichment. Here we use ggrepel to label each TF motif.
ggDo_CF <- ggplot(df_motifsDoClusters, aes(rank, mlog10Padj, color = mlog10Padj)) +
    geom_point(size = 1) +
    ggrepel::geom_label_repel(
        data = df_motifsDoClusters[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF),
        size = 1.5,
        nudge_x = 2,
        color = "black",

    ) + theme_ArchR() +
    ylab("-log10(P-adj) Motif Enrichment") +
    xlab("Rank Sorted TFs Enriched") +
    scale_color_gradientn(colors = paletteContinuous(set = "comet")) +
    ggtitle("Motif enrichment among every cells in CG")


plotPDF(ggUp_CF, ggDo_CF,
        name = "Motifs-Enriched-plot-Condition",
        width = 5, height = 5, ArchRProj = projMulti4, addDOC = FALSE)
rm(ggUp_CF, ggDo_CF)


enrichMotifs <- peakAnnoEnrichment(
    seMarker = markersPeaks,
    ArchRProj = projMulti4,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
)


#motif enrichments across all cell groups
heatmapEM <- plotEnrichHeatmap(enrichMotifs,
                               n = 10, #number of motifs shown per cell group
                               transpose = T)


#ComplexHeatmap::draw(heatmapEM, heatmap_legend_side = "bot", annotation_legend_side = "bot")


plotPDF(heatmapEM, name = "Motifs-Enriched-Marker-Heatmap", width = 15, height = 7, ArchRProj = projMulti4, addDOC = FALSE)
#rm(heatmapEM)


#make sure we have added motif annotations to the project
## if motif not in names(projMulti4@peakAnnotation), then add it
if("Motif" %ni% names(projMulti4@peakAnnotation)){
    projMulti4 <- addMotifAnnotations(ArchRProj = projMulti4, motifSet = "cisbp", name = "Motif")
}


#add set of background peaks
projMulti4 <- addBgdPeaks(projMulti4)


#compute per-cell deviations across all of our motif annotations
#create a deviations matrix in each of our Arrow files called "MotifMatrix"
projMulti4 <- addDeviationsMatrix(
    ArchRProj = projMulti4,
    peakAnnotation = "Motif",
    force = TRUE
)


#to access those deivations
plotVarDev <- getVarDeviations(projMulti4, name = "MotifMatrix", plot = TRUE)


plotVarDev


plotPDF(plotVarDev, name = "Variable-Motif-Deviation-Scores", width = 5, height = 5, ArchRProj = projMulti4, addDOC = FALSE)
rm(plotVarDev)


#useGroups_list list of cluster to use
#bgdGroups_list

i <- 1
chromVAR_condition <- list() #list of dataframes of chromVAR variability for each cluster
chromVAR_ranked <- list() #list of dataframes of chromVAR variability for each cluster
chromVAR_top25_up <- list() #list of dataframes of chromVAR top25 UP (in CF)
chromVAR_top25_do <- list() #list of dataframes of chromVAR top25 DOWN (in CF)
chromVAR_top_abs <- data_frame() #list of dataframes of chromVAR top5 absolute values for motifs to keep for footprinting and more.

for (i in 1:length(useGroups_list)) {
    #print names of group that is running
    print(useGroups_list[i])
    print(bgdGroups_list[i])

    diffmotif <- gsub(" ", "", paste0("diffmotif", "_", useGroups_list[i]))

    #get differential variability between condition in each cluster
    diffmotif <- getMarkerFeatures(
        ArchRProj = projMulti4,
        testMethod = "wilcoxon",
        binarize = F,
        useMatrix = "MotifMatrix",
        useSeqnames = "z",
        groupBy = "CellType_Condition",
        useGroups = useGroups_list[i],
        bgdGroups = bgdGroups_list[i]
    )
    chromVAR_condition[[i]] <- diffmotif

    #get ranked motifs by MeanDiff
    rankedmotif <- gsub(" ", "", paste0("rankedmotif", "_", useGroups_list[i]))
    rankedmotif <- as.data.frame(diffmotif@elementMetadata@listData[["name"]])
    rankedmotif$MeanDiff <- as.numeric(diffmotif@assays@data@listData[["MeanDiff"]][,1])
    colnames(rankedmotif) <- c("name", "MeanDiff")
    rankedmotif <- rankedmotif[order(-rankedmotif$MeanDiff),]
    chromVAR_ranked[[i]] <- rankedmotif

    #get a top50 motifs to run Gene Ontology
    chromVAR_top25_up[[i]] <-  head(rankedmotif,25) #first 25 rows
    chromVAR_top25_do[[i]] <- tail(rankedmotif,25) #last 25 rows

    #get the top10 to keep as motif
    absolute_ranked <- rankedmotif
    absolute_ranked_top <- head(absolute_ranked["name"],10)
    absolute_ranked_top <- cbind(absolute_ranked_top, useGroups_list[i])
    colnames(absolute_ranked_top)[2] <- "cell_origin"

    chromVAR_top_abs <- rbind(chromVAR_top_abs, absolute_ranked_top)

    i <- i+1
}
rm(i, diffmotif, rankedmotif, useGroups_list, bgdGroups_list, chromVAR_condition, absolute_ranked, absolute_ranked_top)



marker_up_CV_enrichr <- data.frame()
marker_down_CV_enrichr <- data.frame()
row_names <- c()
i <- 1
for (i in seq_along(chromVAR_top25_up)) {
    topUP30 <- chromVAR_top25_up[[i]][["name"]][1:25]
    topUP30 <- gsub("_.*", "", topUP30) #remove the characters after the "_" in order to make them compatible for the rest of the code,
    marker_up_CV_enrichr <- rbind(marker_up_CV_enrichr, topUP30)
    row_names <- append(row_names,  paste("Cluster", as.character(i)))
}

i <- 1
for (i in seq_along(chromVAR_top25_up)) {
    topDO30 <- chromVAR_top25_do[[i]][["name"]][1:25]
    topDO30 <- gsub("_.*", "", topDO30) #remove the characters after the "_" in order to make them compatible for the rest of the code,
    marker_down_CV_enrichr <- rbind(marker_down_CV_enrichr, topDO30)
}

i <- 1
my_col_names <- c()
for (i in 1:length(topUP30)) {
    my_col_names <- append(my_col_names, paste("Motif", as.character(i)))
}

rownames(marker_up_CV_enrichr) <- row_names
colnames(marker_up_CV_enrichr) <- my_col_names

rownames(marker_down_CV_enrichr) <- row_names
colnames(marker_down_CV_enrichr) <- my_col_names

# Export dataframe as csv (for report)
write.csv(marker_up_CV_enrichr, "/Users/raphaelmauron/masters_project/doc/archr_multi/marker_up_CV_enrichr.csv", row.names=T)

write.csv(marker_down_CV_enrichr, "/Users/raphaelmauron/masters_project/doc/archr_multi/marker_down_CV_enrichr.csv", row.names=T)

rm(i, topUP30, topDO30, marker_up_CV_enrichr, row_names, my_col_names, chromVAR_top25_up, chromVAR_top25_do)



marker_4_table <- marker_down_CV_enrichr[c("Cluster 1", "Cluster 2", "Cluster 8", "Cluster 9"), ]

new_df <- as.data.frame(t(marker_4_table))
new_df <- head(new_df,20)
new_df
write.csv(new_df, "/Users/raphaelmauron/masters_project/doc/archr_multi/marker_down_CV_enrichr_4_clusters.csv", row.names=T)

marker_4_clusters <- c()
for (i in 1:4) {
    marker_4_clusters <- append(marker_4_clusters, new_df[,i])
}
marker_4_clusters <- sort(unique(marker_4_clusters))

rm(i, new_df)


chromVAR_top_abs$name <- gsub("_.*", "", chromVAR_top_abs$name)
chromVAR_top_abs$cell_origin <- gsub("_x_CF", "", chromVAR_top_abs$cell_origin)
split_df <- split(chromVAR_top_abs, chromVAR_top_abs$cell_origin)

new_cols <- lapply(split_df, function(df) {
    names_vec <- as.character(df$name)
    new_col <- setNames(names_vec, df$cell_origin[1])
    return(new_col)
})


# combine the new columns into a single dataframe
new_df <- as.data.frame(new_cols)


write.csv(new_df, "/Users/raphaelmauron/masters_project/doc/archr_multi/abs_top10_CV_markers.csv", row.names=T)


rm(split_df, new_cols, new_df, n)


chromVAR_top_abs_unique <- unique(chromVAR_top_abs$name)
chromVAR_top_abs_unique <- sub("_.*", "", chromVAR_top_abs_unique)


#to analyze subset of features
# motifs <- c() #instead of using "my_markers" you can specify another list here
### maybe take the marker_motifs_top5
motifs <- marker_4_clusters
markerMotifs <- getFeatures(projMulti4, select = paste(motifs, collapse="|"), useMatrix = "MotifMatrix")
rm(chromVAR_top_abs_unique)


#want to remove a specific marker
markerMotifs <- grep("z:", markerMotifs, value = TRUE) #keeps only z-values
#markerMotifs <- markerMotifs[markerMotifs %ni% "z:SREBF1_22"] #remove this specific marker


#the distribution of chromVAR deviation scores for each cluster
#Notice that we supply the impute weights that we calculated previously during our gene score analyses
#these impute weights allow us to smooth the signal across nearby cells which is helpful in the context of our sparse scATAC-seq data

p <- plotGroups(ArchRProj = projMulti4,
                groupBy = "Conditions" ,#"predictedGroup_Un_3", #"Conditions",
                colorBy = "MotifMatrix",
                name = markerMotifs,
                imputeWeights = getImputeWeights(projMulti4)
)


#cowplot to plot the distributions of all motifs in a single plot
p2 <- lapply(seq_along(p), function(x){
    if(x != 1){
        p[[x]] + guides(color = "none", fill = "none") +
            theme_ArchR(baseSize = 6) +
            theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm")) +
            theme(
                axis.text.y=element_blank(),
                axis.ticks.y=element_blank(),
                axis.title.y=element_blank()
            ) + ylab("")
    }else{
        p[[x]] + guides(color = "none", fill = "none") +
            theme_ArchR(baseSize = 6) +
            theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm")) +
            theme(
                axis.ticks.y=element_blank(),
                axis.title.y=element_blank()
            ) + ylab("")
    }
})
#do.call(cowplot::plot_grid, c(list(nrow = 1, rel_widths = c(2, rep(1, length(p2) - 1))),p2))


plotPDF(p, name = "Plot-Groups-Conditions-Deviations-w-Imputation", width = 5, height = 5, ArchRProj = projMulti4, addDOC = FALSE)
rm(p, p2)


#overlay the z-scores on our UMAP embedding as we've done previously for gene scores
p <- plotEmbedding(
    ArchRProj = projMulti4,
    colorBy = "MotifMatrix",
    name = sort(markerMotifs),
    embedding = "UMAPHarmony_ATAC", #UMAPHarmony", #can change embedding here
    imputeWeights = getImputeWeights(projMulti4)
)


#can plot all of these motif UMAPs using cowplot
p2 <- lapply(p, function(x){
    x + guides(color = FALSE, fill = FALSE) +
        theme_ArchR(baseSize = 6.5) +
        theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
        theme(
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank()
        )
})
#do.call(cowplot::plot_grid, c(list(ncol = 3),p2))


plotPDF(p, name = "Plot-Groups-Deviations-w-Z-score-embedding", width = 5, height = 5, ArchRProj = projMulti4, addDOC = FALSE)
rm(p, p2)


markerRNA <- getFeatures(projMulti4,
                         select = paste(motifs, collapse="|"),
                         useMatrix = "GeneExpressionMatrix")

markerRNA <- markerRNA[markerRNA %in% motifs] #exclude some motifs that are not in "motifs"
sort(markerRNA)
#sort(motifs) #to compare if needed


p <- plotEmbedding(
    ArchRProj = projMulti4,
    colorBy = "GeneExpressionMatrix",
    name = sort(markerRNA),
    embedding = "UMAPHarmony_ATAC",
    imputeWeights = getImputeWeights(projMulti4)
)


#cowplot to plot single features
p2 <- lapply(p, function(x){
    x + guides(color = FALSE, fill = FALSE) +
        theme_ArchR(baseSize = 6.5) +
        theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
        theme(
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank()
        )
})
#do.call(cowplot::plot_grid, c(list(ncol = 3),p2))


plotPDF(p, name = "Plot-Groups-Deviations-w-GeneExpression-embedding", width = 5, height = 5, ArchRProj = projMulti4, addDOC = FALSE)
rm(p, p2)


#getAvailableMatrices(projMulti4)


markerRNA <- getFeatures(projMulti4,
                         select = paste(motifs, collapse="|"),
                         useMatrix = "GeneIntegrationMatrix_db_3")
markerRNA <- markerRNA[markerRNA %in% motifs] #exclude some motifs that are not in "motifs"


p <- plotEmbedding(
    ArchRProj = projMulti4,
    colorBy = "GeneIntegrationMatrix_db_3",
    name = sort(markerRNA),
    embedding = "UMAPHarmony_ATAC",
    continuousSet = "blueYellow",
    imputeWeights = getImputeWeights(projMulti4)
)


p2 <- lapply(p, function(x){
    x + guides(color = FALSE, fill = FALSE) +
        theme_ArchR(baseSize = 6.5) +
        theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
        theme(
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank()
        )
})
#do.call(cowplot::plot_grid, c(list(ncol = 3),p2))


plotPDF(plotList = p,
        name = "Plot-Groups-Deviations-w-linked-gene-expression-of-TF-embedding.pdf",
        ArchRProj = projMulti4,
        addDOC = FALSE, width = 5, height = 5)
rm(p, p2)


#create a GRangesList
motifPositions <- getPositions(projMulti4)


#list of marker to use for the subset
# motifs <- c() #specify another motifs list here if want others than my_markers

markerMotifs <- unlist(lapply(motifs, function(x) grep(x, names(motifPositions), value = TRUE)))
markerMotifs <- markerMotifs[markerMotifs %ni% c("Bhlha15_87", "Ctcf_146", "Ctcfl_820", "Foxd4_833", "Foxl1_354", "Foxl1_355", "Foxl1_356", "Foxl1_357",        "Foxl1_358", "Foxl1_359", "Foxl1_360", "Foxl1_361", "Foxl1_362", "Foxl1_363", "Foxl1_364", "Foxl1_365", "Foxl1_366", "Foxl1_367", "Foxl1_368", "Foxl1_369", "Foxl1_370", "Foxl1_371", "Foxl1_372", "Foxl1_373", "Foxl1_374", "Foxl1_375", "Foxl1_376", "Foxl1_377", "Foxl1_378", "Foxl1_379", "Foxl1_380", "Foxl1_381", "Gm5294_331", "Gm5294_332", "Gm5294_333", "Gm5294_334", "Gm5294_335", "Gm5294_336", "Gm5294_337", "Gm5294_338", "Gm5294_339", "Gm5294_340", "Gm5294_341", "Gm5294_342", "Gm5294_343", "Gm5294_344", "Gm5294_345", "Gm5294_346", "Gm5294_347", "Gm5294_348", "Gm5294_349", "Gm5294_350", "Gm5294_351", "Gm5294_352", "Gm5294_353", "Gm5294", "Jund_135", "Mesp1_58", "Mesp2_57", "Nfia_862", "Nfic_724", "Tcf4_88", "X9430076C15Rik_128", "Zeb1_139", "Zfp105_218")] #those motifs were removed because no connection with neuronal function after UniProt check
#markerMotifs <- markerMotifs[markerMotifs %in% chromVAR_ranked[[1]][["name"]]]


#compute footprints for the subset of marker motifs that we previously selected using the getFootprints() function
#recommended to perform footprinting on a subset of motifs rather than all motifs
#provide the subset of motifs to footprint via the positions parameter
seFoot <- getFootprints(
    ArchRProj = projMulti4,
    positions = motifPositions[markerMotifs], #specify the subset of motifs
    groupBy = "CellType_Condition",
    useGroups = c("Astroependymal_cells_x_CF",
                  "Astroependymal_cells_x_CG",
                  "Cholinergic,_monoaminergic_and_peptidergic_neurons_x_CF",
                  "Cholinergic,_monoaminergic_and_peptidergic_neurons_x_CG",
                  "Telencephalon_interneurons_x_CF",
                  "Telencephalon_interneurons_x_CG",
                  "Telencephalon_projecting_neurons_x_CF",
                  "Telencephalon_projecting_neurons_x_CG"
    ) #groups from "groupBy" to keep for the footprinting
)


#footprints
plotFootprints(
    seFoot = seFoot,
    ArchRProj = projMulti4,
    normMethod = "Subtract",
    plotName = "Footprints-Subtract-Bias",
    addDOC = FALSE,
    smoothWindow = 5
)


#other strategy to normalize: divide footprint signal by the Tn5 bias signal
plotFootprints(
    seFoot = seFoot,
    ArchRProj = projMulti4,
    normMethod = "Divide",
    plotName = "Footprints-Divide-Bias",
    addDOC = FALSE,
    smoothWindow = 5
)


#no normalization
plotFootprints(
    seFoot = seFoot,
    ArchRProj = projMulti4,
    normMethod = "None", #says that no normalization is performed
    plotName = "Footprints-No-Normalization",
    addDOC = FALSE,
    smoothWindow = 5
)


#create TSS insertion profiles without normalization for Tn5 bias
#main difference from our previous analyses is that we specify flank = 2000 to extend these footprints 2000 bp on either side of each TSS
seTSS <- getFootprints(
    ArchRProj = projMulti4,
    positions = GRangesList(TSS = getTSS(projMulti4)), #specify
    groupBy = "CellType_Condition",
    flank = 2000 #extend the footprints on both side of each TSS
)


plotFootprints(
    seFoot = seTSS,
    ArchRProj = projMulti4,
    normMethod = "None",
    plotName = "TSS-No-Normalization",
    addDOC = FALSE,
    flank = 2000,
    flankNorm = 100
)


#add co-accessibility peaks to the object
projMulti4 <- addCoAccessibility(
    ArchRProj = projMulti4,
    reducedDims = "Harmony_ATAC" #can change dimensionality reduction here
)


#retrieve the co-accessibility information
cA <- getCoAccessibility(
    ArchRProj = projMulti4,
    corCutOff = 0.5,
    resolution = 1,
    returnLoops = FALSE #returns dataframe only when = FALSE
)


#cA


#metadata(cA)[[1]]


#resolution = 1 creates loops that connect the center of each peak
cA <- getCoAccessibility(
    ArchRProj = projMulti4,
    corCutOff = 0.5,
    resolution = 1, #when decrease the resolution (= 1000 instead of = 1), can help with overplotting of co-accessibility interactions
    returnLoops = TRUE # return co-accessibility in the form of a loop track
)


#cA[[1]]


# plot every name from the list "motifs"
#motifs_for_co_acc <- motifs[56:68]
#my_markers the ones from RR3
p <- plotBrowserTrack(
    ArchRProj = projMulti4,
    geneSymbol = c(my_markers, motifs), #specify motifs to plot here
    upstream = 50000,
    downstream = 50000,
    loops = getCoAccessibility(projMulti4),
    groupBy = "CellType_Condition",
    useGroups = c("Astroependymal_cells_x_CF",
                  "Astroependymal_cells_x_CG",
                  "Cholinergic,_monoaminergic_and_peptidergic_neurons_x_CF",
                  "Cholinergic,_monoaminergic_and_peptidergic_neurons_x_CG",
                  "Telencephalon_interneurons_x_CF",
                  "Telencephalon_interneurons_x_CG",
                  "Telencephalon_projecting_neurons_x_CF",
                  "Telencephalon_projecting_neurons_x_CG"
    ) #groups from "groupBy" to keep for the footprinting
)


#to open in browser
# grid::grid.newpage()
# grid::grid.draw(p$Dkk3)


plotPDF(plotList = p,
        name = "Plot-Tracks-Marker-Genes-with-CoAccessibility.pdf",
        ArchRProj = projMulti4,
        addDOC = FALSE, width = 7, height = 5)
rm(p)


p <- plotBrowserTrack(
    ArchRProj = projMulti4,
    groupBy = "CellType_Condition",
    useGroups = c("Astroependymal_cells_x_CF",
                  "Astroependymal_cells_x_CG",
                  "Cholinergic,_monoaminergic_and_peptidergic_neurons_x_CF",
                  "Cholinergic,_monoaminergic_and_peptidergic_neurons_x_CG",
                  "Telencephalon_interneurons_x_CF",
                  "Telencephalon_interneurons_x_CG",
                  "Telencephalon_projecting_neurons_x_CF",
                  "Telencephalon_projecting_neurons_x_CG"
    ), #groups from "groupBy" to keep for the footprinting
    geneSymbol = c(my_markers, motifs) ,#markerGenes,
    upstream = 50000,
    downstream = 50000,
    loops = getPeak2GeneLinks(projMulti4) #the links were stored in projMulti4
)


 browser track
# grid::grid.newpage()
# grid::grid.draw(p$Dkk3)


plotPDF(plotList = p,
        name = "Plot-Tracks-Marker-Genes-with-Peak2GeneLinks.pdf",
        ArchRProj = projMulti4,
        addDOC = FALSE, width = 5, height = 5)
rm(p)


p <- plotPeak2GeneHeatmap(ArchRProj = projMulti4,
                          groupBy = "CellType_Condition"
                          # useGroups = c("Astroependymal_cells_x_CF",
                          #               "Astroependymal_cells_x_CG",
                          #               "Cholinergic,_monoaminergic_and_peptidergic_neurons_x_CF",
                          #               "Cholinergic,_monoaminergic_and_peptidergic_neurons_x_CG",
                          #               "Telencephalon_interneurons_x_CF",
                          #               "Telencephalon_interneurons_x_CG",
                          #               "Telencephalon_projecting_neurons_x_CF",
                          #               "Telencephalon_projecting_neurons_x_CG"
                          #               ) #groups from "groupBy" to keep for the footprinting)
)


plotPDF(p,
        name = "Heatmap-coaccessibility-ATAC-RNA.pdf",
        ArchRProj = projMulti4,
        addDOC = FALSE,
        height = 10,
        width = 15)
rm(p)


#retrieve this data (deviant TF motifs)
seGroupMotif <- getGroupSE(ArchRProj = projMulti4,
                           useMatrix = "MotifMatrix",
                           groupBy = "CellType_Condition")


#since it comes from MotifMatrix, has "deviations" & "z" corresponding to row and z-score deviations
#seGroupMotif
#subset this SummarizedExperiment to just the "deviation z-scores"
seZ <- seGroupMotif[rowData(seGroupMotif)$seqnames=="z",]


#identify the maximum delta in z-score between all clusters
#will be helpful in stratifying motifs based on the degree of variation observed across clusters
rowData(seZ)$maxDelta <- lapply(seq_len(ncol(seZ)), function(x){
    rowMaxs(assay(seZ) - assay(seZ)[,x])
}) %>% Reduce("cbind", .) %>% rowMaxs


corGSM_MM <- correlateMatrices(
    ArchRProj = projMulti4,
    useMatrix1 = "GeneExpressionMatrix",
    useMatrix2 = "MotifMatrix",
    reducedDims = "Harmony"
)


#look at the dataframe
#corGSM_MM
#can perform the same analysis using the GeneIntegrationMatrix instead of the GeneScoreMatrix
corGIM_MM <- correlateMatrices(
    ArchRProj = projMulti4,
    useMatrix1 = "GeneIntegrationMatrix_db_3", #do we have that in our case?
    useMatrix2 = "MotifMatrix",
    reducedDims = "Harmony"
)


#corGIM_MM
corGSM_MM$maxDelta <- rowData(seZ)[match(corGSM_MM$MotifMatrix_name, rowData(seZ)$name), "maxDelta"]
corGIM_MM$maxDelta <- rowData(seZ)[match(corGIM_MM$MotifMatrix_name, rowData(seZ)$name), "maxDelta"]


corGSM_MM <- corGSM_MM[order(abs(corGSM_MM$cor), decreasing = TRUE), ]
corGSM_MM <- corGSM_MM[which(!duplicated(gsub("\\-.*","",corGSM_MM[,"MotifMatrix_name"]))), ]
corGSM_MM$TFRegulator <- "NO"
corGSM_MM$TFRegulator[which(corGSM_MM$cor > 0.5 & corGSM_MM$padj < 0.01 & corGSM_MM$maxDelta > quantile(corGSM_MM$maxDelta, 0.75))] <- "YES"
sort(corGSM_MM[corGSM_MM$TFRegulator=="YES",1])


 the positive TF
p <- ggplot(data.frame(corGSM_MM), aes(cor, maxDelta, color = TFRegulator)) +
    geom_point() +
    theme_ArchR() +
    geom_vline(xintercept = 0, lty = "dashed") +
    scale_color_manual(values = c("NO"="darkgrey", "YES"="firebrick3")) +
    xlab("Correlation To Gene Score") +
    ylab("Max TF Motif Delta") +
    scale_y_continuous(
        expand = c(0,0),
        limits = c(0, max(corGSM_MM$maxDelta)*1.05)
    )


plotPDF(p,
        name = "positive-TF-GSM.pdf",
        ArchRProj = projMulti4,
        addDOC = FALSE,
        height = 5,
        width = 5)
rm(p)


#same from the GeneIntegrationMatrix
corGIM_MM <- corGIM_MM[order(abs(corGIM_MM$cor), decreasing = TRUE), ]
corGIM_MM <- corGIM_MM[which(!duplicated(gsub("\\-.*","",corGIM_MM[,"MotifMatrix_name"]))), ]
corGIM_MM$TFRegulator <- "NO"
corGIM_MM$TFRegulator[which(corGIM_MM$cor > 0.5 & corGIM_MM$padj < 0.01 & corGIM_MM$maxDelta > quantile(corGIM_MM$maxDelta, 0.75))] <- "YES"
sort(corGIM_MM[corGIM_MM$TFRegulator=="YES",1])


p <- ggplot(data.frame(corGIM_MM), aes(cor, maxDelta, color = TFRegulator)) +
    geom_point() +
    theme_ArchR() +
    geom_vline(xintercept = 0, lty = "dashed") +
    scale_color_manual(values = c("NO"="darkgrey", "YES"="firebrick3")) +
    xlab("Correlation To Gene Expression") +
    ylab("Max TF Motif Delta") +
    scale_y_continuous(
        expand = c(0,0),
        limits = c(0, max(corGIM_MM$maxDelta)*1.05)
    )


plotPDF(p,
        name = "positive-TF-GIM.pdf",
        ArchRProj = projMulti4,
        addDOC = FALSE,
        height = 5,
        width = 5)
rm(p)


saveArchRProject(ArchRProj = projMulti4, outputDirectory = "Save-ProjMulti4", load = F)


sessionInfo()
