# _CellChat_Execution_EC.R
# source activate cellchat
# R

suppressMessages(library(CellChat))
suppressMessages(library(patchwork))
suppressMessages(library(reticulate))
options(stringsAsFactors = FALSE)

# Setting the ligand-receptor interaction database
CellChatDB <- CellChatDB.mouse

# Importing Scanpy-processed Data
anndata <- import("anndata", convert=FALSE)
anndata_object <- anndata$read_h5ad("/data/Projects/phenomata/01.Projects/11.Vascular_Aging/03.Scanpy/CellChat/test3_forCellChat.h5ad") # scran normalized data matrix input

data.input <- t(py_to_r(anndata_object$X))
rownames(data.input) <- rownames(py_to_r(anndata_object$var))
colnames(data.input) <- rownames(py_to_r(anndata_object$obs))

# access meta data
meta.data <- py_to_r(anndata_object$obs)
meta <- meta.data

# Create a CellChat Object using Data matrix as input
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "celltype")

## levels(cellchat@idents)
## [1] "EC1" "EC4" "EC2" "EC3"

# Preprocessing the expression data for CellChat Analysis

cellchat@DB <- CellChatDB # https://github.com/sqjin/CellChat/issues/96#issuecomment-885818183
cellchat <- subsetData(cellchat)

future::plan("multiprocess", workers = 20)
#Warning message:
#Strategy 'multiprocess' is deprecated in future (>= 1.20.0). Instead, explicitly specify either 'multisession' or 'multicore'. In the current R session, 'multiprocess' equals 'multicore'.

cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)


https://github.com/sqjin/CellChat/issues/149 # ==> EC1, EC2, EC3, EC4 이렇게 ordering될 수 있도록 하는 방법


# (i) Projecting gene expression data on to Protein-protein interaction (PPI) network 
# (https://htmlpreview.github.io/?https://github.com/sqjin/CellChat/blob/master/tutorial/CellChat-vignette.html)

## project gene expression data onto PPI (Optional: when running it, USER should set `raw.use = FALSE` in the function `computeCommunProb()` in order to use the projected data)
## cellchat <- projectData(cellchat, PPI.mouse)
## TO be continued

# (ii) "NOT" Projecting gene expression data on to Protein-protein interaction (PPI) network 

# Inferring of cell-cell communication network

## The number of inferred ligand-receptor pairs clearly depends on the method for calculating the average gene expression per cell group. By default, CellChat uses a statistically robust mean method called ‘trimean’, which produces fewer interactions than other methods. However, we find that CellChat performs well at predicting stronger interactions, which is very helpful for narrowing down on interactions for further experimental validations. In computeCommunProb, we provide an option for using other methods, such as 5% and 10% truncated mean, to calculating the average gene expression. Of note, ‘trimean’ approximates 25% truncated mean, implying that the average gene expression is zero if the percent of expressed cells in one group is less than 25%. To use 10% truncated mean, USER can set type = "truncatedMean" and trim = 0.1. The function computeAveExpr can help to check the average expression of signaling genes of interest, e.g, computeAveExpr(cellchat, features = c("CXCL12","CXCR4"), type =  "truncatedMean", trim = 0.1).

## When analyzing unsorted single-cell transcriptomes, under the assumption that abundant cell populations tend to send collectively stronger signals than the rare cell populations, CellChat can also consider the effect of cell proportion in each cell group in the probability calculation. USER can set population.size = TRUE.

cellchat <- computeCommunProb(cellchat) # This may take about 10 minutes with "future::plan("multiprocess", workers = 10)"
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10) # 이거 일단 보류

cellchat <- computeCommunProbPathway(cellchat)

cellchat <- aggregateNet(cellchat)

groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")


mat <- cellchat@net$weight
par(mfrow = c(2,2), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

#### 

cellchat@netP$pathways # ==> All the signaling pathways showing significant communications can be accessed using this

pathways.show <- c("COLLAGEN") 
pathways.show <- c("CDH5")
pathways.show <- c("THBS")

# Hierarchy plot
# Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
vertex.receiver = seq(1,4) # a numeric vector. 
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver)

# Circle plot
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")

# Chord Diagram
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")

# Heatmap
par(mfrow=c(1,1))
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")

netAnalysis_contribution(cellchat, signaling = pathways.show)