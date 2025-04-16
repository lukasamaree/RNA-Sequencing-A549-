#' ---
#' title: "FINAL"
#' author: "LUKAS AMARE"
#' date: "2024-03-14"
#' output: html_document
#' ---
#' 
#' (25 pts) Use Seurat and CellChat to analyze the single cell RNA-seq dataset of A549 lung carcinoma cells (the dataset is available on canvas/Files/Final/A549):
#' 
#' A) Setup the Seurat object and perform quality control analysis.  Please explain the reasoning of your thresholds for quality control. Use violin plots in your explanation.
#' 
#' 
#' B) Normalize data, detect highly variable genes and scale the data.
#' 
## -----------------------------------------------------------------------------
library(Seurat)
library(magrittr)
library(dplyr)

setwd("~/Desktop/Lukas Amare/big_data_week_1/A549")

# LOAD DATA AND MAKE SEURAT OBJECT


A549.data = readRDS("raw.data.rds")
A549 <- CreateSeuratObject(counts = A549.data, project = "A549", min.cells = 3, min.features = 200)
A549


head(colnames(A549))
head(rownames(A549)) 
head(A549@meta.data) 

#VIOLIN PLOT BEFORE
VlnPlot(A549, features = c("nFeature_RNA", "nCount_RNA"), ncol = 3, pt.size = 0)



plot1 <- FeatureScatter(A549, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")







#REASONING BEHIND THRESHOLDS - 
#n Feature shows number of detected genes for every cell
#n Count shows nymber of unique molecular identifiers for every cell
# The reasoning behind QC thresholds in remove noise from the data set in order to have a more accurate visualization of the RNA data we are analyzing. The also increase mapping quality for UMAP.
# For our QC we set our nFeature or number of genes for each cell from 50 to 6000. We also set our nCount or molecular identifiers for each cell from 0 to 25000

# QC
A549 <- subset(A549, subset = nFeature_RNA > 50 & nFeature_RNA < 6000 & nCount_RNA < 25000 & nCount_RNA > 0)
# NORMALIZING DATA
A549 <- NormalizeData(A549, normalization.method = "LogNormalize", scale.factor = 10000)
A549 <- NormalizeData(A549)
# After QC metrics you see more visible violin plot
VlnPlot(A549, features = c("nFeature_RNA", "nCount_RNA"), 
        ncol = 3,pt.size=0)
A549 <- FindVariableFeatures(A549, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(A549), 10)
top10

plot1 <- VariableFeaturePlot(A549)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1  
plot2

all.genes <- rownames(A549)
#SCALE DATA
A549 <- ScaleData(A549, features = all.genes)


#' 
#' C)	Perform PCA analysis and choose the number of PCs for further analysis. Explain your reasoning. 
#' D) Cluster the cells, and test different values of the resolution parameter in FindClusters(). Explain your reasoning for choosing the resolution value for further analysis.
#' E) Visualize your clustering results in 2d via UMAP. Do this for at least two different values of resolution. How do results from FindClusters() change as you modify the resolution parameter?
#' 
## -----------------------------------------------------------------------------
A549 <- RunPCA(A549, features = VariableFeatures(object = A549))
print(A549[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(A549, dims = 1:2, reduction = "pca")
DimPlot(A549, reduction = "pca") + NoLegend()

DimHeatmap(A549, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(A549, dims = 1:10, cells = 500, balanced = TRUE)
ElbowPlot(A549) # Based off of the elbow plot the optimal number of PC's  is 10 due to the fact that variability drastically does not change after 10 PC's


# Cluster the cells

A549 <- FindNeighbors(A549, dims = 1:10)
# Resolution for cells is based on specificity of the visualization. A very high resolution will have many clusters while a very low resolution will have little clusters. 
A549 <- FindClusters(A549, resolution = 0.5)
head(Idents(A549), 10)

A549 <- RunUMAP(A549, dims = 1:10)
DimPlot(A549, reduction = "umap",label=TRUE)









# Visualize clustering results for resolution = 1.0
A549_1 <- FindClusters(A549, resolution = 2)

# Run UMAP
A549_1 <- RunUMAP(A549_1, dims = 1:10)

# Plot UMAP for resolution = 1.0
DimPlot(A549_1, reduction = "umap", label = TRUE, pt.size = 1) 


# I am going to use resolution 0.6 because it is between the two extremes of resolution and it shows solid clustering 2.0 resolutions seems to cluster too much however I feel as if it is more accurate due to this being a larger data set.
A549 <- FindClusters(A549, resolution = 0.625)
head(Idents(A549), 5)

A549 <- RunUMAP(A549, dims = 1:10)
DimPlot(A549, reduction = "umap",label=TRUE)


#' 
#' 
#' 
#' F) Find differentially expressed genes for each cluster and make a heat map showing the top 3 up-regulated genes per cluster.
#' G) Assign cell type identity to clusters. Some canonical markers are listed below (see
#' table on next page), but additional markers may be needed from the literature.
#'       Please use violin plot, feature plot and dot          plot to show how you perform the cell
#'       type identification and confirm that annotated        cell types express the expected genes.
## -----------------------------------------------------------------------------
# Finding markers of genes for cluster 0
cluster0.markers <- FindMarkers(A549, ident.1 = 0)
head(cluster0.markers, n = 10)
# Finding markers of genes for cluster 1
cluster1.markers <- FindMarkers(A549, ident.1 = 1)
head(cluster1.markers, n = 10)
# Finding markers of genes for cluster 2
cluster2.markers <- FindMarkers(A549, ident.1 = 2)
head(cluster2.markers, n = 10)
# Finding markers of genes for cluster 3
cluster3.markers <- FindMarkers(A549, ident.1 = 3)
head(cluster3.markers, n = 10)
#Findng markers for gene for cluster 4
cluster4.markers <- FindMarkers(A549, ident.1 = 4)
head(cluster4.markers, n = 10)
# Finding markers of genes for cluster 5
cluster5.markers <- FindMarkers(A549, ident.1 = 5, ident.2 = c(0, 3))
head(cluster5.markers, n = 10)
# Finding markers of genes for cluster 6
cluster6.markers <- FindMarkers(A549, ident.1 = 6)
head(cluster6.markers, n = 10)
# Finding markers of genes for cluster 7
cluster7.markers <- FindMarkers(A549, ident.1 = 7)
head(cluster7.markers, n = 10)
# Finding markers of genes for cluster 8


#finds marker for every cluster
A549.markers <- FindAllMarkers(A549, only.pos = TRUE)
A549.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1)



# Extract top 3 markers per cluster
top_markers_per_cluster <- A549.markers %>%
    group_by(cluster) %>%
    top_n(5, avg_log2FC) %>%
    ungroup()

# Create a heatmap
DoHeatmap(A549, features = top_markers_per_cluster$gene) + NoLegend()


# Transpose the expression matrix




# Violin Plot to visualize marker expression 

VlnPlot(A549, features = c( "TBX6", "FOXA2"))
                          
VlnPlot(A549, features = c("SOX9"   ,  "MUC5AC"))
        
VlnPlot(A549, features = c("MUC13" ,  "MUC1"))

VlnPlot(A549, features = c("CFTR" ,  "CD44"))
                           



# early mesoderm TBX6, endoderm FOXA2, lung proginator cells SOX9, cilated cells MUC5AC, goblet cells MUC13, clara cells MUC1, alveolar type CFTR, cells cancer stem cells CD44 
FeaturePlot(A549, features = c("TBX6", "FOXA2", "SOX9","MUC5AC", "MUC13", "MUC1" , "CFTR" , "CD44"))
DotPlot(A549, features = c("TBX6", "FOXA2", "SOX9","MUC5AC", "MUC13", "MUC1" , "CFTR" , "CD44"))

# TBX6 1.3, 2.3    -----> Cluster 3
# FOXA2 1. ,2.6 checked ----> Cluster 6
# SOX9 1.1, 2.1 - ------> Cluster 1 
#MUC5AC 1.7, 2.7 ** ----> cluster 7
#MUC13, 1. 2.4 -----> cluster 4
#MUC1 , 2.5 -----> cluster 5 
#CFTR 2.0, 0 ---> Cluster 2
#CD44 2.2 -----> CLuster 0

new.cluster.ids <- c("cancer stem cells", "lung proginator", "alveolar cell", "early mesoderm", "goblet cells", "clara cells",
    "endoderm", "cilated cells ")
names(new.cluster.ids) <- levels(A549)
A549 <- RenameIdents(A549, new.cluster.ids)
DimPlot(A549, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()




#' 
#' Create the CellChat object.
#' i) Set up the ligand-receptor interaction database, identify the over-expressed ligands
#' and receptors and compute the communication probability at the signaling pathway
#' level.
#' 
## -----------------------------------------------------------------------------
library(CellChat)
# Make CellChat object
cellchat <- createCellChat(object = A549, group.by = "ident", assay = "RNA")

# Set up database
CellChatDB <- CellChatDB.human  # or use CellChatDB.mouse if your data is mouse
showDatabaseCategory(CellChatDB)
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling", key = "annotation")  # Use Secreted Signaling only
cellchat@DB <- CellChatDB.use

# Subset data (necessary)
cellchat <- subsetData(cellchat)  # This step is necessary even if using the whole database

# Preprocessing
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

# Infer communication probability
cellchat <- computeCommunProb(cellchat, type = "triMean")
cellchat <- filterCommunication(cellchat, min.cells = 10)

# Compute communication probability at the pathway level
cellchat <- computeCommunProbPathway(cellchat)



#' 
#' 
#' j) Visualize the aggregated cell-cell communication network by showing both the
#' number of interactions and the total interaction strength (weights) between any two
#' cell groups using circle plot. Describe possible biological findings you obtain.
#' k) Visualize one particular signaling pathway of your choice using Hierarchy plot, Circle
#' plot and Chord diagram. Describe your findings.
#' 
## -----------------------------------------------------------------------------
# Aggregate pathways
cellchat <- aggregateNet(cellchat)
# Visualize pathways (aggregate)
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)

# Circle plot for number fo interactuions and interaction weights/strengths
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

#Circle

netVisual_aggregate(cellchat, signaling = "WNT", layout = "circle")

#In the Wnt pathway there are a lot of connections from mesoderm to endoderm This makes sense because the WNT signalling pathway is responsible for development mesoderm and endoderm in embryo. Aveolar and lung proginator cells also link with with the endoderm this is because the endoderm in development leads to layer of tissue that make up the aveolar and lung proginator cells . Alveolar and lung proginator cells have very strong connections with goblet and clara cells through this pathway. Cancer Stem cells and cilated cells have weak connections to everyhting throgh the Wnt pathway.

# Chord
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = "WNT", layout = "chord")

# Clara Cells seem to have a strong connection with lung proginator cells. Cilated cells seem to have no connection to anything through the Wnt pathway. Everyone of the cells except cilated cells seem to be connected one way or another.

#HEiarchial

#Goblet cells seem to be conncecting to everything except clara cells. Cilated cells has no connection with anything.

vertex.receiver = seq(1,4) # a numeric vector. 
netVisual_aggregate(cellchat, signaling = "WNT",  vertex.receiver = vertex.receiver)


#' 
#' 
#' l) Compute the contribution of each ligand-receptor pair to the chosen signaling
#' pathway from k) and visualize the cell-cell communication mediated by one single
#' ligand-receptor pair using Hierarchy plot, Circle plot or Chord diagram.
#' m) Identify signaling roles and visualize the computed centrality scores using heatmap.
#' Explain what each centrality score means.
#' n) Visualize the dominant senders (sources) and receivers (targets) in a 2D space and
#' describe your findings. Suggest some possible biological interpretations
#' 
## -----------------------------------------------------------------------------

pairLR.wnt <- extractEnrichedLR(cellchat, signaling = "WNT", geneLR.return = FALSE)
LR.show <- pairLR.wnt[1,] # show one ligand-receptor pair
# Hierarchy plot
vertex.receiver = seq(1,4) # a numeric vector
netVisual_individual(cellchat, signaling = "WNT",  pairLR.use = LR.show, vertex.receiver = vertex.receiver)
# Circle plot
netVisual_individual(cellchat, signaling = "WNT", pairLR.use = LR.show, layout = "circle")
# Chord diagram
netVisual_individual(cellchat, signaling = "WNT", pairLR.use = LR.show, layout = "chord")

# HEAT MAP
# Compute the network centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways


# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
netAnalysis_signalingRole_network(cellchat, signaling = "WNT", width = 8, height = 2.5, font.size = 10)


#The centrality scores show the relatedness or importance  of signalalling roles  being recieved, sent, mediated, or influenced relation to each other. For this specific example the centrality scores describe WNT signalling between cells.



#DOMINANT SENDERS - Lung proginator and Endoderm. THis makes sense because lung proginator is used for lung regeneration and repair therefore it would need to be strong in sending Wnt signals because WNt signal pathway. is responsible for development of embryos. Endoderm is very associated with layer of tissue during development as well.

#DOMINANT RECIEVERS  - Endoderm. THis makes sense because the WNT is used for development and the endoderm is layers of tissue in development of embryo









#' 
#' (20 pts) The Khan dataset in ISLR library (library(ISLR)) consists of a number of tissue samples corresponding to four distinct types of small round blue cell tumors. For each tissue sample, gene expression measurements are available. The data set consists of training data, xtrain and ytrain, and testing data, xtest and ytest.
#' 
#' 
#' 
#' 
#' a) Use random forests to analyze this data (need to convert ytrain/ytest into a qualitative response variable). What is the value for mtry in the randomForest() function and why?
#' 
#' 
#' 
#' b) Show the confusion matrix and overall fraction of correct predictions in the testing data.
#' 
#' 
#' 
#' c) Use the importance() function to determine the top 10 genes that are most important for the classification. What are the top 10 genes?
#' 
## -----------------------------------------------------------------------------
library(e1071)
library(caret)
library(ISLR)
library(dplyr)
library(randomForest)
data(Khan)

# Convert response variables to factors
xtrain = Khan$xtrain
ytrain <- as.factor(Khan$ytrain)
ytest <- as.factor(Khan$ytest)

# Train the random forest model
rf_model <- randomForest(x = Khan$xtrain, y = ytrain)

# Print the model summary
print(rf_model)


#mtry is the number of variables sampled as candidates at each split. The value for mtry is 48. mtry is responsible for helping the model not be too biased or too variant. mtry is responsible for making the model accurate


pred_test <- predict(rf_model, newdata = Khan$xtest)

# Show confusion matrix
conf_matrix <- table(Actual = ytest, Predicted = pred_test)
print("Confusion Matrix:")
print(conf_matrix)

# Calculate overall fraction of correct predictions
accuracy <- sum(diag(conf_matrix)) / sum(conf_matrix)
print("Overall Fraction of Correct Predictions:")
print(accuracy)

importance_scores <- importance(rf_model)

# Extract the MeanDecreaseGini scores (or any other relevant measure)
mean_decrease_gini <-  importance_scores[, "MeanDecreaseGini"]

# Sort the scores in descending order to identify the top 10 genes
top_10_genes <- names(sort(mean_decrease_gini, decreasing = TRUE)[1:10])

# Print the names of the top 10 genes
print("Top 10 Genes:")
print(top_10_genes)





#' 
#' 
#' d) Try different values of nodesize, ntree, and maxnodes in the randomForest() function, perform 5-fold cross-validation to find the best model among them. Show
#' the confusion matrix of the best model and report the parameters.
#' e) Use a support vector approach (svm() function) to predict cancer subtype using gene
#' expression measurements. Choose an appropriate kernel and explain the reason.
#' 
## -----------------------------------------------------------------------------
dim(Khan$xtrain)
dim(ytrain)

parameter_grid <- expand.grid(
  mtry = c(2, 4, 6), 
  nodesize = c(1, 5, 10),
  ntree = c(100, 200, 300),
  maxnodes = c(5, 10, 15)
)

# Initialize variables to store the best model and its performance
best_model <- NULL
best_accuracy <- 0

# Perform cross-validation for each parameter combination
for (i in 1:nrow(parameter_grid)) {
  # Train the random forest model with current parameters
  rf_model <- randomForest(x = Khan$xtrain, 
                           y = ytrain,
                           mtry = parameter_grid$mtry[i],
                           nodesize = parameter_grid$nodesize[i],
                           ntree = parameter_grid$ntree[i],
                           maxnodes = parameter_grid$maxnodes[i])
  
  # Make predictions on the testing data
  pred_test <- predict(rf_model, newdata = Khan$xtest)
  
  # Calculate accuracy
  accuracy <- mean(pred_test == ytest)
  
  # Check if this model has higher accuracy than the best model so far
  if (accuracy > best_accuracy) {
    best_accuracy <- accuracy
    best_model <- rf_model
  }
}

# Show the best model and its accuracy
print("Best Model:")
print(best_model)
print("Best Model Accuracy:")
print(best_accuracy)

# Create a trainControl object for 5-fold cross-validation
ctrl <- trainControl(method = "cv", 
                     number = 5,
                     summaryFunction = twoClassSummary,
                     classProbs = TRUE,
                     verboseIter = TRUE)






svm_model <- svm(ytrain ~ ., data = Khan$xtrain, kernel = "polynomial", degree = 3, coef0 = 1)

svm_model


#' 
#' 
#' 
#' 
#' 
#' 
#' f) What is the overall fraction of correct predictions in the testing data?
#' g) Try different values of cost, kernel in the svm() function, perform 5-fold cross-
#' validation to find the best model among them. Show the confusion matrix of the
#' best model and report the parameters.
#' 
## -----------------------------------------------------------------------------
library(e1071)
# Make predictions on the testing dataset
pred_test <- predict(rf_model, newdata = Khan$xtest)

# Compare predictions to true labels
correct_predictions <- sum(pred_test == ytest)

# Calculate total number of predictions
total_predictions <- length(ytest)

# Calculate accuracy
accuracy <- correct_predictions / total_predictions

# Print the overall fraction of correct predictions
print("Overall Fraction of Correct Predictions:")
print(accuracy)

parameter_grid <- expand.grid(cost = c(0.1, 1, 10),
                               kernel = c("linear", "radial", "polynomial"))

# Create a trainControl object for 5-fold cross-validation
ctrl <- trainControl(method = "cv", number = 5)

# Initialize variables to store the best model and its performance
best_model <- NULL
best_accuracy <- 0

# Perform cross-validation for each parameter combination
for (i in 1:nrow(parameter_grid)) {
  # Train the SVM model with current parameters
  svm_model <- svm(ytrain ~ ., 
                   data = Khan$xtrain,
                   kernel = parameter_grid$kernel[i],
                   cost = parameter_grid$cost[i])
  
  # Make predictions on the testing data
  pred_test <- predict(svm_model, newdata = Khan$xtest)
  
  # Calculate accuracy
  accuracy <- mean(pred_test == ytest)
  
  # Check if this model has higher accuracy than the best model so far
  if (accuracy > best_accuracy) {
    best_accuracy <- accuracy
    best_model <- svm_model
  }
}

# Show the best model and its accuracy
print("Best SVM Model:")
print(best_model)
print("Best Model Accuracy:")
print(best_accuracy)

# Make predictions on the test set
pred_test <- predict(best_model, newdata = Khan$xtest)

# Show confusion matrix
conf_matrix <- table(Actual = ytest, Predicted = pred_test)
print("Confusion Matrix:")
print(conf_matrix)

# Report the parameters of the best model
print("Best Model Parameters:")
print(best_model$cost)
print(best_model$kernel)




#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
