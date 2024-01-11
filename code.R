if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("graph", "Rgraphviz", "RBGL", "GO.db", "impute", "preprocessCore", "bnlearn"))
install.packages('tidyverse')
install.packages("WGCNA")
install.packages("gRain")
install.packages("caret")
library(tidyverse)     
library(magrittr)     
library(WGCNA)  
library(readr)
library(dplyr)
library(bnlearn)
library(graph)
library(Rgraphviz)
library(RBGL)
library(gRain)
library(ggplot2)
library(caret)



# Read the data files
data1 <- readr::read_delim("C:/Users/ALY HISSAM/Downloads/mrna_cont.tsv", delim = "\t")
data2 <- readr::read_delim("C:/Users/ALY HISSAM/Downloads/labels_df.tsv", delim = "\t")
# Merge the data frames based on the Sample Identifier
merged_data <- merge(data1, data2, by = "Sample Identifier") # Including sample name and ICB response.
# Filter the merged data frame to include only responders (ICB Response = 1)
responders_data <- filter(merged_data, `ICB Response` == 1)
responders_data  # with sample name and response
selected_columns <- responders_data[, 3:ncol(responders_data)-1]
selected_columns # without sample name or ICB response


allowWGCNAThreads()  
powers = c(c(1:10), seq(from = 12, to = 20, by = 2)) # vector of 'powers' ranging from 1 to 10 and then from 12 to 20 in steps of 2.
sft = pickSoftThreshold(selected_columns, powerVector = powers, verbose = 5)
par(mfrow = c(1,2));  # Set the plotting area to 1 row and 2 columns
cex1 = 0.9;  # Set the character expansion size for plot text
# Scale Independence Plot
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2], xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit, signed R^2", main = "Scale independence")
text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2], labels = powers, cex = cex1, col = "red")
abline(h = 0.90, col = "red")  # Add a horizontal line at R^2 = 0.90
# Mean Connectivity Plot
plot(sft$fitIndices[, 1], sft$fitIndices[, 5], xlab = "Soft Threshold (power)", ylab = "Mean Connectivity", type = "n", main = "Mean connectivity")
text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers, cex = cex1, col = "red")




picked_power = 6
# Resolve namespace conflict by temporarily overriding the default cor function.
temp_cor <- cor
cor <- WGCNA::cor  # Force it to use WGCNA cor function
# Construct the network using block wise Modules function
netwk <- blockwiseModules(
  selected_columns,  # Input data matrix (gene expression data)
  power = picked_power,  # Soft-threshold power
  networkType = "signed",  # Type of network (signed or unsigned)
  TOMType = "unsigned",
  maxBlockSize = 16067
)
# Restore the original cor function to its namespace
cor <- temp_cor

#####
mergedColors = labels2colors(netwk$colors)
# Plot the dendrogram with module colors
plotDendroAndColors(
  netwk$dendrograms[[1]],  # The first dendrogram in the network analysis results
  mergedColors[netwk$blockGenes[[1]]],  # Colors corresponding to the first block of genes
  "Module colours",  # Title of the plot
  dendroLabels = FALSE,  # Whether to show labels on the dendrogram
  hang = 0.03,  # Adjusts how the labels and tips of the dendrogram hang
  addGuide = TRUE,  # Add a color guide
  guideHang = 0.05  # Adjusts the position of the color guide
)

# Create a data frame with gene identifiers and their corresponding module colors
module_df <- data.frame(
  gene_id = names(netwk$colors),  # Extract gene names from the network object
  colors = labels2colors(netwk$colors)  # Convert module labels to colors
)
# Save module_df to CSV file
write.csv(module_df, file = "C:/Users/ALY HISSAM/OneDrive/Desktop/AI - Project/Gene_Color_Coded.csv", row.names = TRUE)


merg<- merged_data[, 3:ncol(merged_data)-1]
sam <-merged_data[, 1]
MEs0 <- moduleEigengenes(merg, mergedColors)
eigengenes_with_sample_id <- cbind("Sample Identifier" = merged_data[, 1], MEs0$eigengenes)
Final_Eigengenes <- merge(eigengenes_with_sample_id, data2, by = "Sample Identifier")
write.csv(Final_Eigengenes, file = "Final_Eigengenes_16067max.csv", row.names = TRUE)

# Final_Eigengenes <- readr::read_delim("C:/Users/ALY HISSAM/OneDrive/Desktop/AI - Project/Final_Eigengenes_16067max.csv", delim = ",")
# Final_Eigengenes<- Final_Eigengenes[, 2:ncol(Final_Eigengenes)]
# Discretizing Eigengenes.
Discretized_Eigengenes<-discretize(Final_Eigengenes[,3:ncol(Final_Eigengenes)-1], method = 'hartemink', breaks = 3, ibreaks = 60, idisc = "quantile")
# Adding 'Sample Identifier' and 'ICB Response' to Final_Discretized_Eigengenes.
Final_Discretized_Eigengenes <- Discretized_Eigengenes
Final_Discretized_Eigengenes$'ICB Response' <- Final_Eigengenes$'ICB Response'
Final_Discretized_Eigengenes$'ICB Response' <- factor(Final_Discretized_Eigengenes$'ICB Response')
Final_Discretized_Eigengenes_with_Samples <- Final_Discretized_Eigengenes %>%
  mutate('Sample Identifier' = Final_Eigengenes$'Sample Identifier') %>%
  select('Sample Identifier', everything())

# 50 genes
bestgenes <- readr::read_delim("C:/Users/ALY HISSAM/OneDrive/Desktop/AI - Project/mrna_categorized_labelled_filtered_stringent.tsv", delim = "\t")
b<-bestgenes[,2:ncol(bestgenes)]
b_discrete_factors <- data.frame(lapply(b, as.factor))
# First, identify the position of the column 'ICB.Response'
col_index <- which(names(b_discrete_factors) == "ICB.Response")
# Now, rename the column
names(b_discrete_factors)[col_index] <- "ICB Response"


b_discrete_factors <- Final_Discretized_Eigengenes
# Final_Eigengenes = bestgenes
bl = data.frame(from = "ICB Response", to = colnames(b_discrete_factors))
bl = bl[-which(bl$to == "ICB Response"),]
dag2 <- boot.strength(b_discrete_factors, R = 50, algorithm = "hc", algorithm.args = list(score = "bde", iss = 10, blacklist=bl))
dag2
Consensus_Network <- averaged.network(dag2, threshold = 0.3)
Consensus_Network
undirected.arcs(Consensus_Network)
has_cycle <- function(net, from, to) {
  test_net <- set.arc(net, from, to, check.cycles = FALSE)
  !acyclic(test_net)
}
undirected_arcs <- undirected.arcs(Consensus_Network)
for (i in 1:nrow(undirected_arcs)) {
  from <- undirected_arcs[i, "from"]
  to <- undirected_arcs[i, "to"]
    if (!has_cycle(Consensus_Network, from, to)) {
      Consensus_Network <- set.arc(Consensus_Network, from, to)
  } else if (!has_cycle(Consensus_Network, to, from)) {
    Consensus_Network <- set.arc(Consensus_Network, to, from)
  }
}
Consensus_Network
score(Consensus_Network, b_discrete_factors, type = "bde", iss = 10)
graphviz.plot(Consensus_Network)
# Assuming your data is in a data frame called 'data'
set.seed(123) # For reproducibility
sample_size <- nrow(b_discrete_factors)
train_indices <- sample(sample_size, sample_size * 0.8) # 80% for training
test_indices <- setdiff(1:sample_size, train_indices)   # Remaining for testing
# Splitting data into training and testing sets
train_data <- b_discrete_factors[train_indices, ]
test_data <- b_discrete_factors[test_indices, ]
# Fit Bayesian network model on the training set
fitted = bn.fit(Consensus_Network, train_data, method = "bayes")
# Predict on the test set
predicted_values <- predict(fitted, node = "ICB Response", data = test_data, method = "bayes-lw", prob = TRUE)
predicted_values
# Extract probabilities for class '1'
prob_class_1 <- attr(predicted_values, "prob")[2, ]
# Convert these probabilities to predicted labels based on a threshold (e.g., 0.5)
predicted_labels <- ifelse(prob_class_1 > 0.5, 1, 0)
predicted_labels
# Actual values - ensuring they are numeric
actual_values <- as.numeric(as.character(test_data$`ICB Response`))
actual_values
# For F1 Score
confusion_matrix <- confusionMatrix(as.factor(predicted_labels), as.factor(actual_values))
f1_score <- confusion_matrix$byClass['F1']
print(paste("F1 Score:", f1_score))
# For Accuracy
accuracy <- sum(predicted_labels == actual_values, na.rm = TRUE) / sum(!is.na(predicted_labels))
print(paste("Accuracy:", accuracy))
# Create the confusion matrix
cm <- confusionMatrix(as.factor(predicted_labels), as.factor(actual_values))
# Convert the confusion matrix to a data frame for ggplot2
cm_table <- as.data.frame.matrix(cm$table)
# Add a 'Prediction' column to the data frame
cm_table$Prediction <- rownames(cm_table)
cm_table_long <- reshape2::melt(cm_table, id.vars = 'Prediction')
# Rename the columns appropriately
names(cm_table_long) <- c("Prediction", "Reference", "Freq")
# Plotting the confusion matrix using ggplot2
ggplot(data = cm_table_long, aes(x = Reference, y = Prediction, fill = Freq)) +
  geom_tile() +
  geom_text(aes(label = sprintf("%0.0f", Freq)), vjust = 1) +
  scale_fill_gradient(low = "white", high = "steelblue") +
  theme_minimal() +
  labs(fill = "Frequency", title = "Confusion Matrix")

