# Install and load the bnlearn package
if (!require(bnlearn)) {
  install.packages("bnlearn")
}
library(bnlearn)

# Read the CSV file
Final_Eigengenes <- read.csv("D:/UST/Year 4/Fall/Artificial Intelligence Techniques/Project/Data/Final_Eigengenes.csv", stringsAsFactors = FALSE)

Final_Eigengenes = as.data.frame(Final_Eigengenes)

samples <- Final_Eigengenes[, 2]

label <- Final_Eigengenes[, ncol(Final_Eigengenes)]

data_selected <- Final_Eigengenes[, -c(1, 2, ncol(Final_Eigengenes))]

cat = discretize(data_selected, method = "hartemink", breaks = 3)

# Combine the dataframes back together
data_final <- cbind(samples, cat, label)

colnames(data_final)[1] = 'Sample Identifier'

colnames(data_final)[ncol(data_final)] = 'ICB Response'

write.table(data_final, "Data/hartmink_cat_eigengenes.tsv", sep='\t', row.names=FALSE, quote=FALSE)
