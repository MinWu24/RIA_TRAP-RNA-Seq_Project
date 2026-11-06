
# Load libraries
library(clusterProfiler)
library(org.Ce.eg.db)

# Complete GO enrichment workflow for C. elegans
# ---------------------------------------------------
# Step 1: Prepare your data
# Replace these with your actual gene lists
  my_genes_of_interest <- read.csv("list.csv",sep=" ",header = TRUE)  # Your list A
  my_background_genes <- read.csv("background.csv",sep=",",header = TRUE)  # Your list B

my_genes_of_interest <- my_genes_of_interest$GeneID
my_background_genes <- my_background_genes$Gene.ID
# -------------------------------------------------
# Step 1: Prepare your data
# Replace these with your actual gene lists

# Step 2: Convert IDs to Entrez format if needed
convert_to_entrez <- function(gene_ids, id_type = "WORMBASE") {
  mapped <- mapIds(org.Ce.eg.db,
                   keys = gene_ids,
                   column = "ENTREZID",
                   keytype = id_type,
                   multiVals = "first")
  return(na.omit(mapped))
}

genes_entrez <- convert_to_entrez(my_genes_of_interest)
background_entrez <- convert_to_entrez(my_background_genes)

# Step 3: Perform enrichment for all three GO ontologies
# Biological Process
ego_bp <- enrichGO(gene = genes_entrez,
                   universe = background_entrez,
                   OrgDb = org.Ce.eg.db,
                   ont = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.05,
                   readable = TRUE)

# Cellular Component
ego_cc <- enrichGO(gene = genes_entrez,
                   universe = background_entrez,
                   OrgDb = org.Ce.eg.db,
                   ont = "CC",
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.05,
                   readable = TRUE)

# Molecular Function
ego_mf <- enrichGO(gene = genes_entrez,
                   universe = background_entrez,
                   OrgDb = org.Ce.eg.db,
                   ont = "MF",
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.05,
                   readable = TRUE)

# Step 4: Combine results
results <- list(
  Biological_Process = as.data.frame(ego_bp),
  Cellular_Component = as.data.frame(ego_cc),
  Molecular_Function = as.data.frame(ego_mf)
)

# Step 5: Save all results
for (ont in names(results)) {
  if (nrow(results[[ont]]) > 0) {
    write.csv(results[[ont]], 
              file = paste0("GO_", gsub("_", "", ont), "_results.csv"), 
              row.names = FALSE)
  }
}
