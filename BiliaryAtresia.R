# Load necessary libraries
library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(viridis)
library(vegan)
library(limma)
library(clusterProfiler)
library(org.Hs.eg.db)
library(forcats)
set.seed(123)

# ======================== Data Loading ======================== #
# Read clinical information and proteome data
clinical_info <- read.csv("data/clinical_info.csv") 
proteome_data <- read.csv("data/proteome_data.csv", row.names = 1, check.names = FALSE)

# ======================== PCA Analysis ======================== #
p <- prcomp(t(proteome_data))
clinical_info <- clinical_info %>%
  mutate(PC1 = p$x[,1], PC2 = p$x[,2],
         Group = factor(Group, levels = c("A", "B", "C", "D")),
         Sex = factor(Sex, levels = c("Male", "Female")))

# PCA Plots (Sex, Age, Group, AST levels)
pca_plots <- list(
  ggplot(clinical_info, aes(PC1, PC2, fill = Sex)) +
    geom_point(shape = 21, size = 3) + theme_classic() +
    xlab("PC1 (34.5%)") + ylab("PC2 (14.1%)") +
    scale_fill_brewer(palette = "Set1", direction = -1) +
    geom_text_repel(aes(label = Name), size = 2, box.padding = 0.15),
  
  ggplot(clinical_info, aes(PC1, PC2, fill = Month)) +
    geom_point(shape = 21, size = 3) + theme_classic() +
    scale_fill_viridis(direction = -1, name = "Age (month)") +
    xlab("PC1 (34.5%)") + ylab("PC2 (14.1%)"),
  
  ggplot(clinical_info, aes(PC1, PC2, fill = Group)) +
    geom_point(shape = 21, size = 3) + theme_classic() +
    scale_fill_manual(values = c("#7CAE00", "#F8766D", "#C77CFF", "#00BFC4")) +
    xlab("PC1 (34.5%)") + ylab("PC2 (14.1%)"),
  
  ggplot(clinical_info, aes(PC1, PC2, fill = AST)) +
    geom_point(shape = 21, size = 3) + theme_classic() +
    scale_fill_viridis(option = "B", direction = -1, name = "AST (U/L)") +
    xlab("PC1 (34.5%)") + ylab("PC2 (14.1%)")
)

# ======================== Statistical Test (PERMANOVA) ======================== #
dist_matrix <- dist(p$x)
adonis2(dist_matrix ~ Sex + Month + Group + AST, data = clinical_info, method = "euclidean", by = "margin")

# ======================== Differential Expression Analysis ======================== #
design <- model.matrix(~ Sex + Month + Group, data = clinical_info)

# Define contrast groups
groups <- c("B", "C")
contrasts <- makeContrasts(
  contrasts = setNames(paste0("Group", groups), paste0("Group", groups, "_vs_GroupA")), 
  levels = design
)

# limma analysis
fit <- lmFit(proteome_data, design) %>%
  contrasts.fit(contrasts) %>%
  eBayes()

# ======================== Volcano Plot & GSEA Analysis ======================== #
plot_list <- list()  # Store volcano plots
gsea_plot_list <- list()  # Store GSEA results

for (group in groups) {
  contrast_name <- paste0("Group", group)
  
  # Differential Expression Analysis
  results <- topTable(fit, coef = contrast_name, adjust.method = "BH", sort.by = "P", n = Inf) %>%
    arrange(P.Value)
  
  # Select Top Up/Downregulated Genes for Labeling
  selected_genes <- results %>%
    filter(logFC > 0) %>% arrange(P.Value) %>% head(6) %>%
    bind_rows(results %>% filter(logFC < 0) %>% arrange(P.Value) %>% head(6)) %>%
    rownames()
  
  results$label <- ifelse(rownames(results) %in% selected_genes, rownames(results), "")
  
  # Volcano Plot
  plot_list[[contrast_name]] <- ggplot(results) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey")+
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey")+
    geom_point(aes(logFC, -log10(P.Value), fill = -log10(P.Value) * sign(logFC), size = (-log10(P.Value))^2),
               shape = 21, color = "black", stroke = 0.3) +
    geom_text_repel(aes(logFC, -log10(P.Value), label = label), size = 3, max.overlaps = 99) +
    scale_size(range = c(0.1, 4)) +
    scale_fill_gradientn(colors = c("blue", "white", "red"), breaks = c(-20, 0, 20), limits = c(-20, 20), oob = scales::squish) +
    theme_classic(base_size = 16) +
    ggtitle(paste("Group", group, "vs Group A")) +
    xlab("logFC") + ylab(bquote(~-log[10] ~ italic(P) ~ "-value")) +
    theme(aspect.ratio=1,legend.position="none", plot.title = element_text(face="bold", hjust=0.5))
  
  # ======================== GSEA Analysis ======================== #
  results$score <- -log10(results$P.Value) * sign(results$logFC)
  
  rnk <- results %>%
    rownames_to_column("Gene") %>%
    dplyr::select(Gene, score) %>%
    dplyr::arrange(desc(score))
  
  # Convert Gene Symbols to ENTREZ IDs
  entrez_ID <- bitr(rnk$Gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db", drop = TRUE) %>% distinct(SYMBOL, .keep_all = TRUE) 
  rnk.set <- merge(rnk, entrez_ID, by.x = "Gene", by.y = "SYMBOL") %>% dplyr::arrange(desc(score))
  ranks <- setNames(rnk.set$score, rnk.set$ENTREZID)
  
  # Run GSEA
  ego <- gseGO(
    geneList = ranks,
    OrgDb = org.Hs.eg.db,
    ont = "BP",
    minGSSize = 10,
    maxGSSize = 500,
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    by = "fgsea",
    eps = 0
  )
  
  # Extract Top Pathways
  topPathways <- rbind(
    ego@result %>% filter(p.adjust < 0.05 & NES > 0) %>% dplyr::arrange(pvalue) %>% head(10),
    ego@result %>% filter(p.adjust < 0.05 & NES < 0) %>% dplyr::arrange(pvalue) %>% head(10)
  ) %>%
    mutate(metric = sign(NES) * -log10(p.adjust))
  
  # GSEA Bar Plot
  gsea_plot_list[[contrast_name]] <- ggplot(topPathways, aes(y = fct_reorder(Description, metric), x = metric)) +
    geom_bar(stat = "identity", aes(fill = NES)) +
    theme_bw() +
    scale_fill_gradient2(high = "#832424", low = "#3A3A98") +
    theme(panel.grid.major.y = element_blank(), axis.ticks.y = element_blank(), plot.title = element_text(face = "bold", hjust = 0.5)) +
    ggtitle(paste("Group", group, "vs Group A")) +
    xlab(bquote(~-log[10] ~ "adjusted" ~ italic(P) ~ "-value")) +
    ylab("") +
    geom_vline(xintercept = 0)
}

# ======================== Plot Outputs ======================== #
pca_plots

plot_list[["GroupB"]]
plot_list[["GroupC"]]

gsea_plot_list[["GroupB"]]
gsea_plot_list[["GroupC"]]
