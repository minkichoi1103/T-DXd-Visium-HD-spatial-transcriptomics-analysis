library(decoupleR)
library(limma)
library(tidyverse)


## progeny 
net <- decoupleR::get_progeny(organism = 'human', 
                              top = 500)

net


expr <- read.csv('../data/pseudobulk_cancer.csv', row.names = 1)
sample_acts <- decoupleR::run_mlm(mat = expr, 
                                  net = net, 
                                  .source = 'source', 
                                  .target = 'target',
                                  .mor = 'weight', 
                                  minsize = 5)
sample_acts

write_csv(sample_acts, '../data/pathway_activity_decoupler.csv')


##HER2 positive

### DEG analysis : limma ###
expr <- read.csv('../data/pseudobulk_cancer.csv', row.names = 1)
metadata <- read.csv('../data/results_df_spatial_vessel_linkerenz.csv')
metadata <- metadata[,c('batch', 'clinical_type', 'HER2_type', 'sample_id')]
rownames(metadata) <- metadata$batch
metadata <- metadata[colnames(expr), ]


metadata <- metadata[metadata$HER2_type == 'Pos', ]
expr <- expr[, colnames(expr) %in% rownames(metadata)]


# expr: log-transformed expression matrix (genes x samples)
# design: model.matrix(~clinical_type)

# 1) create the design matrix
design <- model.matrix(~ 0 + clinical_type, data=metadata)
colnames(design) <- c("clinical_type_resistant", "clinical_type_sensitive")

# 2) estimate the inter-sample correlation using batch (= patient) information
#    the block argument is the patient ID vector (i.e., metadata$batch)
corfit <- duplicateCorrelation(expr, design=design, block=metadata$sample_id)
corfit$consensus

fit <- lmFit(expr, design, block=metadata$sample_id, correlation=corfit$consensus)
fit <- eBayes(fit)

# contrast for clinical_type_resistant vs. clinical_type_sensitive
contrast.matrix <- makeContrasts(
  resistant_vs_sensitive = clinical_type_resistant - clinical_type_sensitive,
  levels = design
)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

results <- topTable(fit2, coef="resistant_vs_sensitive", n=Inf)

results %>%
  rownames_to_column(var = "GeneID") %>%  # "GeneID" : put your desired column name here
  write_csv('../data/DEG_results_HER2_positive.csv')

# genes with logFC > 0 in 'results' are more highly expressed in resistant


# Run mlm (use limma results)
deg <- results$t %>% as.data.frame()  ## Result based on t-values from R limma; reflects the model structure used in limma.  
rownames(deg) <- rownames(results)
contrast_acts <- decoupleR::run_mlm(mat  =deg, 
                                    net = net, 
                                    .source = 'source', 
                                    .target = 'target',
                                    .mor = 'weight', 
                                    minsize = 5)
contrast_acts




# Bar Plot about Pathway activity 
colors <- rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")[c(2, 10)])
title <- "Pathways score : HER2 positive"
p <- ggplot2::ggplot(data = contrast_acts, 
                     mapping = ggplot2::aes(x = stats::reorder(source, score), 
                                            y = score)) + 
  ggplot2::geom_bar(mapping = ggplot2::aes(fill = score),
                    color = "black",
                    stat = "identity") +
  ggplot2::scale_fill_gradient2(low = colors[1], 
                                mid = "whitesmoke", 
                                high = colors[2], 
                                midpoint = 0) + 
  ggplot2::theme_minimal() +
  ggplot2::theme(axis.title = element_text(face = "bold", size = 12),
                 axis.text.x = ggplot2::element_text(angle = 45, 
                                                     hjust = 1, 
                                                     size = 10, 
                                                     face = "bold"),
                 axis.text.y = ggplot2::element_text(size = 10, 
                                                     face = "bold"),
                 panel.grid.major = element_blank(), 
                 panel.grid.minor = element_blank()) +
  ggplot2::xlab("Pathways") + 
  ggplot2::ggtitle(title) + theme(plot.title = element_text(hjust = 0.5))
p




## Identify genes contributing to increased EGFR pathway activity. 
pathway <- 'EGFR'

df <- net %>%
  dplyr::filter(source == pathway) %>%
  dplyr::arrange(target) %>%
  dplyr::mutate(ID = target, 
                color = "3") %>%
  tibble::column_to_rownames('target')

inter <- sort(dplyr::intersect(rownames(deg), rownames(df)))

df <- df[inter, ]

df['t_value'] <- deg[inter, ]

df <- df %>%
  dplyr::mutate(color = dplyr::if_else(weight > 0 & t_value > 0, '1', color)) %>%
  dplyr::mutate(color = dplyr::if_else(weight > 0 & t_value < 0, '2', color)) %>%
  dplyr::mutate(color = dplyr::if_else(weight < 0 & t_value > 0, '2', color)) %>%
  dplyr::mutate(color = dplyr::if_else(weight < 0 & t_value < 0, '1', color))

colors <- rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")[c(2, 10)])

title <- sprintf('%s : HER2 positive', pathway)

p <- ggplot2::ggplot(data = df, 
                     mapping = ggplot2::aes(x = weight, 
                                            y = t_value, 
                                            color = color)) + 
  ggplot2::geom_point(size = 2.5, 
                      color = "black") + 
  ggplot2::geom_point(size = 1.5) +
  ggplot2::scale_colour_manual(values = c(colors[2], colors[1], "grey")) +
  ggrepel::geom_label_repel(mapping = ggplot2::aes(label = ID)) + 
  ggplot2::theme_minimal() +
  ggplot2::theme(legend.position = "none") +
  ggplot2::geom_vline(xintercept = 0, linetype = 'dotted') +
  ggplot2::geom_hline(yintercept = 0, linetype = 'dotted') +
  ggplot2::ggtitle(title) + theme(plot.title = element_text(hjust = 0.5)) + 

# [fine-tune font sizes here]
ggplot2::theme(
  # 1. Title size and alignment (size: size, face: bold, hjust: center)
  plot.title = ggplot2::element_text(size = 20,  hjust = 0.5),
  
  # 2. Axis-label size (text size for x='weight', y='t_value')
  axis.title = ggplot2::element_text(size = 16),
  
  # 3. Axis-tick text size (the numbers on the axes)
  axis.text = ggplot2::element_text(size = 12)
)

ggsave(filename = "../figures/HER2_postive_EGFR.png", 
       plot = p, 
       dpi = 400, 
       width = 8, height = 6, units = "in", bg = "white") # adjust width and height as needed


# 1) compute padj & order the y-axis
df_plot <- contrast_acts %>%
  mutate(padj   = p.adjust(p_value, method = "BH"),
         source = fct_reorder(source, score))   # sort by score

# 2) Bubble plot
p <- ggplot(df_plot,
       aes(x = score,
           y = source,
           size = -log10(padj),    # point size: significance
           colour = score)) +      # color: score (positive/negative)
  geom_point() +
  scale_colour_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                         midpoint = 0, name = "score") +
  scale_size(range = c(4, 12), name = "P adj",
             guide = guide_legend(order = 1)) +
  labs(x = "score", y = NULL) +
  theme_classic(base_size = 20) +
  theme(legend.position = "right",
        axis.ticks.y = element_blank())

# 2. Render the plot (for inspection)
print(p)

# [change 2] Save the image (dpi = 400)
ggsave(filename = "../figures/HER2_postive_bubbleplot.png", 
       plot = p, 
       dpi = 400, 
       width = 8, height = 6, units = "in") # adjust width and height as needed






#### 2. HER2 low
### DEG analysis using limma ###
expr <- read.csv('../data/pseudobulk_cancer.csv', row.names = 1)
metadata <- read.csv('../data/results_df_spatial_vessel_linkerenz.csv')
metadata <- metadata[,c('batch', 'clinical_type', 'HER2_type', 'sample_id')]
rownames(metadata) <- metadata$batch
metadata <- metadata[colnames(expr), ]


metadata <- metadata[metadata$HER2_type == 'Low', ]
expr <- expr[, colnames(expr) %in% rownames(metadata)]


# expr: log-transformed expression matrix (genes x samples)
# design: model.matrix(~clinical_type)

# 1) create the design matrix
design <- model.matrix(~ 0 + clinical_type, data=metadata)
colnames(design) <- c("clinical_type_resistant", "clinical_type_sensitive")

# 2) estimate the inter-sample correlation using batch (= patient) information
#    the block argument is the patient ID vector (i.e., metadata$batch)
corfit <- duplicateCorrelation(expr, design=design, block=metadata$sample_id)
corfit$consensus

fit <- lmFit(expr, design, block=metadata$sample_id, correlation=corfit$consensus)
fit <- eBayes(fit)

# contrast for clinical_type_resistant vs. clinical_type_sensitive
contrast.matrix <- makeContrasts(
  resistant_vs_sensitive = clinical_type_resistant - clinical_type_sensitive,
  levels = design
)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

results <- topTable(fit2, coef="resistant_vs_sensitive", n=Inf)

results %>%
  rownames_to_column(var = "GeneID") %>%  # "GeneID" : put your desired column name here
  write_csv('../data/DEG_results_HER2_low.csv')

# genes with logFC > 0 in 'results' are more highly expressed in resistant


# Run mlm (use limma results)
deg <- results$t %>% as.data.frame()  ## Result based on t-values from R limma; reflects the model structure used in limma.  
rownames(deg) <- rownames(results)
contrast_acts <- decoupleR::run_mlm(mat  =deg, 
                                    net = net, 
                                    .source = 'source', 
                                    .target = 'target',
                                    .mor = 'weight', 
                                    minsize = 5)
contrast_acts



# Bar Plot about Pathway activity 
colors <- rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")[c(2, 10)])
title <- "Pathways score : HER2 low"
p <- ggplot2::ggplot(data = contrast_acts, 
                     mapping = ggplot2::aes(x = stats::reorder(source, score), 
                                            y = score)) + 
  ggplot2::geom_bar(mapping = ggplot2::aes(fill = score),
                    color = "black",
                    stat = "identity") +
  ggplot2::scale_fill_gradient2(low = colors[1], 
                                mid = "whitesmoke", 
                                high = colors[2], 
                                midpoint = 0) + 
  ggplot2::theme_minimal() +
  ggplot2::theme(axis.title = element_text(face = "bold", size = 12),
                 axis.text.x = ggplot2::element_text(angle = 45, 
                                                     hjust = 1, 
                                                     size = 10, 
                                                     face = "bold"),
                 axis.text.y = ggplot2::element_text(size = 10, 
                                                     face = "bold"),
                 panel.grid.major = element_blank(), 
                 panel.grid.minor = element_blank()) +
  ggplot2::xlab("Pathways") + 
  ggplot2::ggtitle(title) + theme(plot.title = element_text(hjust = 0.5))
p



## Identify genes contributing to increased JAK-STAT pathway activity. 
pathway <- 'JAK-STAT'

df <- net %>%
  dplyr::filter(source == pathway) %>%
  dplyr::arrange(target) %>%
  dplyr::mutate(ID = target, 
                color = "3") %>%
  tibble::column_to_rownames('target')

inter <- sort(dplyr::intersect(rownames(deg), rownames(df)))

df <- df[inter, ]

df['t_value'] <- deg[inter, ]

df <- df %>%
  dplyr::mutate(color = dplyr::if_else(weight > 0 & t_value > 0, '1', color)) %>%
  dplyr::mutate(color = dplyr::if_else(weight > 0 & t_value < 0, '2', color)) %>%
  dplyr::mutate(color = dplyr::if_else(weight < 0 & t_value > 0, '2', color)) %>%
  dplyr::mutate(color = dplyr::if_else(weight < 0 & t_value < 0, '1', color))

colors <- rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")[c(2, 10)])

title <- sprintf('%s : HER2_low', pathway)

p <- ggplot2::ggplot(data = df, 
                     mapping = ggplot2::aes(x = weight, 
                                            y = t_value, 
                                            color = color)) + 
  ggplot2::geom_point(size = 2.5, 
                      color = "black") + 
  ggplot2::geom_point(size = 1.5) +
  ggplot2::scale_colour_manual(values = c(colors[2], colors[1], "grey")) +
  ggrepel::geom_label_repel(mapping = ggplot2::aes(label = ID)) + 
  ggplot2::theme_minimal() +
  ggplot2::theme(legend.position = "none") +
  ggplot2::geom_vline(xintercept = 0, linetype = 'dotted') +
  ggplot2::geom_hline(yintercept = 0, linetype = 'dotted') +
  ggplot2::ggtitle(title) + theme(plot.title = element_text(hjust = 0.5))

p


# 1) compute padj & order the y-axis
df_plot <- contrast_acts %>%
  mutate(padj   = p.adjust(p_value, method = "BH"),
         source = fct_reorder(source, score))   # sort by score

# 2) Bubble plot
p <- ggplot(df_plot,
       aes(x = score,
           y = source,
           size = -log10(padj),    # point size: significance
           colour = score)) +      # color: score (positive/negative)
  geom_point() +
  scale_colour_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                         midpoint = 0, name = "score") +
  scale_size(range = c(4, 12), name = "P adj",
             guide = guide_legend(order = 1)) +
  labs(x = "score", y = NULL) +
  theme_classic(base_size = 20) +
  theme(legend.position = "right",
        axis.ticks.y = element_blank())

# [change 2] Save the image (dpi = 400)
ggsave(filename = "../figures/HER2_low_bubbleplot.png", 
       plot = p, 
       dpi = 400, 
       width = 8, height = 6, units = "in") # adjust width and height as needed


deg_her2_pos <- read.csv('../data/DEG_results_HER2_positive.csv')
deg_her2_low <- read.csv('../data/DEG_results_HER2_low.csv')

deg_her2_pos$'HER2 IHC' <- 'Positive'
deg_her2_low$'HER2 IHC' <- 'Low'

combined_deg <- bind_rows(deg_her2_pos, deg_her2_low)

combined_deg %>% write_csv('../data/Supplementary_table2.csv')
