#Mouse DESeq2
library(DESeq2)
library(biomaRt)
library(dplyr)

ms_read <- read.csv("salmon.count_LDL_Ms.csv", row.names = 1, check.names = FALSE)
samples <- c("Ctrl_1", "Ctrl_2", "Ctrl_3", "WT_1", "WT_2", "WT_3")
ms_read <- ms_read[, samples, drop = FALSE]
colnames(ms_read) <- c("Vehicle_1", "Vehicle_2", "Vehicle_3", "LDL_1", "LDL_2", "LDL_3")
sample_info <- data.frame(row.names = colnames(ms_read), condition = factor(c("Vehicle", "Vehicle", "Vehicle", "LDL", "LDL", "LDL"), levels = c("Vehicle", "LDL")))
dds <- DESeqDataSetFromMatrix(countData = round(as.matrix(ms_read)), colData = sample_info, design = ~ condition)
dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", "LDL", "Vehicle"))
res_df <- as.data.frame(res)
res_df$ensembl_gene_id <- rownames(res_df)
mart <- useEnsembl(biomart = "genes", dataset = "mmusculus_gene_ensembl")
anno <- getBM(attributes = c("ensembl_gene_id", "mgi_symbol"), filters = "ensembl_gene_id", values = res_df$ensembl_gene_id, mart = mart)
anno <- anno %>%
  filter(mgi_symbol != "") %>%
  distinct(ensembl_gene_id, .keep_all = TRUE)
res_df <- res_df %>%
  left_join(anno, by = "ensembl_gene_id")
deseq_mouse <- res_df %>%
filter(
!is.na(mgi_symbol),
mgi_symbol != ""
) %>%
group_by(mgi_symbol) %>%
slice_min(order_by = padj,n = 1,with_ties = FALSE
) %>%
ungroup()
write.csv(deseq_mouse,"mouse_DESeq2_LDL_vs_Vehicle_symbol.csv", row.names = FALSE)

##Comparison of human and mouse RNA-seq
library(biomaRt)
library(dplyr)
library(biomaRt)
library(ggplot2)
library(pheatmap)
deseq_human <- read.csv("human_DESeq2_LDL_vs_Vehicle.csv", check.names = FALSE)
human_all <- deseq_human %>%
  dplyr::filter(
    !is.na(hgnc_symbol),
    hgnc_symbol != ""
  ) %>%
  dplyr::arrange(
    hgnc_symbol,
    padj,
    pvalue,
    dplyr::desc(baseMean)
  ) %>%
  dplyr::group_by(hgnc_symbol) %>%
  dplyr::slice_head(n = 1) %>%
  dplyr::ungroup() %>%
  dplyr::transmute(
    human_symbol = hgnc_symbol,
    human_ensembl_id = ensembl_gene_id,
    baseMean_human = baseMean,
    log2FoldChange_human = log2FoldChange,
    pvalue_human = pvalue,
    padj_human = padj
  )
mouse_all <- deseq_mouse %>%
  dplyr::filter(
    !is.na(mgi_symbol),
    mgi_symbol != ""
  ) %>%
  dplyr::arrange(
    mgi_symbol,
    padj,
    pvalue,
    dplyr::desc(baseMean)
  ) %>%
  dplyr::group_by(mgi_symbol) %>%
  dplyr::slice_head(n = 1) %>%
  dplyr::ungroup() %>%
  dplyr::transmute(
    mouse_symbol = mgi_symbol,
    mouse_ensembl_id = ensembl_gene_id,
    baseMean_mouse = baseMean,
    log2FoldChange_mouse = log2FoldChange,
    pvalue_mouse = pvalue,
    padj_mouse = padj)
ortholog_all <- getLDS(
  attributes = c("mgi_symbol", "ensembl_gene_id"),
  filters = "mgi_symbol",
  values = unique(deseq_mouse$mgi_symbol),
  mart = mouse_mart,
  attributesL = c("hgnc_symbol", "ensembl_gene_id"),
  martL = human_mart,
  uniqueRows = TRUE
)
colnames(ortholog_all) <- c(
  "mouse_symbol",
  "mouse_ensembl_id",
  "human_symbol",
  "human_ensembl_id"
)

ortholog_1to1 <- ortholog_all %>%
  filter(
    !is.na(mouse_symbol),
    mouse_symbol != "",
    !is.na(human_symbol),
    human_symbol != ""
  ) %>%
  distinct() %>%
  group_by(mouse_symbol) %>%
  filter(n_distinct(human_symbol) == 1) %>%
  ungroup() %>%
  group_by(human_symbol) %>%
  filter(n_distinct(mouse_symbol) == 1) %>%
  ungroup()
plot_df <- ortholog_1to1 %>%
  dplyr::select(
    mouse_symbol,
    mouse_ensembl_id,
    human_symbol,
    human_ensembl_id
  ) %>%
  dplyr::distinct() %>%
  dplyr::inner_join(
    mouse_all,
    by = c(
      "mouse_symbol" = "mouse_symbol",
      "mouse_ensembl_id" = "mouse_ensembl_id"
    )
  ) %>%
  dplyr::inner_join(
    human_all,
    by = c(
      "human_symbol" = "human_symbol",
      "human_ensembl_id" = "human_ensembl_id"
    )
  ) %>%
  dplyr::filter(
    is.finite(log2FoldChange_mouse),
    is.finite(log2FoldChange_human)
  )

plot_df <- plot_df %>%
  dplyr::mutate(
    Category = dplyr::case_when(
      !is.na(pvalue_mouse) &
        !is.na(pvalue_human) &
        pvalue_mouse < 0.05 &
        pvalue_human < 0.05 &
        log2FoldChange_mouse > 1 &
        log2FoldChange_human > 1 ~ "Up in both",

      !is.na(pvalue_mouse) &
        !is.na(pvalue_human) &
        pvalue_mouse < 0.05 &
        pvalue_human < 0.05 &
        log2FoldChange_mouse < -1 &
        log2FoldChange_human < -1 ~ "Down in both",

!is.na(pvalue_mouse) &
!is.na(pvalue_human) &
pvalue_mouse < 0.05 &
pvalue_human < 0.05 &
abs(log2FoldChange_mouse) > 1 &
abs(log2FoldChange_human) > 1 &
log2FoldChange_mouse * log2FoldChange_human < 0 ~ "Opposite",

      TRUE ~ "Other"
    )
  ) %>%
  dplyr::mutate(
    Category = factor(
      Category,
      levels = c(
        "Other",
        "Opposite",
        "Up in both",
        "Down in both"
      )
    )
  )

plot_df %>%
  dplyr::count(Category)
plot_df %>%
  dplyr::filter(Category %in% c("Up in both", "Down in both")) %>%
  dplyr::count(Category)

ggplot(
  plot_df,
  aes(
    x = log2FoldChange_mouse,
    y = log2FoldChange_human
  )
) +
  geom_point(
    data = plot_df %>%
      dplyr::filter(Category == "Other"),
    aes(color = Category),
    size = 2.5,
    alpha = 0.5
  ) +
  geom_point(
    data = plot_df %>%
      dplyr::filter(Category != "Other"),
    aes(color = Category),
    size = 3.5,
    alpha = 0.9
  ) +
  geom_vline(
    xintercept = c(-1,1),
    linetype = "dashed",
    linewidth = 0.7,
    color = "black"
  ) +
  geom_hline(
    yintercept = c(-1,1),
    linetype = "dashed",
    linewidth = 0.7,
    color = "black"
  ) +
  scale_color_manual(
    values = c(
      "Other" = "grey75",
      "Opposite" = "black",
      "Up in both" = "#D32F2F",
      "Down in both" = "#1976D2"
    ),
    labels = c(
      "Other" = "Other genes",
      "Opposite" = "Opposite direction",
      "Up in both" = "Up in both",
      "Down in both" = "Down in both"
    )
  ) +

  labs(
    x = expression(Mouse~log[2]~fold~change),
    y = expression(Human~log[2]~fold~change),
    color = NULL
  ) +

  theme_classic(base_size = 16) +

  theme(
    axis.title = element_text(size = 35),
    axis.text = element_text(size = 30),
    legend.text = element_text(size = 25),
    legend.position = "right"
  )

heat_df_2 <- human_gene_tpm %>%
  dplyr::filter(hgnc_symbol %in% genes_2) %>%
  dplyr::select(
    hgnc_symbol,
    Vehicle_1,
    Vehicle_2,
    Vehicle_3,
    LDL_1,
    LDL_2,
    LDL_3
  )
setdiff(genes_2, heat_df_2$hgnc_symbol)
heat_mat_2 <- heat_df_2 %>%
  as.data.frame()
rownames(heat_mat_2) <- heat_mat_2$hgnc_symbol
heat_mat_2$hgnc_symbol <- NULL
heat_mat_2 <- as.matrix(heat_mat_2)
heat_log_2 <- log2(heat_mat_2 + 1)
heat_z_2 <- t(scale(t(heat_log_2)))
heat_z_2 <- heat_z_2[
  apply(heat_z_2, 1, function(x) all(is.finite(x))),
  ,
  drop = FALSE
]
heat_z_2[heat_z_2 > 3] <- 3
heat_z_2[heat_z_2 < -3] <- -3
gene_order_2 <- genes_2[
  genes_2 %in% rownames(heat_z_2)
]
heat_z_2 <- heat_z_2[
  gene_order_2,
  ,
  drop = FALSE

anno_col_2 <- data.frame(
  Group = factor(
    c(
      "Vehicle", "Vehicle", "Vehicle",
      "LDL", "LDL", "LDL"
    ),
    levels = c("Vehicle", "LDL")
  )
)
rownames(anno_col_2) <- colnames(heat_z_2)

ann_colors_2 <- list(
  Group = c(
    Vehicle = "white",
    LDL = "black"
  )
)
pheatmap(
  heat_z_2,
  color = colorRampPalette(
    c("blue", "white", "red")
  )(100),
  breaks = seq(
    -3,
    3,
    length.out = 101
  ),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  annotation_col = anno_col_2,
  annotation_colors = ann_colors_2,
  show_colnames = FALSE,
  show_rownames = TRUE,
  annotation_names_col = FALSE,
  annotation_legend = FALSE,
  border_color = "white",
  cellwidth = 22,
  cellheight = 22,
  fontsize = 14,
  fontsize_row = 15,
  main = "Wound healing genes"
)

plot_df <- deseq_gene %>%
  dplyr::filter(
    !is.na(padj),
    !is.na(log2FoldChange)
  ) %>%
  dplyr::mutate(
    padj_plot = pmax(padj, .Machine$double.xmin),
    logP = -log10(padj_plot),

    significance = dplyr::case_when(
      padj < 0.05 & log2FoldChange > 1  ~ "Up",
      padj < 0.05 & log2FoldChange < -1 ~ "Down",
      TRUE                              ~ "NS"
    )
  )
volcano_plot <- ggplot(
  plot_df,
  aes(x = log2FoldChange, y = logP)
) +
  geom_point(
    aes(color = significance),
    size = 2,
    alpha = 0.5
  ) +
  geom_vline(
    xintercept = c(-1, 1),
    linetype = "dashed",
    linewidth = 0.6
  ) +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed",
    linewidth = 0.6
  ) +
  geom_text_repel(
    data = top10_up,
    aes(label = hgnc_symbol),
    color = "red",
    size = 4,
    fontface = "bold",
    max.overlaps = Inf,
    box.padding = 0.5,
    point.padding = 0.3,
    min.segment.length = 0,
    segment.color = "#D55E00",
    show.legend = FALSE
  ) +
  geom_text_repel(
    data = top10_down,
    aes(label = hgnc_symbol),
    color = "blue",
    size = 4,
    fontface = "bold",
    max.overlaps = Inf,
    box.padding = 0.5,
    point.padding = 0.3,
    min.segment.length = 0,
    segment.color = "#0072B2",
    show.legend = FALSE
  ) +
  scale_color_manual(
    values = c(
      Up   = "red",
      Down = "blue",
      NS   = "grey75"
    ),
    breaks = c("Up", "Down", "NS"),
    labels = c(
      Up = "Upregulated",
      Down = "Downregulated",
      NS = "Not significant"
    )
  ) +
  labs(
    x = expression(log[2] ~ fold ~ change),
    y = expression(-log[10] ~ adjusted ~ italic(P)),
    color = NULL
  ) +
  theme_classic(base_size = 15) +
  theme(
    axis.title = element_text(size = 17),
    axis.text = element_text(size = 14),
    legend.text = element_text(size = 13),
    legend.position = "right"
  )
volcano_plot