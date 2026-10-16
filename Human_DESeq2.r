#Salmon on Mac terminal
for R1 in *_1.fastq.gz
do
    sample=${R1%_1.fastq.gz}
    R2="${sample}_2.fastq.gz"
    if [ ! -f "$R2" ]; then
        echo "Skipping ${sample}: R2 not found"
        continue
    fi
    echo "Processing ${sample}..."
    fastp \
        -i "$R1" \
        -I "$R2" \
        -3 \
        -o "out_${sample}_1.fastq.gz" \
        -O "out_${sample}_2.fastq.gz" \
        -h "report_${sample}.html" \
        -j "report_${sample}.json" \
        -q 15 \
        -n 10 \
        -t 1 \
        -T 1 \
        -l 20 \
        -w 16 \
        -f 1 \
        -F 1
done

wget ftp://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz
salmon index -t Homo_sapiens.GRCh38.cdna.all.fa.gz -i Homo_sapiens.GRCh38_salmon_index
INDEX_human="Homo_sapiens.GRCh38_salmon_index"
samples=(
Vehicle_1
Vehicle_2
Vehicle_3
LDL_1
LDL_2
LDL_3
)

for sample in "${samples[@]}"
do
    echo "Running ${sample}..."

    salmon quant \
        -i "$INDEX_human" \
        -l IU \
        -1 "out_${sample}_1.fastq.gz" \
        -2 "out_${sample}_2.fastq.gz" \
        --validateMappings
        -o "${sample}_salmon_quant"
done

#Tximport and DESeq2 on R
library(tximport)
library(DESeq2)
library(biomaRt)
library(dplyr)
library(readr)
salmon_dirs <- list.files(".", pattern = "_salmon_quant$", full.names = TRUE)
samples <- sub("_salmon_quant$", "", basename(salmon_dirs))
salmon.files <- file.path(salmon_dirs, "quant.sf")
names(salmon.files) <- samples
sample_order <- c("Vehicle_1", "Vehicle_2", "Vehicle_3", "LDL_1", "LDL_2", "LDL_3")
salmon.files <- salmon.files[sample_order]
stopifnot(!anyNA(salmon.files), all(file.exists(salmon.files)))
mart <- useEnsembl(biomart = "genes",dataset = "hsapiens_gene_ensembl",　version = 110)
tx2gene <- getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id"),mart = mart)
colnames(tx2gene) <- c("TXNAME", "GENEID")
txi <- tximport(salmon.files, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)
write.csv(txi$counts, "salmon.count_LDL_Hm.csv")
write.csv(txi$abundance, "salmon.tpm_LDL_Hm.csv")

sample_info <- data.frame(row.names = sample_order, condition = factor(c("Vehicle", "Vehicle", "Vehicle", "LDL", "LDL", "LDL"), levels = c("Vehicle", "LDL")))
dds <- DESeqDataSetFromTximport(txi = txi, colData = sample_info, design = ~ condition)
keep <- rowSums(counts(dds) >= 10) >= 3
dds <- DESeq(dds)
res <- results(dds,contrast = c("condition", "LDL", "Vehicle"))
res_df <- as.data.frame(res)
res_df$ensembl_gene_id <- rownames(res_df)
anno <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), filters = "ensembl_gene_id", values = res_df$ensembl_gene_id, mart = mart)
anno <- anno %>%
  filter(hgnc_symbol != "") %>%
  distinct(ensembl_gene_id, .keep_all = TRUE)
res_df <- res_df %>%
  left_join(anno, by = "ensembl_gene_id")
deseq_human <- res_df %>%
  mutate(DE = case_when(!is.na(padj) & padj < 0.05 & log2FoldChange > 1 ~ "Up", !is.na(padj) & padj < 0.05 & log2FoldChange < -1 ~ "Down", TRUE ~ "NC"))
write.csv(deseq_human, "human_DESeq2_LDL_vs_Vehicle.csv", row.names = FALSE)

#VolcanoPlot
library(ggplot2)
library(ggrepel)
plot_df <- deseq_human %>%
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
TRUE                              ~ "NS"))
top10_up <- plot_df %>%
dplyr::filter(
significance == "Up",
!is.na(hgnc_symbol),
hgnc_symbol != ""
) %>%
dplyr::arrange(padj) %>%
head(10)
top10_down <- plot_df %>%
dplyr::filter(
significance == "Down",
!is.na(hgnc_symbol),
hgnc_symbol != ""
) %>%
dplyr::arrange(padj) %>%
head(10)

volcano_plot <- ggplot(
  plot_df,
  aes(
    x = log2FoldChange,
    y = logP
  )
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
    segment.color = "red",
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
    segment.color = "blue",
    show.legend = FALSE
  ) +
  scale_color_manual(
    values = c(
      "Up" = "red",
      "Down" = "blue",
      "NS" = "grey75"
    ),
    breaks = c("Up", "Down", "NS"),
    labels = c(
      "Up" = "Upregulated",
      "Down" = "Downregulated",
      "NS" = "Not significant"
    ),
    drop = FALSE
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
    legend.position = "right")
volcano_plot

#Heatmap
library(dplyr)
library(tidyr)
library(ggplot2)
library(pheatmap)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(GO.db)
TPM <- read.csv("salmon.tpm_LDL_Hm.csv", row.names = 1, check.names = FALSE)
samples <- c("Vehicle_1", "Vehicle_2", "Vehicle_3", "LDL_1", "LDL_2", "LDL_3")
TPM <- TPM[, samples, drop = FALSE]
gene_anno <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = rownames(TPM),
  columns = "SYMBOL",
  keytype = "ENSEMBL"
) %>%
  filter(!is.na(SYMBOL), SYMBOL != "") %>%
  distinct(ENSEMBL, .keep_all = TRUE)
TPM$ENSEMBL <- rownames(TPM)
TPM_symbol <- TPM %>%
  inner_join(gene_anno, by = "ENSEMBL")
GO_ids <- c("GO:0008284", "GO:0042060", "GO:0030198","GO:0006869")
GO_genes <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = GO_ids,
  columns = "ENSEMBL",
  keytype = "GOALL"
) %>%
  filter(!is.na(ENSEMBL)) %>%
  distinct(GOALL, ENSEMBL)
colData <- data.frame(condition = factor(c("Vehicle", "Vehicle", "Vehicle","LDL", "LDL", "LDL"), levels = c("Vehicle", "LDL")), row.names = samples)
heatmap_colors <- colorRampPalette(c("blue", "white", "red"))(100)
heatmap_breaks <- seq(-3, 3, length.out = 101)
for (go_id in GO_ids) {
  genes <- GO_genes %>%
    filter(GOALL == go_id) %>%
    pull(ENSEMBL)
  heatmap_df <- TPM_symbol %>%
    filter(ENSEMBL %in% genes) %>%
    select(SYMBOL, all_of(samples)) %>%
    group_by(SYMBOL) %>%
    summarise(across(all_of(samples), mean), .groups = "drop")
  heatmap_mat <- as.matrix(
    heatmap_df[, samples, drop = FALSE])
  rownames(heatmap_mat) <- heatmap_df$SYMBOL
  heatmap_mat <- heatmap_mat[
    apply(heatmap_mat, 1, sd) > 0,
    ,drop = FALSE]
  if (nrow(heatmap_mat) < 2) {
    message("Skipping ", go_id, ": insufficient genes")
    next}
  log_mat <- log2(heatmap_mat + 1)
  z_mat <- t(scale(t(log_mat)))
  z_mat[z_mat > 3] <- 3
  z_mat[z_mat < -3] <- -3
  p <- pheatmap(z_mat, color = heatmap_colors, breaks = heatmap_breaks, annotation_col = colData, cluster_cols = FALSE, cluster_rows = TRUE, show_rownames = TRUE, scale = "none", main = go_id, fontsize_row = 7, border_color = NA)
  pheatmap(
    z_mat,
    color = heatmap_colors,
    breaks = heatmap_breaks,
    annotation_col = colData,
    cluster_cols = FALSE,
    cluster_rows = TRUE,
    show_rownames = TRUE,
    scale = "none",
    main = go_id,
    fontsize_row = 7,
    border_color = NA,
    filename = paste0(
      "Heatmap_",
      gsub(":", "_", go_id),
      ".pdf"
    ),
    width = 7,
    height = max(5, nrow(z_mat) * 0.15)
  )
}