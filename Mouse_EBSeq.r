#Salmon on Mac terminal
conda install bioconda::fastp
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

conda install bioconda::salmon
wget ftp://ftp.ensembl.org/pub/release-110/fasta/mus_musculus/cdna/Mus_musculus.GRCm39.cdna.all.fa.gz
salmon index -t Mus_musculus.GRCm39.cdna.all.fa.gz -i Mus_musculus.GRCm39_salmon_index
INDEX_mouse="Mus_musculus.GRCm39_salmon_index"
samples=(
Ctrl_1
Ctrl_2
Ctrl_3
WT_1
WT_2
WT_3
KO_1
KO_2
KO_3
)

for sample in "${samples[@]}"
do
    echo "Running ${sample}..."

    salmon quant \
        -i "$INDEX_mouse" \
        -l IU \
        -1 "out_${sample}_1.fastq.gz" \
        -2 "out_${sample}_2.fastq.gz" \
        --validateMappings
        -o "${sample}_salmon_quant"
done

#Tximport on R
library(tximport)
library(readr)
library(biomaRt)
salmon_dirs <- list.files(".", pattern = "_salmon_quant$", full.names = TRUE)
samples <- sub("_salmon_quant$", "", basename(salmon_dirs))
salmon.files <- file.path(salmon_dirs, "quant.sf")
names(salmon.files) <- samples
mart <- useEnsembl(biomart = "genes",dataset = "mmusculus_gene_ensembl", version = 110)
tx2gene <- getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id"),mart = mart)
colnames(tx2gene) <- c("TXNAME", "GENEID")
txi <- tximport(salmon.files,type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)
count_table <- as.data.frame(txi$counts)
write.csv(count_table, file = "salmon.count_LDL_Ms.csv", row.names = TRUE)
tpm_table <- as.data.frame(txi$abundance)
write.csv(tpm_table, file = "salmon.tpm_LDL_Ms.csv", row.names = TRUE)

#EBSeq on R
library(EBSeq)
library(biomaRt)
library(tidyr)
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(clusterProfiler)
library(org.Mm.eg.db)
library(pheatmap)
Mm_counts <- read.csv("salmon.count_LDL_Ms.csv", row.names = 1, check.names = FALSE)
samples <- c("Ctrl1","Ctrl2","Ctrl3","WT1","WT2","WT3","KO1","KO2","KO3")
EBseq_input <- as.matrix(Mm_counts[, samples])
keep <- rowMeans(EBseq_input) >= 1
EBseq_input <- EBseq_input[keep, , drop = FALSE]
EBseq_input <- round(EBseq_input)
Conditions=c("Ctrl","Ctrl","Ctrl","WT","WT","WT","KO","KO","KO")  
PosParti=GetPatterns(Conditions)
PosParti
MultiSize = MedianNorm(EBseq_input)
MultiOut = EBMultiTest(EBseq_input, NgVector=NULL,Conditions= Conditions, AllParti=PosParti, sizeFactors=MultiSize, maxround=5)
MultiPP = GetMultiPP(MultiOut)
MultiFC = GetMultiFC(MultiOut)
table(MultiPP$MAP)
head(MultiFC$PostFCMat)
FC<-as.data.frame(MultiFC$PostFCMat)
FC$id <- rownames(FC)
PP<-as.data.frame(MultiPP$PP)
PP$id <- rownames(PP)
MAP <- data.frame(id = names(MultiPP$MAP), MAP = as.character(MultiPP$MAP))
FC_PP_MAP <- FC %>%
  left_join(PP, by = "id") %>%
  left_join(MAP, by = "id")
pattern_cols <- colnames(PP)[colnames(PP) != "id"]
FC_PP_MAP$PP <- apply(FC_PP_MAP[, pattern_cols, drop = FALSE], 1, max)
FC_PP_MAP <- FC_PP_MAP %>%
  rename(ensembl_gene_id = id)
MultiPP$MAP == "pattern2" | MultiPP$MAP == "pattern5"
Mm_TPM <- read.csv("salmon.tpm_LDL_Ms.csv", row.names = 1, check.names = FALSE)
Mm_TPM$ensembl_gene_id <- rownames(Mm_TPM)
FC_PP_MAP_TPM <- left_join(FC_PP_MAP, Mm_TPM, by = "ensembl_gene_id")
gene_symbol <- AnnotationDbi::select(org.Mm.eg.db, keys = unique(FC_PP_MAP_TPM$ensembl_gene_id),columns = "SYMBOL",keytype = "ENSEMBL") %>%
  distinct(ENSEMBL, .keep_all = TRUE)
FC_PP_MAP_TPM <- FC_PP_MAP_TPM %>%
    left_join(
    gene_symbol,
    by = c("ensembl_gene_id" = "ENSEMBL")) %>%
    rename(external_gene_name = SYMBOL)

#plot
FC_PP_MAP_TPM_up <- FC_PP_MAP_TPM %>%
  filter(DE == "LDL_high") %>%
  slice_max(order_by = PP, n = 10, with_ties = FALSE)
FC_PP_MAP_TPM_down <- FC_PP_MAP_TPM %>%
  filter(DE == "LDL_low") %>%
  slice_max(order_by = PP, n = 10, with_ties = FALSE)
FC_PP_MAP_TPM_top10 <- c(FC_PP_MAP_TPM_up$external_gene_name, FC_PP_MAP_TPM_down$external_gene_name)
FC_PP_MAP_TPM <- FC_PP_MAP_TPM %>%
  mutate(
    label = if_else(
      !is.na(external_gene_name) &
        external_gene_name %in% FC_PP_MAP_TPM_top10,
      "T", ""
    ),
    DE = factor(DE, levels = c("NC", "LDL_high", "LDL_low")),
    FDR = 1 - PP
  ) %>%
  arrange(DE != "NC")
volcano_EBsea <- ggplot(FC_PP_MAP_TPM,aes(x = log2(1/CtrkOverWT),y = log2(1/KOOverWT),size = FDR, color = DE)) +
  geom_point(alpha = 0.5) +
  labs(x = "log2 (FC_WT/Ctrl)",　y = "log2 (FC_WT/KO)") +
  scale_color_manual(values = c(
      "NC" = "grey", "LDL_high" = "red", "LDL_low" = "blue")) +
  scale_size(range = c(3, 0.3), limits = c(0.6, 0), breaks = seq(0, 0.6, by = 0.1)) +
  geom_hline(yintercept = log2(c(1.5, 1/1.5)), linewidth = 0.2, color = "grey") +
  geom_vline(xintercept = log2(c(1.5, 1/1.5)), linewidth = 0.2, color = "grey") +
  theme_classic() +
  geom_text_repel(data = subset(FC_PP_MAP_TPM, label == "T"), aes(label = external_gene_name), max.overlaps = Inf)
plot(volcano_EBsea)
FC_PP_MAP_TPM %>% count(DE)

#mouse_GO
WT_DE<- filter(FC_PP_MAP_TPM,DE == "WT_hight"|DE == "WT_low")
MmGeneList<-WT_DE$ensembl_gene_id

#GO_all
GO_all <- enrichGO(MmGeneList,
                   OrgDb=org.Mm.eg.db,
                   keyType = "ENSEMBL",
                   ont="all", 　　　　　　　　
                   pvalueCutoff=1,
                   pAdjustMethod="BH",
                   qvalueCutoff=1, 
                   readable=TRUE) 
GO_all_txt<-data.frame(GO_all)

GO_all_txt<-GO_all_txt%>%mutate(neglog10_pvalue=-log10(pvalue))
GO_all_txt$Description <- reorder(GO_all_txt$Description, GO_all_txt$neglog10_pvalue)
GO_BP<-GO_all_txt %>% dplyr::filter(ONTOLOGY =="BP")
GO_BP_top10p <- GO_BP[order(GO_BP$pvalue)[1:10], ]
GO_CC<-GO_all_txt %>% dplyr::filter(ONTOLOGY =="CC")
GO_CC_top10p <- GO_CC[order(GO_CC$pvalue)[1:10], ]
GO_MF<-GO_all_txt %>% dplyr::filter(ONTOLOGY =="MF")
GO_MF_top10p <-GO_MF[order(GO_MF$pvalue)[1:10], ]
GO_combined_Mm<-rbind(GO_BP_top10p,GO_CC_top10p,GO_MF_top10p)
GO_combined_Mm$ONTOLOGY <- factor(GO_combined_Mm$ONTOLOGY , levels = c("BP", "MF", "CC"))
GO_combined_Mm$Description <- reorder(GO_combined_Mm$Description, GO_combined_Mm$neglog10_pvalue)
range(GO_combined_Mm$neglog10_pvalue)
range(GO_combined_Mm$Count)
range(GO_combined_Mm$RichFactor)
GO_Mm<-ggplot(GO_combined_Mm, aes(x = neglog10_pvalue, y = Description)) +
  geom_point(aes(size = Count, color = RichFactor)) +
  scale_color_continuous(low = "blue", high = "red", limits = c(0,0.4),breaks = seq(0,0.4,by=0.1)) + #richfactor
  scale_size(range = c(0, 4),limits = c(0,10), breaks = seq(2,10,2))+ #count
  scale_x_continuous(limits = c(0, 6),breaks = seq(0,6,1))+ #pvalue
  facet_grid(ONTOLOGY ~ ., scales = "free_y") +
  theme_bw()
GO_Mm

#heatmap
WT_DE_TPM<-WT_DE[,c(1,14:23)]
rownames(WT_DE_TPM)<-WT_DE_TPM$external_gene_name
colData<- data.frame(row.names=colnames(WT_DE_TPM[,2:10]),condition=c("Ctrl","Ctrl","Ctrl","WT","WT","WT","KO","KO","KO"))
colData$condition <- factor(colData$condition, levels = c("Ctrl","WT","KO"))
GO0008284<-read.table("GO_term_summary_GO0008284_20260702_205206.txt",header =TRUE)
GO0008284 <- GO0008284 %>% distinct(Symbol, .keep_all = TRUE)
WT_DE_TPM_GO0008284 <-WT_DE_TPM[WT_DE_TPM$external_gene_name %in% GO0008284$Symbol,]
heatap_GO0008284<-pheatmap(WT_DE_TPM_GO0008284[2:10], show_rownames = T, scale = "row", cluster_cols = F,colorRampPalette(c("blue", "white","red"))(100),breaks = seq(-3, 3, length.out = 100),annotation_col = colData)

GO1903053<-read.table("GO_term_summary_ GO1903053_20260723_214018.txt",header =TRUE)
GO1903053 <- GO1903053 %>% distinct(Symbol, .keep_all = TRUE)
WT_DE_TPM_GO1903053 <-WT_DE_TPM[WT_DE_TPM$external_gene_name %in% GO1903053$Symbol,]
heatap_GO1903053<-pheatmap(WT_DE_TPM_GO1903053[2:10], show_rownames = T, scale = "row", cluster_cols = F,colorRampPalette(c("blue", "white","red"))(100),breaks = seq(-3, 3, length.out = 100),,annotation_col = colData)

GO0006869<-read.table("GO_term_summary_GO0006869_20260723_214944.txt",header =TRUE)
GO0006869<- GO0006869 %>% distinct(Symbol, .keep_all = TRUE)
WT_DE_TPM_GO0006869 <-WT_DE_TPM[WT_DE_TPM$external_gene_name %in% GO0006869$Symbol,]
heatap_GO0006869<-pheatmap(WT_DE_TPM_GO0006869 [2:10], show_rownames = T, scale = "row", cluster_cols = F,colorRampPalette(c("blue", "white","red"))(100),breaks = seq(-3, 3, length.out = 100),,annotation_col = colData)

GO0042060<-read.table("GO_term_summary_GO0042060_20260723_215342.txt",header =TRUE)
GO0042060<- GO0042060 %>% distinct(Symbol, .keep_all = TRUE)
WT_DE_TPM_GO0042060 <-WT_DE_TPM[WT_DE_TPM$external_gene_name %in% GO0042060$Symbol,]
heatap_GO0006869<-pheatmap(WT_DE_TPM_GO0042060 [2:10], show_rownames = T, scale = "row", cluster_cols = F,colorRampPalette(c("blue", "white","red"))(100),breaks = seq(-3, 3, length.out = 100),,annotation_col = colData)

ALL_TPM<-FC_PP_MAP_TPM[,c(1,14:23)]
GO00224024<-read.table("GO_term_summary_GO00224024_20260702_232509.txt",header =TRUE)
GO00224024<- GO00224024 %>% distinct(Symbol, .keep_all = TRUE)
GO00224024_Ki67_Top2a_CCND_CCNE_CDK<- GO00224024 %>%filter(str_starts(Symbol, "^(Mki67|Top2a|Ccnd|Ccne|Cdk)"))
rownames(GO00224024_Ki67_Top2a_CCND_CCNE_CDK)<-GO00224024_Ki67_Top2a_CCND_CCNE_CDK$external_gene_name
WT_DE_TPM_GO00224024_Ki67_Top2a_CCND_CCNE_CDK <-ALL_TPM[ALL_TPM$external_gene_name %in% GO00224024_Ki67_Top2a_CCND_CCNE_CDK$Symbol,]
rownames(WT_DE_TPM_GO00224024_Ki67_Top2a_CCND_CCNE_CDK)<-WT_DE_TPM_GO00224024_Ki67_Top2a_CCND_CCNE_CDK$external_gene_name
heatap_GO00224024 <-pheatmap(WT_DE_TPM_GO00224024_Ki67_Top2a_CCND_CCNE_CDK[2:10], show_rownames = T, scale = "row", cluster_cols = F,colorRampPalette(c("blue", "white","red"))(100),breaks = seq(-3, 3, length.out = 100),annotation_col = colData,cutree_rows = 6)