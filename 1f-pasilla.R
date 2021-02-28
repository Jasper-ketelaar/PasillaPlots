library("tibble")
library("ggplot2")
library("matrixStats")
library('pasilla')
library('DESeq2')
library("dplyr")

fn <- system.file("extdata", "pasilla_gene_counts.tsv",
                  package = "pasilla", mustWork = TRUE)
counts <- as.matrix(read.csv(fn, sep = "\t", row.names = "gene_id"))
sf <- estimateSizeFactorsForMatrix(counts)

## Figure 8.2
ggplot(tibble(
  `size factor` = sf,
  `sum` = colSums(counts)), aes(x = `size factor`, y = `sum`)) +
  geom_point()


# Gene_id matrix over the size factor matrix
ncounts <- counts / matrix(
  sf, byrow = TRUE, ncol = ncol(counts), nrow = nrow(counts)
)
uncounts <- ncounts[, grep("^untreated", colnames(ncounts)),
                      drop = FALSE]

## Figure 8.3
ggplot(tibble(
  mean = rowMeans(uncounts),
  var = rowVars(uncounts)),
       aes(x = log(mean), y = log(var))) +
  geom_hex() +
  coord_fixed() +
  theme(legend.position = "none") +
  geom_abline(slope = 1:2, color = c("forestgreen", "red"))


annotationFile <- system.file(
  "extdata",
  "pasilla_sample_annotation.csv",
  package = "pasilla", mustWork = TRUE
)

pasillaSampleAnno <- readr::read_csv(annotationFile)
pasillaSampleAnno <- mutate(
  pasillaSampleAnno,
  condition = factor(condition, levels = c("untreated", "treated")),
  type = factor(sub("-.*", "", type), levels = c("single", "paired"))
)

mt <- match(colnames(counts), sub("fb$", "", pasillaSampleAnno$file))
stopifnot(!any(is.na(mt)))
pasilla <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = pasillaSampleAnno[mt,],
  design = ~condition
)
pdseq <- DESeq(pasilla)
pasRes <- results(pdseq)

## Figure 8.4
ggplot(as(pasRes, "data.frame"), aes(x = pvalue)) +
  geom_histogram(binwidth = 0.01, fill = "Royalblue", boundary = 0)

## Figure 8.5
plotMA(pdseq, ylim = c(-2, 2))

pas_rlog <- rlogTransformation(pasilla)

## Figure 8.6
plotPCA(pas_rlog, intgroup = c("condition", "type")) + coord_fixed()


library("pheatmap")
select <- order(rowMeans(assay(psa_rlog)), decreasing = TRUE)[1:30]

## Figure 8.7
pheatmap(assay(psa_rlog)[select,],
         scale = "row",
         annotation_col = as.data.frame(
           colData(psa_rlog)[, c("condition", "type")]))
