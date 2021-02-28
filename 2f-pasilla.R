library("tibble")
library("ggplot2")
library("matrixStats")
library('pasilla')
library('DESeq2')
library("pheatmap")
library("vsn")
library("dplyr")

fn <- system.file("extdata", "pasilla_gene_counts.tsv",
                  package = "pasilla", mustWork = TRUE)
counts <- as.matrix(read.csv(fn, sep = "\t", row.names = "gene_id"))

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
pasRes <- DESeq(pasilla)

rho <- function(x, s)
  ifelse(abs(x) < s, x^2 / 2,  s * abs(x) - s^2 / 2)

df <- tibble(
  x = seq(-7, 7, length.out = 100),
  parabola = x^2 / 2,
  Huber = rho(x, s = 2))

## Figure 8.8
ggplot(reshape2::melt(df, id.vars = "x"),
       aes(x = x, y = value, col = variable)) + geom_line()

pasillaTwoFactor <- pasilla
design(pasillaTwoFactor) <- formula(~type + condition)
pasillaTwoFactor <- DESeq(pasillaTwoFactor)
pasRes2F <- results(pasillaTwoFactor)
resType <- results(pasillaTwoFactor,
                   contrast = c("type", "single", "paired"))
trsf <- function(x) ifelse(is.na(x), 0, (-log10(x))^(1 / 6))
ggplot(as(pasRes2F, "data.frame"), aes(x = pvalue)) +
  geom_histogram(binwidth = 0.01, fill = "Royalblue", boundary = 0)
## Figure 8.9
ggplot(tibble(pOne = pasRes$pvalue,
              pTwo = pasRes2F$pvalue),
       aes(x = trsf(pOne), y = trsf(pTwo))) +
  geom_hex(bins = 75) +
  coord_fixed() +
  xlab("Single factor analysis (condition)") +
  ylab("Two factor analysis (type + condition)") +
  geom_abline(col = "orange")


## Figure 8.10


vsp <- varianceStabilizingTransformation(pasilla)
j <- 1

##Figure 8.11
ggplot(tibble(
  x = assay(pasilla)[, j],
  VST = assay(vsp)[, j],
  log2 = log2(assay(pasilla)[, j])) %>%
         reshape2::melt(id.vars = "x"),
       aes(x = x, y = value, col = variable)) +
  geom_line() +
  xlim(c(0, 600)) +
  ylim(c(0, 9)) +
  xlab("counts") +
  ylab("transformed")

rlp <- rlogTransformation(pasilla)

msd <- function(x)
  meanSdPlot(x, plot = FALSE)$gg +
    ylim(c(0, 1)) +
    theme(legend.position = "none")

## Figure 8.12
gridExtra::grid.arrange(
  msd(log2(counts(pasilla, normalized = TRUE) + 1)) +
    ylab("sd(log2)"),
  msd(assay(vsp)) + ylab("sd(vst)"),
  msd(assay(rlp)) + ylab("sd(rlog)"),
  ncol = 3
)

par(mfrow = c(4, 1), mar = c(2, 2, 1, 1))

myMA <- function(h, v, theta = 0.5) {
  plotMA(pasilla, lfcThreshold = theta, altHypothesis = h,
         ylim = c(-2.5, 2.5))
  abline(h = v * theta, col = "dodgerblue", lwd = 2)
}

## Figure 8.13
gridExtra::grid.arrange(
  myMA("greaterAbs", c(-1, 1)),
  myMA("lessAbs", c(-1, 1)),
  myMA("greater", 1),
  myMA("less", -1),
  nrow = 4
)