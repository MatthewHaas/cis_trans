# 24 November 2019
# making plots for other important barley genes
# WD: /filer-dg/agruppen/DG/haas/rna_cold

# Read in data
load("170508_expression_data.Rdata")

# Read in gene info
fread("/filer-dg/agruppen/DG/mascher/gene_annotation_160429/160504_gene_info_with_header.tsv") -> gene_info

# Define accessions
accessions=c("Morex", "Barke", "Igri", "BCC131", "HOR1969", "FT11", "FT67", "FT581", "FT279")
accessions = as.data.table(accessions)

pdf("out.pdf", width=12, height=5)
layout(matrix(c(1,2), nrow=1), widths=c(6,6))
par(oma=c(0,0,6,0))

# VRS1
expr_mat[gene == "HORVU4Hr1G007040"] -> x
accessions[, sample_index := as.integer(NA)]
accessions[, sample_index := c(1:9)]
setnames(accessions, c("accession", "sample_index"))
accessions[x, on="accession"]->x

par(mar=c(4,4,3,1))
t="control"
x[, plot(ylim=c(-1, 20), xlab="sample", ylab="log2(transcript per million reads mapped)", xlim=c(1,9), main=t, type='n', yaxt='n', 1, las=1, xaxt='n')]
abline(col="gray", lwd=2, h=-1)
axis(2, las=1)
legend("topleft", horiz=T, pch=16, col=1:2, c("Parent", "Hybrid"))
x[generation == "parent" & sample_index >= 1 & treatment == t, points(sample_index, col=rgb(0,0,0, 0.5), pch=19, log2(0.5+count))]
x[generation == "hybrid" & sample_index >= 1 & treatment == t, points(sample_index, col=rgb(1,0,0, 0.5), pch=19, log2(0.5+count))]
axis(1, 1:9, c("Morex", "Barke", "Igri", "BCC131", "HOR1969", "FT11", "FT67", "FT581", "FT279"), cex=0.5)
#axis(1, 1:9, accessions[sample_index %in% 1:9, accession], cex.lab=0.5)

par(mar=c(4,4,3,1))
t="cold"
x[, plot(ylim=c(-1, 20), xlab="sample", ylab="", xlim=c(1,9), main=t, type='n', yaxt='n', 1, las=1, xaxt='n')]
abline(col="gray", lwd=2, h=-1)
axis(2, las=1)
x[generation == "parent" & sample_index >= 1 & treatment == t, points(sample_index, col=rgb(0,0,0, 0.5), pch=19, log2(0.5+count))]
x[generation == "hybrid" & sample_index >= 1 & treatment == t, points(sample_index, col=rgb(1,0,0, 0.5), pch=19, log2(0.5+count))]
axis(1, 1:9, c("Morex", "Barke", "Igri", "BCC131", "HOR1969", "FT11", "FT67", "FT581", "FT279"), cex=0.5)
#axis(1, 1:9, accessions[sample_index %in% 1:9, accession], cex.lab=0.5)
title(outer=T, paste0("Homeobox-leucine zipper protein family (VRS1)", "\n", "HORVU2Hr1G092290"))

# INT-C
expr_mat[gene == "HORVU4Hr1G007040"] -> x
accessions[, sample_index := as.integer(NA)]
accessions[, sample_index := c(1:9)]
setnames(accessions, c("accession", "sample_index"))
accessions[x, on="accession"]->x

par(mar=c(4,4,3,1))
t="control"
x[, plot(ylim=c(-1, 20), xlab="sample", ylab="log2(transcript per million reads mapped)", xlim=c(1,9), main=t, type='n', yaxt='n', 1, las=1, xaxt='n')]
abline(col="gray", lwd=2, h=-1)
axis(2, las=1)
legend("topleft", horiz=T, pch=16, col=1:2, c("Parent", "Hybrid"))
x[generation == "parent" & sample_index >= 1 & treatment == t, points(sample_index, col=rgb(0,0,0, 0.5), pch=19, log2(0.5+count))]
x[generation == "hybrid" & sample_index >= 1 & treatment == t, points(sample_index, col=rgb(1,0,0, 0.5), pch=19, log2(0.5+count))]
axis(1, 1:9, c("Morex", "Barke", "Igri", "BCC131", "HOR1969", "FT11", "FT67", "FT581", "FT279"), cex=0.5)
#axis(1, 1:9, accessions[sample_index %in% 1:9, accession], cex.lab=0.5)

par(mar=c(4,4,3,1))
t="cold"
x[, plot(ylim=c(-1, 20), xlab="sample", ylab="", xlim=c(1,9), main=t, type='n', yaxt='n', 1, las=1, xaxt='n')]
abline(col="gray", lwd=2, h=-1)
axis(2, las=1)
x[generation == "parent" & sample_index >= 1 & treatment == t, points(sample_index, col=rgb(0,0,0, 0.5), pch=19, log2(0.5+count))]
x[generation == "hybrid" & sample_index >= 1 & treatment == t, points(sample_index, col=rgb(1,0,0, 0.5), pch=19, log2(0.5+count))]
axis(1, 1:9, c("Morex", "Barke", "Igri", "BCC131", "HOR1969", "FT11", "FT67", "FT581", "FT279"), cex=0.5)
#axis(1, 1:9, accessions[sample_index %in% 1:9, accession], cex.lab=0.5)
title(outer=T, paste0("Transcription factor TEOSINTE BRANCHED 1 (INT-C)", "\n", "HORVU4Hr1G007040"))

# VRN3
pdf("191124_VRN3_candidate_expression_profile.pdf", width=20, height=5)
layout(matrix(c(1,2), nrow=1), widths=c(10,10))
par(oma=c(0,0,6,0))

expr_mat[gene == "HORVU7Hr1G024610"]->x
accessions[, sample_index := as.integer(NA)]
accessions[, sample_index := c(1:9)]
setnames(accessions, c("accession", "sample_index"))
accessions[x, on="accession"]->x


par(mar=c(4,4,3,1))
t="control"
x[, plot(ylim=c(-1, 20), xlab="sample", ylab="log2(transcript per million reads mapped)", xlim=c(1,9), main=t, type='n', yaxt='n', 1, las=1, xaxt='n')]
abline(col="gray", lwd=2, h=-1)
axis(2, las=1)
legend("topleft", horiz=T, pch=16, col=1:2, c("Parent", "Hybrid"))
x[generation == "parent" & sample_index >= 1 & treatment == t, points(sample_index, col=rgb(0,0,0, 0.5), pch=19, log2(0.5+count))]
x[generation == "hybrid" & sample_index >= 1 & treatment == t, points(sample_index, col=rgb(1,0,0, 0.5), pch=19, log2(0.5+count))]
axis(1, 1:9, c("Morex", "Barke", "Igri", "BCC131", "HOR1969", "FT11", "FT67", "FT581", "FT279"), cex=0.5)
#axis(1, 1:9, accessions[sample_index %in% 1:9, accession], cex.lab=0.5)

par(mar=c(4,4,3,1))
t="cold"
x[, plot(ylim=c(-1, 20), xlab="sample", ylab="", xlim=c(1,9), main=t, type='n', yaxt='n', 1, las=1, xaxt='n')]
abline(col="gray", lwd=2, h=-1)
axis(2, las=1)
x[generation == "parent" & sample_index >= 1 & treatment == t, points(sample_index, col=rgb(0,0,0, 0.5), pch=19, log2(0.5+count))]
x[generation == "hybrid" & sample_index >= 1 & treatment == t, points(sample_index, col=rgb(1,0,0, 0.5), pch=19, log2(0.5+count))]
axis(1, 1:9, c("Morex", "Barke", "Igri", "BCC131", "HOR1969", "FT11", "FT67", "FT581", "FT279"), cex=0.5)
#axis(1, 1:9, accessions[sample_index %in% 1:9, accession], cex.lab=0.5)
title(outer=T, paste0("FLOWERING LOCUS T 1 (VRN3)", "\n", "HORVU7Hr1G024610"))
dev.off()
