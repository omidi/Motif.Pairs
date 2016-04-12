library(gplots)
directory = "/home/somidi/scratch/Tissue.Specificity/Motif.Pairs/Mutual_information/liver.peaks.fg.rhyth.tiss.Liver.bed/"
files = dir(path=directory)

motif_names = gsub('.closest.bed.filt$', '', files)
M = matrix(0, nr=length(files), nc=length(files))
dimnames(M) = list(files, files)

for(file in files) {
  df = read.table(paste(directory, file, sep=''),
                  row.names = 1)
  M[, file] = df[rownames(M),1]
}

dimnames(M) = list(motif_names, motif_names)
# heatmap.2(M, Rowv = False, Colv=False, dendrogram = "none")

max.vec = apply(M, 1, max)
M2 = sweep(M, 1, max.vec, FUN='/')
M3 = M2[names(max.vec[max.vec>0]), names(max.vec[max.vec>0])]
my_palette <- colorRampPalette(c("white", "white", "black"))(n = 1000)
pdf("mutual_information_matrix.pdf", height = 12, width = 12)
heatmap.2(M3, Rowv = FALSE, Colv = FALSE, dendrogram = "none", trace="none", col=my_palette, cexRow=.42, cexCol=.42)
dev.off()

