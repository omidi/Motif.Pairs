library(gplots)
directory = c("/home/somidi/scratch/Tissue.Specificity/Motif.Pairs/Mutual_information/nonliver.peaks.bg.rhyth.tiss.Liver.bed/", 
              "/home/somidi/scratch/Tissue.Specificity/Motif.Pairs/Mutual_information/liver.peaks.bg.flat.bed/",
              "/home/somidi/scratch/Tissue.Specificity/Motif.Pairs/Mutual_information/liver.peaks.fg.rhyth.tiss.Liver.bed/")

Mat = list()
index = 1
for(dir in directory) {
  files = dir(path=dir)
  motif_names = gsub('.closest.bed.filt$', '', files)
  # M = matrix(0, nr=length(files), nc=length(files))
  # dimnames(M) = list(files, files)
  # for(file in files) {
  #   df = read.table(paste(dir, file, sep=''),
  #                   row.names = 1)
  #   M[, file] = df[rownames(M),1]
  # }
  
  Mz = matrix(0, nr=length(files), nc=length(files))
  dimnames(Mz) = list(files, files)
  for(file in files) {
    df = read.table(paste(dir, file, sep=''),
                    row.names = 1)
    Mz[, file] = df[rownames(Mz),2]
  }
  dimnames(Mz) = list(motif_names, motif_names)
  filter = rownames(Mz)[! is.na(diag(Mz))]
  Mz = Mz[filter, filter]
  Mat[[index]] = Mz
  index = index + 1
}

len = c(1824, 1243, 506)

# motif_names = gsub('.closest.bed.filt$', '', files)
# M = matrix(0, nr=length(files), nc=length(files))
# dimnames(M) = list(files, files)
# 
# for(file in files) {
#   df = read.table(paste(directory, file, sep=''),
#                   row.names = 1)
#   M[, file] = df[rownames(M),1]
# }
# 
# dimnames(M) = list(motif_names, motif_names)
# # heatmap.2(M, Rowv = False, Colv=False, dendrogram = "none")
# 
# max.vec = apply(M, 1, max)
# # M2 = sweep(M, 1, max.vec, FUN='/')
# # M2 = M2[rownames(M), colnames(M)]
# # M3 = M2[names(max.vec[max.vec>0]), names(max.vec[max.vec>0])]
# my_palette <- colorRampPalette(c("white", "white", "black"))(n = 1000)
# # # pdf("mutual_information_matrix.pdf", height = 12, width = 12)
# # heatmap.2(M, Rowv = FALSE, Colv = FALSE, dendrogram = "none", trace="none", col=my_palette, cexRow=.42, cexCol=.42)
# # dev.off()
# 
# motif_names = gsub('.closest.bed.filt$', '', files)
# Mz = matrix(0, nr=length(files), nc=length(files))
# dimnames(Mz) = list(files, files)
# 
# for(file in files) {
#   df = read.table(paste(directory, file, sep=''),
#                   row.names = 1)
#   Mz[, file] = df[rownames(Mz),2]
# }
# 
# dimnames(Mz) = list(motif_names, motif_names)
# 
# motif_names = gsub('.closest.bed.filt$', '', files)
# Mf = matrix(0, nr=length(files), nc=length(files))
# dimnames(Mf) = list(files, files)
# 
# for(file in files) {
#   df = read.table(paste(directory, file, sep=''),
#                   row.names = 1)
#   Mf[, file] = df[rownames(Mz),3] / df[file, 3]
# }
# 
# dimnames(Mf) = list(motif_names, motif_names)
# # 
# 
# filter = rownames(Mz)[! is.na(diag(Mz))]
# Mz = Mz[filter, filter]
# Mf = Mf[filter, filter]