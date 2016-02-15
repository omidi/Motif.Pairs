from imports import *

sitecount_file = argv[1]
motif_index_file = argv[2]

motif_index = defaultdict()
with open(motif_index_file) as inf:
    for line in inf:
        motif_index[line.split()[0]] = int(line.split()[1])
motif_names = motif_index.keys()
num_motifs = len(motif_index.keys())

mat = np.matrix(np.zeros(num_motifs**2)).reshape(num_motifs,num_motifs)

with open(sitecount_file) as inf:
    for rec in csv.reader(inf, delimiter='\t'):
        counts = np.array(map(float, rec[1:]))
        for pairs in combinations(np.arange(num_motifs), 2):
            mat[pairs] = mat[pairs[::-1]] = counts[pairs[0]] + counts[pairs[1]]