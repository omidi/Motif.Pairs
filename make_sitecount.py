from imports import *

tfbs_directory = "data/closest_bed"
filesnames = [path.join(tfbs_directory, f) for f in
         listdir(tfbs_directory) if re.search('\.bed$', f)]
motif_names = [re.sub('\.closest\.bed$', "", path.basename(f)) for f in filesnames]

motif_files = defaultdict()
for m, f in zip(motif_names, filesnames):
    motif_files[m] = f

with open("motifs_index", "w") as outf:
    for i, m in enumerate(motif_names):
        outf.write('\t'.join([
            m,
            str(i) + '\n',
        ]))

dhs_file = "data/dhs_signal_windows500.chr.sorted.closest.mat"
dhs_sitecounts = defaultdict()
with open(dhs_file) as inf:
    x = 1
    for dhs in csv.reader(inf, delimiter='\t'):
        if not np.mean(map(float, dhs[3:9])):
            continue
        dhs_sitecounts["mm10_%s:%s-%s" % (dhs[0], dhs[1], dhs[2])] = \
            np.zeros(len(motif_names), dtype=np.float32)

for index, motif in enumerate(motif_names):
    with open(motif_files[motif], 'r') as inf:
        for rec in csv.reader(inf, delimiter='\t'):
            dhs_sitecounts[rec[3].split(";")[-1]][index] += float(rec[4])


for dhs_region, sitecounts in dhs_sitecounts.items():
    print '\t'.join([
        dhs_region,
        '\t'.join(map(str, np.round(sitecounts, 5)))
    ])
