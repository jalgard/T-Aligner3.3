import sys
from collections import defaultdict

regions_data = defaultdict(lambda : [])
with open(sys.argv[1], 'r') as good_regions_file:
    for line in good_regions_file:
        toks = line.rstrip().split('\t')
        if len(toks) > 2:
            for region in toks[2:]:
                c1, c2 = region.split(' ')[0].split(':')
                for p in range(int(c1) - 10, int(c2) + 10):
                    regions_data[toks[0]].append(p)

# name the regions
regions = defaultdict(lambda : defaultdict(lambda : 'Maxi'))
midofcircle = 600

for mini in regions_data:
    for pos in regions_data[mini]:
        if pos < midofcircle:
            regions[mini][pos] = mini + '-L'
        else:
            regions[mini][pos] = mini + '-R'

# function to cast coordinate to t-less space
def MapOnRef(coord, ref_name, ref_bank):
    seq = ref_bank[ref_name]
    pos = 0
    cur = 0
    while cur < coord:
        if seq[cur] != 'T':
            pos += 1
        cur += 1
    return pos

# load fasta files with mRNA sequences
ref_bank = {}
last_entry = ''

with open(sys.argv[2], 'r') as ifile:
    for line in ifile:
        if line[0] == '>':
            last_entry = line[1:].rstrip()
            ref_bank[last_entry] = ''
        else:
            ref_bank[last_entry] += line.rstrip().upper()


gene_name_given = True
gene_name = sys.argv[3]
if gene_name == '_':
    gene_name_given = False

grnas_seen = []
mrna_coverage_raw = defaultdict(lambda : [])

output_prefix = sys.argv[4]
main_ofile = open(output_prefix + '_by_gene_table.txt', 'w')

for line in sys.stdin:
    toks = line.rstrip().split('\t')
    target_name  = toks[0]
    source_name  = toks[2]
    target_begin = int(toks[1])
    source_begin = int(toks[3])
    target_seq   = toks[4]
    source_seq   = toks[5]
    align_gu     = int(toks[6])
    align_mm     = int(toks[7])
    if source_begin in regions_data[source_name] or 'H10_complete_maxicircle' in source_name :
        # simple output
        #sys.stdout.writelines('{}\n'.format(line))
        # region-centered output
        if gene_name_given is False:
            gene_name = target_name.upper()
        begin = MapOnRef(target_begin, target_name, ref_bank)
        end   = MapOnRef(target_begin + len(target_seq), target_name, ref_bank)
        # edited region should contain more then 10 A/G/C
        if end - begin < 7:
            continue
        grna = gene_name + ':' + regions[source_name][source_begin+5] + ':' + str(begin) + ':' + str(end)
        if grna not in grnas_seen:
            main_ofile.writelines('{}\t{}\t{}\t{}\t{}\t{}\n'.format(gene_name, regions[source_name][source_begin+1], begin, end, target_seq, source_seq))
            mrna_coverage_raw[target_name].append([begin, end, gene_name, regions[source_name][source_begin+1], source_seq])
            grnas_seen.append(grna)

main_ofile.close()
mrna_cov_ofile = open(output_prefix + '_mrna_cov.txt', 'w')

mrnas_coverage_data = []
for mrna in mrna_coverage_raw:
    mrna_coverage_raw[mrna].sort(key = lambda x : x[0])
    mrnas_coverage_data.append([mrna, mrna_coverage_raw[mrna]])

mrnas_coverage_data.sort(key = lambda x : len(x[1]), reverse=True)

for mrna in mrnas_coverage_data:
    mrna_cov_ofile.writelines('{}'.format(mrna[0]))
    for grna in mrna[1]:
        mrna_cov_ofile.writelines('\t{},{},{},{},{}'.format(grna[0], grna[1], grna[2], grna[3], grna[4]))
    mrna_cov_ofile.writelines('\n')


mrna_cov_ofile.close()

#cat ./results_mm3_gu13_ml20/ND8_ORFs.txt | python filter_grna_finder_hits.py filtered_regions.txt ND8_alt_polyAmix_merged_filtered_sl10_besthits1_rev_199_217_54_ORFs_mrna.fasta ND8 | sort -nk 3,4 | uniq |awk '{if(length($5 > 25)) print $0}' | cut -f2 | sort | uniq | wc -l
