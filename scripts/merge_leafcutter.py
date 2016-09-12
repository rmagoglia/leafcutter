"""merge_leafcutter.py

Takes output files from leafcutter as well as a bed file of known splice sites
and generates a single file containing relevant info about each splicing 
cluster. Clusters that contain splice sites from multiple genes are discarded.
Clusters containing undocumented splice sites are retained. Run in python3.

"""

import argparse
import gzip
from collections import defaultdict


def get_splice_sites(exon_file):
	"""Build dictionary of splice sites from exon file

	Returns a dictionary of dictionaries:
    splice_dict = {
                'chrom1' : {
                            pos1 : gene,
                            pos2 : gene,
                           },
                'chrom2' : {...},
                ...
                }

	"""

	splice_dict = defaultdict(dict)

	with open(exon_file, 'r') as exons:
		for line in exons:
			chrom, pos1, pos2, strand, gene = line.strip().split()
			splice_dict[chrom][pos1] = gene
			splice_dict[chrom][pos2] = gene

	return splice_dict


def get_cluster_stats(stat_file):
	"""Build dictionary of cluster statistics

	Returns a dictionary:
	stat_dict = {
				'cluster1' : 'cluster_statistics'
				'cluster2' : '...'
				}

	"""

	stat_dict = {}

	with open(stat_file, 'r') as stats:
		stats.readline()
		for line in stats:
			cluster, status, loglr, df, p = line.strip().split('\t')
			cluster = cluster.split(':')[-1]
			stat_dict[cluster] = '\t'.join([status, loglr, df, p])
            
	return stat_dict


def get_genes(chrom, pos1, pos2, splice_dict):
	"""Takes genomic coordinates of a splice site and returns the gene
	or genes corresponding to those splice sites, if any.

	"""

	if chrom in splice_dict:
		if pos1 in splice_dict[chrom]:
			gene1 = splice_dict[chrom][pos1]
		else:
			gene1 = "none"
		if pos2 in splice_dict[chrom]:
			gene2 = splice_dict[chrom][pos2]
		else:
			gene2 = "none"
	else:
		gene1 = "none"
		gene2 = "none"

	return gene1, gene2


def merge_info(out_file, counts_file, splice_dict, stat_dict, is_verbose):
	"""Writes output text file containing the following columns:
	
	chromosome cluster gene status loglr df p-value

	"""

	out = open(out_file, 'w')

	with gzip.open(counts_file, 'rt') as counts_file:
		header = counts_file.readline()
		firstline = counts_file.readline()
		prevChrom, pos1, pos2, prevClus = firstline.strip().split()[0].split(':')
		gene1, gene2 = get_genes(prevChrom, pos1, pos2, splice_dict)
		genes = [gene1, gene2]

		if is_verbose:
			stats = stat_dict[prevClus]
			out.write('\t'.join([prevChrom, pos1, pos2, prevClus, 
							gene1, gene2]) + '\t' + stats + '\n')

		for line in counts_file:
			line =line.strip().split()[0]
			chrom, pos1, pos2, cluster = line.split(':')
			gene1, gene2 = get_genes(chrom, pos1, pos2, splice_dict)
				
			if not is_verbose:
				if cluster == prevClus:
					genes.append(gene1)
					genes.append(gene2)

				else:
					geneSet = set(genes)
					if "none" in geneSet:
						geneSet.remove("none")
					if len(geneSet) == 1:
						gene = geneSet.pop()
						stats = stat_dict[prevClus]
						out.write('\t'.join([prevChrom, prevClus, 
							gene]) + '\t' + stats + '\n')

					genes = [gene1, gene2]
					prevChrom = chrom
					prevClus = cluster

			else:
				stats = stat_dict[cluster]
				out.write('\t'.join([chrom, pos1, pos2, cluster, 
							gene1, gene2]) + '\t' + stats + '\n')
	out.close()


if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument("counts_file", help="Leafcutter counts.gz file.")
	parser.add_argument("stat_file", help=("Cluster statistics file: "
						"leafcutter_ds_cluster_significance.txt"))
	parser.add_argument("exon_file", help=("File of known splice sites; format: "
						"chrom pos1 pos2 strand gene "
						"(no header line)"))
	parser.add_argument("out_file", help="Output cluster file.")
	parser.add_argument("-v", "--verbose", action="store_true", 
						default=False, dest="is_verbose",
						help=("Rather than writing one line per cluster, "
						"write one line per splicing event, with the "
						"corresponding gene(s) for each splice site."))

	options = parser.parse_args()

	SPLICE_DICT = get_splice_sites(options.exon_file)

	STAT_DICT = get_cluster_stats(options.stat_file)

	merge_info(options.out_file, options.counts_file, SPLICE_DICT, STAT_DICT, options.is_verbose)

