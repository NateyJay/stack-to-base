#!/usr/bin/env python3


import sys
from collections import Counter
import argparse
from pprint import pprint
from subprocess import Popen, PIPE
from math import floor

from os.path import isfile, dirname


script_dir = dirname(__file__).rstrip("/.").replace(" ", "\ ")


parser = argparse.ArgumentParser()
parser.add_argument('-r', '--results_file', 
	type=str,
	required=True,
	help='shortstack results file')
parser.add_argument('-g', '--genome_file', 
	type=str,
	required=True,
	help='genome used for shortstack run')
parser.add_argument('-p', '--mirbase_hairpin_file', 
	type=str,
	default = f"{script_dir}/miRBase_files/hairpin.fa",
	required=False,
	help='mirbase hairpin sequences in fasta format')
parser.add_argument('-m', '--mirbase_mature_file', 
	type=str,
	default = f"{script_dir}/miRBase_files/mature.fa",
	required=False,
	help='mirbase mature sequences in fasta format')
parser.add_argument('-o', '--output_name', 
	type=str,
	required=True,
	help='to be filled...')
parser.add_argument('-s', '--species_prefix', 
	type=str,
	required=False,
	default = False,
	help='3 character species prefix if you wish to pre-filter the mirbase hairpin')



args = parser.parse_args()
results_file = args.results_file
genome_file  = args.genome_file
hairpin_file = args.mirbase_hairpin_file
mature_file  = args.mirbase_mature_file
output_name  = args.output_name
species      = args.species_prefix






if species:
	filtered_hp_file = f"{output_name}.filtered_hps.fa"

	if not isfile(filtered_hp_file):

		with open(filtered_hp_file, 'w') as outf:
			with open(hairpin_file, 'r') as f:
				to_out = False
				for line in f:
					if line[0] == '>':
						if f">{species}" in line:
							to_out = True
							print(line.strip(), file=outf)
						else:
							to_out = False

					elif to_out:
						print(line.strip(), file=outf)

	hairpin_file = filtered_hp_file


# print(hairpin_file)


def blat(genome_file, hairpin_file):
	out_file = f"{output_name}.psl"

	if not isfile(out_file):
		call = f"blat {genome_file} {hairpin_file} -out=psl -q=rna {out_file}" 

		p = Popen(call, shell=True, stdout=PIPE, encoding='utf-8')

		out, err = p.communicate()


	return(out_file)



print()
raw_psl = blat(genome_file, hairpin_file)

with open(raw_psl, 'r') as f:
	lines = f.readlines()


header = ['match', 'mis-match', 'rep.match', 'Ns', 'Qgapcount', 'Qgapbases', 'Tgapcount', 'Tgapbases', 'strand', 'Qname', 'Qsize', 'Qstart', 'Qend', 'Tname', 'Tsize', 'Tstart', 'Tend', 'blockcount', 'blockSizes', 'qStarts', 'tStarts']

lines = [l.strip().split() for l in lines[5:]]


lines = [l + [int(l[10]) - int(l[0])] for l in lines]

lines.sort(key=lambda x: float(x[-1]), reverse=False)


# for i, h in enumerate(header):
# 	print(i, h)


hp_bed = f"{output_name}.hps.bed"

with open(hp_bed, 'w') as outf:
	found_queries = set()
	for line in lines:
		qname  = line[9]
		tname  = line[13]
		tstart = int(line[15])
		tend   = int(line[16])
		strand = line[8]

		score = line[-1]

		if qname not in found_queries:
			print(tname, tstart, tend, qname, score, strand, sep='\t', file=outf)

		found_queries.add(qname)



results_bed = f"{output_name}.results.bed"
with open(results_bed, 'w') as outf:
	with open(results_file, 'r') as f:
		for line in f:
			if 'DicerCall' not in line:
				line = line.strip().split("\t")

				locus = line[0]
				name = line[1]

				chrom = locus.split(":")[0]
				start, stop = locus.split(":")[1].split("-")


				print(chrom, start, stop, name, sep='\t', file=outf)


# bedtools intersect -a sly.hps.bed -b sly.results.bed -wb -f .5


def intersect(results_bed, hp_bed):

	call = f"bedtools intersect -a {hp_bed} -b {results_bed} -wb -f .5"

	p = Popen(call.split(), stdout=PIPE, stderr=PIPE, encoding='utf-8')

	out, err = p.communicate()

	out = out.strip().split("\n")

	out = [o.split("\t") for o in out]

	d = {}
	for o in out:
		miRNA = o[3]
		cluster = o[9]

		try:
			d[cluster]
		except KeyError:
			d[cluster] = miRNA

	return(d)


blat_d = intersect(results_bed, hp_bed)


def faidx(locus, file):

	# if not isfile(file+".fai"):
	# 	print("indexing genome...")
	# 	call = f"samtools faidx {file}"

	# 	p = Popen(call.split(), stdout=PIPE, stderr=PIPE, encoding='utf-8')

	# 	out, err = p.communicate()

	# 	print(out)
	# 	print(err)

	call = f"samtools faidx {file} {locus}"

	p = Popen(call.split(), stdout=PIPE, stderr=PIPE, encoding='utf-8')

	out, err = p.communicate()

	if err != "":
		print(call)
		sys.exit(err)

	if out == '':
		sys.exit(f"no sequence found for:\n{call}")

	seq = "".join(out.split("\n")[1:])

	return(seq)


def counter_check(i, total, interval):
	chunk = round(total * (interval/100))
	if i % chunk == 0:
		print(" ", round(i / chunk) * interval)

	elif i == total-1:
		print('  100')

shortstack_fasta = f"{output_name}.results.fa"


with open(results_file, 'r') as f:
	lines = f.readlines()

if isfile(shortstack_fasta):
	print(f"{shortstack_fasta} found, skipping fasta formation...")

else:
	print("extracting fasta sequences from results and genome using faidx...")
	with open(shortstack_fasta, 'w') as outf:
		for i, line in enumerate(lines):

			counter_check(i, len(lines), 5)


			if "DicerCall" not in line:

				line = line.strip().split("\t")

				cluster = line[1]
				# print(cluster)
				locus   = line[0]


				seq = faidx(locus, genome_file)
				# print(seq)
				# sys.exit()

				print(">" + cluster, file=outf)
				print(seq, file=outf)



def best_blast_hit(shortstack_fasta, hp_fasta):
	header = 'qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen qcovs qcovhsp'
	call = f'blastn -query {shortstack_fasta} -subject {hp_fasta} -outfmt 6'# {header}'

	print(call)
	p = Popen(call, shell=True, stdout=PIPE, stderr=PIPE, encoding='utf-8')

	out, err = p.communicate()

	if err != "":
		sys.exit(err)


	blast_file = f"{output_name}.blast.txt"
	with open(blast_file, 'w') as outf:
		print(header.replace(" ", "\t"), file=outf)
		outf.write(out)

	out = out.strip().split("\n")
	out = [o.split("\t") for o in out]

	out.sort(key=lambda x: float(x[11]), reverse=True)



	d = {}
	best_d = {}

	found_clusters = set()
	found_miRNAs = set()

	for o in out:
		cluster, miRNA = o[:2]

		if cluster not in found_clusters:

			d[cluster] = miRNA

		found_clusters.add(cluster)
		found_miRNAs.add(miRNA)

	return(d)

blast_d = best_blast_hit(shortstack_fasta, hairpin_file)



mature_d = {}
with open(mature_file.replace("\ ", " "), 'r') as f:
	for line in f:
		if line[0] == '>':
			header = line.strip().split()[0].lstrip(">")
		else:
			seq = line.strip()

			mature_d[seq] = header


output_file = f"{output_name}.updated_results.txt"
with open(output_file, 'w') as outf:
	with open(results_file, 'r') as f:
		for line in f:
			if "DicerCall" in line:
				line = line.strip().split("\t")

				line.insert(2, "hp_blat_hit")
				line.insert(3, "hp_blast_hit")
				line.insert(4, "mature_hit")

				print("\t".join(line), file=outf)
			else:
				line = line.strip().split('\t')

				try:
					blat_hit = blat_d[line[1]]
				except KeyError:
					blat_hit = '-'	

				try:
					blast_hit = blast_d[line[1]]
				except KeyError:
					blast_hit = '-'

				try:
					mature_hit = mature_d[line[8]]
				except KeyError:
					mature_hit = '-'

				line.insert(2, blat_hit)
				line.insert(3, blast_hit)
				line.insert(4, mature_hit)


				print("\t".join(line), file=outf)


















