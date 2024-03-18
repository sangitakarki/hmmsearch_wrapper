#!/usr/bin/env python
import sys, os, re, shlex, subprocess, numpy, itertools, argparse, time
from collections import defaultdict
import operator

# This is a generic wrapper for hmmsearch. If provided a file of proteins (.faa) the script will run an hmmsearch against a specified .hmm database and output 1) a table of the proteins with hits, and their best-hit HMM, 2) a fasta file of all the proteins that had hits,

inputfile = sys.argv[1]
hmmdb = sys.argv[2]
outputfile = open(sys.argv[3], "w")
CPUS = 4


# run HMMER3
def run_hmmer(input_file, cpus):
	output_file = re.sub(".faa", ".hmmout", input_file)
	cmd = "hmmsearch --cut_nc --cpu " + str(CPUS) + " --domtblout " + output_file + " " + hmmdb + " " + input_file
	print(cmd)
	cmd2 = shlex.split(cmd)
	#subprocess.call(cmd2, stdout=open("out.txt", "w"), stderr=open("err.txt", "w"))
	#os.remove("out.txt")
	return output_file


# define function for parsing HMMER3 output
def parse_hmmout(hmmout):
	input = open(hmmout, "r")
	final_dict = defaultdict(int)
	hit_dict = defaultdict(list)
	bit_dict = defaultdict(list)
	outputfile.write("Protein\tnum_domains\thits\tbits\tcoords\n")

	done_dict = {}

	for line in input.readlines():
		if line.startswith("#"):
			pass
		else:

			newline = re.sub('\s+', '\t', line)
			list1 = newline.split('\t')
			ids = list1[0]
			hit = list1[3]

			hit_ids = ids+"|"+hit

			print(newline)
			start = int(list1[17])
			end = int(list1[18])
			score = float(list1[7])
			domain_evalue = float(list1[11])

			if hit_ids in done_dict:
				pass
			else:			
				done_dict[hit_ids] = hit_ids
				hit_dict[ids].append(hit)
				bit_dict[ids].append(score)

		# print(ids, hit, start, end, coords, score)
	return hit_dict, bit_dict

hmmout = run_hmmer(inputfile, 6)
outputfile = open(re.sub(".hmmout", ".parsed", hmmout), "w")
hit_dict, bit_dict = parse_hmmout(hmmout)

outputfile.write("QUERY\tHIT\tBITSCORE\n")
for i in hit_dict:
	hits = ";".join(hit_dict[i])
	bits = ";".join([str(n) for n in bit_dict[i]])
	outputfile.write(str(i) +"\t"+ str(hits) +"\t"+ str(bits) +"\n")


