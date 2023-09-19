from hmm_alignment import *
from Bio import SeqIO
import random
import matplotlib.pyplot as plt
import numpy as np
from Bio import pairwise2
import csv

ref = "GAGCTGTACAAGTAAGCGGCCGCCGGACTGCACGACAGTCGACGATGGAGTCGACACGACTCGCGCATACGATGGAGTCGACTACAGTCGCTACGACGATGGAGTCGCGAGCGCTATGAGCGACTATGGAGTCGATACGATACGCGCACGCTATGGAGTCGAGAGCGCGCTCGTCGACTATGGAGTCGCGACTGTACGCACACGCGATGGAGTCGATAGTATGCGTACACGCGATGGAGTCGAGTCGAGACGCTGACGATATGGAGTCGATACGTAGCACGCAGACGATGGGAGCTAGAATTCTAACTAGAGCTCGCTGATCAGCCTCGACTGTGCCTTCTAGTTGC"

data = SeqIO.parse("./data/nhej/SRR11311824/SRR11311824_trim.fastq","fastq")
counter = {}
for record in data:
    t = record.seq._data.decode()
    if t not in counter:
        counter[t] = 1
    else:
        counter[t] += 1
        
counter_filtered = {key: value for key, value in counter.items() if value >= 2 and 
                    key.count("N")<2 and key!=ref and
                    len(key) < round(1.5 * len(ref))}
top_500_sequences = dict(sorted(counter_filtered.items(), key=lambda x: x[1], reverse=True)[:500])
top_sequences = [key for key, value in top_500_sequences.items() for _ in range(value)]

import crispr_hmm

model = crispr_hmm.hmm_model(ref)
model.estimate_param(top_sequences,ncores=5)

filtered_reads = [key for key, value in counter_filtered.items() for _ in range(value)]

NW_result = []
random.seed(10)
for t in filtered_reads:
    align = pairwise2.align.globalms(ref,
                    t,
                    2,
                    -2,
                    -10,
                    -0.1,
                    penalize_end_gaps=1)
    random_alignment = random.sample(align,1)
    NW_result.append(random_alignment[0][:3])
    
with open("./output/NW_carlin_result.csv", 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(["aligned_ref","aligned_read","score"])
    for row in NW_result:
        writer.writerow(row)
        
HMM_result = []

for r in model.viterbi(filtered_reads,ncores=10):
    HMM_result.append(r)

with open("./output/HMM_carlin_result.csv", 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(["aligned_ref","aligned_read","score"])
    for row in HMM_result:
        writer.writerow(row)

