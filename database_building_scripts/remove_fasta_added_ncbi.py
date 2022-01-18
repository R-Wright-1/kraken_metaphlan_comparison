#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 19 15:02:31 2021

@author: robynwright
"""

import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import gzip
from multiprocessing import Pool
import pickle
from multiprocessing import freeze_support
import pandas as pd


genome_folder = 'ncbi_genomes_all/'
new_folder = 'RefSeqV205_Complete/fasta_renamed/'
genomes = os.listdir(genome_folder)
all_genomes = []
added = []

for gen in genomes:
    if not os.path.exists(new_folder+gen.replace('.fna', '.tax.fna')):
        all_genomes.append(gen)
    else:
        added.append(gen)

print(len(added))
print(len(all_genomes))

for gen in added:
    os.system('rm '+genome_folder+gen)