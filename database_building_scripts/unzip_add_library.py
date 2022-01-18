#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 20 21:31:04 2021

@author: robynwright
"""

import os
from multiprocessing import Pool
from multiprocessing import freeze_support
import pandas as pd
import argparse

parser = argparse.ArgumentParser(description="A script that will gunzip fasta files, add them to the Kraken2 library and then gzip them again")
parser.add_argument('--genome_folder', dest='genomes', default='fasta_renamed/',
                    help="Folder that contains the fasta files that you'd like to add. These files can be either gzipped or not\
                        but it is expected that they'll be a fasta file\
                        By default this will be fasta_renamed/, the output of the rename_fasta_headers.py script")
parser.add_argument('--database', dest='database', default='Kraken2_GTDB/',
                    help="Name of the Kraken2 database that already contains the taxonomy folder, by default this is Kraken2_GTDB")
parser.add_argument('--already_added', dest='already_added', default=None,
                    help="A text file containing a list of genomes that have already been added to the library and should be skipped \
                        if they are found in the folder containing the fasta files to be added. By default this is None")
parser.add_argument('--processors', dest='proc', default=1,
                    help="Number of processors to use to add files to the library")

args = parser.parse_args()

folder = args.genomes
database = args.database
already_added = args.already_added
processors = args.processors

genomes = os.listdir(folder)

with open(database+'genomes_added.txt', 'w') as f:
    for gen in genomes:
        f.write(gen+'\n')

if already_added != None:
    already_added = pd.read_csv(already_added, header=None, index_col=0, sep='\n')
    already_added = set(list(already_added.index.values))
    
    new_genomes = []
    for gen in genomes:
        if gen not in already_added:
            new_genomes.append(gen)
    genomes = new_genomes

def unzip_add(gen):
    os.system('gunzip '+folder+gen)
    os.system('kraken2-build --add-to-library '+folder+gen.replace('.gz', '')+' --db '+database)
    os.system('gzip '+folder+gen.replace('.gz', ''))
    return

def run_multiprocessing(func, i, n_processors):
    with Pool(processes=n_processors) as pool:
        return pool.map(func, i)

def main():

    print('Starting processing')
    run_multiprocessing(unzip_add, genomes, processors)


if __name__ == "__main__":
    freeze_support()   # required to use multiprocessing
    main()