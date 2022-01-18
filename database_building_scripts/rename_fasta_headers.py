#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 13 11:04:35 2021

@author: robynwright
"""

import pandas as pd
import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import argparse
import gzip
from multiprocessing import Pool
from multiprocessing import freeze_support

parser = argparse.ArgumentParser(description="A script that will rename fasta headers for files based on what is required \
                                 by Kraken2 for database building. It requires a folder full of genomes for renaming \
                                     and a list containing genome accessions and taxid's to use to rename them.\
                                         By default these are taken from the output of download_GTDB_latest.py and are gtdb_genomes_reps/ \
                                             and db_samples.tsv")
parser.add_argument('--genome_folder', dest='genomes', default='gtdb_genomes_reps',
                    help="Folder that contains the fasta files that you'd like to add. These files can be either gzipped or not\
                        but it is expected that they'll be a fasta file\
                        By default this will be gtdb_genomes_reps, the output of the download_GTDB_latest.py script")
parser.add_argument('--genome_list', dest='accession', default='db_samples.tsv',
                    help="Tab-delimited file containing at least two columns - accession in the first column and a second column named 'taxid'\
                        where the accession matches the genome name before an underscore (if present, and without an extension)\
                        By default this is db_samples.tsv, the output from either of the previous scripts")
parser.add_argument('--log_file', dest='log_file', default='logfile_rename_fasta.txt',
                    help="File to write the log to")
parser.add_argument('--processors', dest='proc', default=1,
                    help="Number of processors to use to rename genome files")

                                
args = parser.parse_args()
genomes = args.genomes
accession = args.accession
log_file = args.log_file
n_processors =args.proc

def write_log(message):
    with open(log_file, 'w+') as f:
        f.write(message+'\n')

accession = pd.read_csv(accession, sep='\t', header=0, index_col=0)
acc_dict = {}
for acc in accession.index.values:
  acc_dict[acc] = accession.loc[acc, 'taxid']
print('Made a dictionary with the accession names and taxonomy IDs\n')

genome_folder = genomes+'/'
new_folder = 'fasta_renamed/'
genomes = os.listdir(genome_folder)
count, count2, count3 = 0, 0, 0
logfile = []
os.system('mkdir fasta_renamed/')
new_genomes = []
for genome in genomes:
    tn = genome.replace('.fna', '.tax.fna').replace('.fasta', '.tax.fasta')
    if not os.path.exists(new_folder+tn) and not os.path.exists(new_folder+tn+'.gz'):
        new_genomes.append(genome)
genomes = new_genomes
print('Starting to rename genomes. We have '+str(len(genomes))+' genomes to rename\n')

def rename_genome(genome):
    print(genome)
    if '.gz' not in genome: 
        os.system('gzip '+genome_folder+genome)
        genome = genome+'.gz'
    tn = genome.replace('.fna', '.tax.fna').replace('.fasta', '.tax.fasta')
    if os.path.exists(new_folder+tn):
        write_log('Already had '+new_folder+tn+" so didn't change any fasta headers for "+genome_folder+genome)
        return
    elif os.path.exists(new_folder+tn.replace('.gz', '')):
        write_log('Already had '+new_folder+tn+" so didn't change any fasta headers for "+genome_folder+genome)
        os.system('gzip '+new_folder+tn.replace('.gz', ''))
        return
    try:
        if '_' in genome: 
            acc = genome.split('_')
            acc = acc[0]+'_'+acc[1]
        else: acc = genome
        taxid = acc_dict[acc]
        new_records = []
        with gzip.open(genome_folder+genome, "rt") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                newid = record.id+"|kraken:taxid|"+str(taxid)
                newseq = SeqRecord(record.seq, id=newid, description=record.description)
                new_records.append(newseq)
        with gzip.open(new_folder+tn, "wt") as file:
            SeqIO.write(new_records, file, "fasta")
    except:
        write_log("Couldn't add "+genome)
    return


def run_multiprocessing(func, i, n_processors):
    with Pool(processes=n_processors) as pool:
        return pool.map(func, i)

def main():

    print('Starting processing')
    run_multiprocessing(rename_genome, genomes, int(n_processors))


if __name__ == "__main__":
    freeze_support()   # required to use multiprocessing
    main()



print('Finished! Look at the logfile for any unsuccessful renamings')