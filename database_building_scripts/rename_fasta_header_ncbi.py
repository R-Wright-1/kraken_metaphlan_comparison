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

domains = ['bacteria','archaea','fungi','protozoa','viral','vertebrate_mammalian','invertebrate','vertebrate_other','plant']
accession = []

for dom in domains:
    dom = pd.read_csv('assembly_lists/'+dom+'_assembly_summary.txt', index_col=0, header=1, sep='\t')
    dom = dom.loc[:, ['taxid', 'ftp_path']]
    accession.append(dom)
    
accession = pd.concat(accession)

acc_dict = {}
for acc in accession.index.values:
    acc_dict[accession.loc[acc, 'ftp_path'].split('/')[-1]] = accession.loc[acc, 'taxid']

with open('ncbi_taxid.dict', 'wb') as f:
    pickle.dump(acc_dict, f)

# with open('ncbi_taxid.dict', 'rb') as f:
#         acc_dict = pickle.load(f)

print(len(acc_dict))

genome_folder = 'ncbi_genomes_all/'
new_folder = 'RefSeqV205_Complete/fasta_renamed/'
genomes = os.listdir(genome_folder)
all_genomes = []

for gen in genomes:
    if not os.path.exists(new_folder+gen.replace('.fna', '.tax.fna')):
        all_genomes.append(gen)
        
print(len(genomes), len(all_genomes))
genomes = all_genomes

def rename_genome(gen):
    with open('ncbi_taxid.dict', 'rb') as f:
        acc_dict = pickle.load(f)
    try:
        gen_name = gen.replace('_genomic.fna.gz', '')
        tn = gen.replace('.fna', '.tax.fna')
        if os.path.exists(new_folder+tn):
            with open('logfile_ncbi_refseq.txt', 'w+') as f:
                f.write("Already have "+gen+" saved in the new folder")
            return
        if gen_name not in acc_dict:
            with open('logfile_ncbi_refseq.txt', 'w+') as f:
                f.write("Couldn't gen the taxid for "+gen)
            return
        print(gen_name, tn)
        taxid = acc_dict[gen_name]
        new_records = []
        with gzip.open(genome_folder+gen, "rt") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                newid = record.id+"|kraken:taxid|"+str(taxid)
                newseq = SeqRecord(record.seq, id=newid, description=record.description)
                new_records.append(newseq)
        with gzip.open(new_folder+tn, "wt") as file:
            SeqIO.write(new_records, file, "fasta")
    except:
        with open('logfile_ncbi_refseq.txt', 'w+') as f:
            f.write("Couldn't rename "+gen)
    return

def run_multiprocessing(func, i, n_processors):
    with Pool(processes=n_processors) as pool:
        return pool.map(func, i)

def main():

    n_processors =12
    print('Starting processing')
    run_multiprocessing(rename_genome, genomes, n_processors)


if __name__ == "__main__":
    freeze_support()   # required to use multiprocessing
    main()