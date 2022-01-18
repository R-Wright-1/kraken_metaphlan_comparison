#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 21 12:06:32 2021

@author: robynwright
"""

import os
from multiprocessing import Pool
import pickle
from multiprocessing import freeze_support
import pandas as pd

# domains = ['archaea', 'bacteria', 'fungi', 'invertebrate', 'plant', 'protozoa', 'vertebrate_mammalian', 'vertebrate_other', 'viral']
# accession = []

# for dom in domains:
#     dom = pd.read_csv('assembly_lists/'+dom+'_assembly_summary.txt', index_col=0, header=1, sep='\t')
#     dom = dom.loc[:, ['taxid', 'ftp_path']]
#     accession.append(dom)
    
# accession = pd.concat(accession)

# acc_dict = {}
# for acc in accession.index.values:
#     acc_dict[accession.loc[acc, 'ftp_path'].split('/')[-1]] = accession.loc[acc, 'taxid']

# with open('ncbi_taxid.dict', 'wb') as f:
#     pickle.dump(acc_dict, f)

with open('ncbi_taxid.dict', 'rb') as f:
    acc_dict = pickle.load(f)

print(len(acc_dict))

choco = pd.read_csv('taxid_metaphlan.csv', index_col=0, header=0)
ids_adding = set(list(choco.index.values))
print(len(ids_adding))
ids_adding = [str(id) for id in ids_adding]

folder = 'fasta_renamed_refseq_complete/'
database = 'kraken_chocophlan/'
genomes = os.listdir(folder)
genomes_adding = []
no_taxid = []

print(len(genomes))

for gen in genomes:
    gen_name = gen.replace('_genomic.tax.fna.gz', '')
    if gen_name in acc_dict:
        taxid = acc_dict[gen_name]
        if str(taxid) in ids_adding:
            genomes_adding.append(gen)
    else:
        if not '_' in gen:
            no_taxid.append(gen)
            continue
        gen_name2 = gen.split('_')
        gen_name2 = gen_name2[0]+'_'+gen_name2[1]
        if gen_name2 in acc_dict:
            taxid = acc_dict[gen_name2]
            if taxid in ids_adding:
                genomes_adding.append(gen)
        else:
            no_taxid.append(gen)

print(len(genomes_adding))

genomes = genomes_adding

with open(database+'genomes_added.txt', 'w') as f:
    for gen in genomes:
        f.write(gen+'\n') 

def unzip_add(gen):
    os.system('gunzip '+folder+gen)
    os.system('kraken2-build --add-to-library '+folder+gen.replace('.gz', '')+' --db '+database)
    os.system('gzip '+folder+gen.replace('.gz', ''))
    print(gen)
    return

def run_multiprocessing(func, i, n_processors):
    with Pool(processes=n_processors) as pool:
        return pool.map(func, i)

def main():

    n_processors =12
    print('Starting processing')
    run_multiprocessing(unzip_add, genomes, n_processors)


if __name__ == "__main__":
    freeze_support()   # required to use multiprocessing
    main()