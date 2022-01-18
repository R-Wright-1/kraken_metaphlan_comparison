#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 14 12:03:11 2021

@author: robynwright

"""

import os
import pandas as pd
import argparse
import sys
from multiprocessing import Pool
from multiprocessing import freeze_support

parser = argparse.ArgumentParser(description='This script is to download all sequences from a certain domain that are present in NCBI refseq (ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/).\n\
                                 This requires the additional packages pandas and Biopython\nAs long as the script is able to download the assembly summary file, then it will create a log_file that tells you about whether each sequence was downloaded or not\n\
                                Re-running it with additional domains will by default just add these to what you already have')
parser.add_argument('--domain', dest='domain', default='all',
                    help='pick which domain to download, deafult is all. Otherwise just give the domains that you want as a list separated by commas: fungi,protozoa,viral,vertebrate_mammalian,invertebrate,vertebrate_other,plant')
parser.add_argument('--complete', dest='complete', default=False, 
                    help="choose whether to only download complete genomes, or all genomes. Default is False, meaning all genomes are downloaded")
parser.add_argument('--folder', dest='folder', default='ncbi_genomes', 
                    help="name of the folder to download the genomes to. If this already exists, the genomes will be added to it. By default this is ncbi_genomes")
parser.add_argument('--download_genomes', dest='download_genomes', default=True,
                    help="False if you don't want to download the genomes. Default is True")
parser.add_argument('--log_file', dest='log_file', default='logfile_download_genomes.txt',
                    help="File to write the log to")
parser.add_argument('--processors', dest='proc', default=1,
                    help="Number of processors to use to rename genome files")

args = parser.parse_args()

complete = args.complete
domains = args.domain
folder = args.folder
download_genomes = args.download_genomes
log_file = args.log_file
n_processors =args.proc

if ',' in domains: domains = domains.split(',')
elif 'all' in domains: domains = ['bacteria','archaea','fungi','protozoa','viral','vertebrate_mammalian','invertebrate','vertebrate_other','plant']
else: domains = [domains]

wd = os.getcwd()+"/"
assembly_summaries = []

for domain in domains:
    try:
        if not os.path.exists(str(domain)+"_assembly_summary.txt"):
            os.system("wget -q ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/"+str(domain)+"/assembly_summary.txt")
            os.system("mv assembly_summary.txt "+str(domain)+"_assembly_summary.txt")
        else:
            print("Already got "+str(domain)+"_assembly_summary.txt")
        summary = pd.read_csv(str(domain)+"_assembly_summary.txt", sep='\t', header=1, index_col=0)
        summary = summary.loc[:, ['taxid', 'species_taxid', 'organism_name', 'assembly_level', 'ftp_path']]
        summary['Domain'] = str(domain)
        assembly_summaries.append(summary)
    except:
        print("Unable to download assembly_summary.txt for "+domain)
        

assembly_summaries = pd.concat(assembly_summaries)
print('Joined the assembly summaries\n')

try:
    if not os.path.exists('rankedlineage.dmp'):
        os.system('wget -q https://ftp.ncbi.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz')
        os.system('tar -xf new_taxdump.tar.gz')
    full_lineage = pd.read_csv('rankedlineage.dmp', sep='|', header=None, index_col=0)
    print('Got the full lineage from the current NCBI taxdump\n')
    other_files = ["typeoftype.dmp", "typematerial.dmp", "taxidlineage.dmp", "nodes.dmp", "names.dmp", "merged.dmp",
"host.dmp", "gencode.dmp", "division.dmp", "fullnamelineage.dmp", "delnodes.dmp", "citations.dmp", 'new_taxdump.tar.gz']
    for f in other_files:
        if os.path.exists(f):
            os.system('rm '+f)
    print('Removed the other files from the taxdump folder')
except:
    print("Couldn't get the full lineage from the current NCBI taxdump\n")
    sys.exit()

assembly_summaries['Taxonomy'] = ''

prefix = ['d__', 'k__', 'p__', 'c__', 'o__', 'f__', 'g__', 's__', 's1__']
wd = os.getcwd()
if not os.path.exists(folder):
    os.system('mkdir '+folder)

for acc in assembly_summaries.index.values:
    taxid = assembly_summaries.loc[acc, 'taxid']
    if taxid in full_lineage.index.values:
        taxonomy = list(full_lineage.loc[taxid, :].values)
        tax_string = ''
        taxonomy.reverse()
        count = 0
        new_tax = []
        for tax in taxonomy:
            if isinstance(tax, str):
                new_tax.append(tax.replace('\t', ''))
        taxonomy = new_tax
        for t in range(len(taxonomy)):
            if t != 0: tax_string += ';'
            if taxonomy[t] == '':
                if t == 1 and assembly_summaries.loc[acc, 'Domain'] == 'protozoa': 
                    taxonomy[t] = 'Protista'
                elif t < 7:
                    taxonomy[t] = previous
                else:
                    taxonomy[t] = assembly_summaries.loc[acc, 'organism_name']
            tax_string += prefix[t]+taxonomy[t]
            previous = taxonomy[t]
        assembly_summaries.loc[acc, 'Taxonomy'] = tax_string

assembly_summaries = assembly_summaries.rename(columns={'taxid':'NCBI_taxid', 'species_taxid':'NCBI_species_taxid'})
assembly_summaries.to_csv('summary_to_download.csv')
print('Got phylogeny and file paths for all domains to include. This is saved as summary_to_download.csv\n')

def write_log(message):
    with open(log_file, 'w+') as f:
        f.write(message+'\n')

taxonomy = [[], []]
os.chdir(folder)
def download_genomes(acc):
    if complete:
            if assembly_summaries.loc['assembly_level'] != 'Complete Genome':
                write_log("Didn't get "+acc+" because it wasn't complete")
                return
    try:
            ftp_path = assembly_summaries.loc[acc, 'ftp_path']
            aname = ftp_path.split('/')[-1]+'_genomic.fna.gz'
            if not os.path.exists(aname):
                print(aname)
                ftp_path = ftp_path+'/'+aname
                cmd = 'wget -q '+ftp_path
                os.system(cmd)
                if os.path.exists(aname):
                    taxonomy[0].append(acc), taxonomy[1].append(assembly_summaries.loc[acc, 'Taxonomy'])
                else:
                    write_log("Couldn't get "+acc)
            else:
                print(aname)
                taxonomy[0].append(acc), taxonomy[1].append(assembly_summaries.loc[acc, 'Taxonomy'])
                write_log("Already had "+acc+" as file "+aname+" so didn't download it again")
    except:
            write_log("Couldn't get "+acc)
    return

if not download_genomes:
    print("Finished running everything. The genomes haven't been downloaded because you didn't want to download them.")
    sys.exit()
    
def run_multiprocessing(func, i, n_processors):
    with Pool(processes=n_processors) as pool:
        return pool.map(func, i)

def main():

    print('Starting processing')
    run_multiprocessing(download_genomes, assembly_summaries.index.values, int(n_processors))


if __name__ == "__main__":
    freeze_support()   # required to use multiprocessing
    main()



taxonomy = pd.DataFrame(taxonomy).transpose()
taxonomy.to_csv('genomes.tsv', sep='\t', index=False, header=False)



print("Finished running everything. The genomes should be downloaded in "+folder+" and the list of these and their taxonomy is in genomes.tsv\
      \nA log of any genomes that couldn't be downloaded is in logfile.txt")