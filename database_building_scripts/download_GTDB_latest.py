#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 12 15:22:33 2021

@author: robynwright
"""

import pandas as pd
import numpy as np
import os
import argparse
import sys
import subprocess

parser = argparse.ArgumentParser(description='This script will download the latest GTDB representative genome files and build a taxonomy structure that Kraken2 can use.\n \
                                         The only option that you can provide here is to not download the genome files (if all you want is to build the taxonomy from the latest release)\n \
                                         You can do this by running python download_GTDB_latest.py --download_genomes False\n \
                                         Otherwise, just run python download_GTDB_latest.py\n \
                                         If a file called genomes.tsv exists in the directory that you run this script from, it will try to add the accession and taxonomy for these\
                                         This should be a tab delimited file with the first column being the genome accession name and the second being the taxonomy in the format d__;k__;p__;c__;o__;f__;g__;s__;s1__\
                                         i.e. in domain/kingdom/phylum/class/order/family/genus/species/strain format and including the prefixes for each level with two underscores.')
parser.add_argument('--download_genomes', dest='download_genomes', default=True,
                    help="False if you don't want to download the GTDB genomes. Default is True")
args = parser.parse_args()
download_genomes = args.download_genomes

try:
    print('Getting the GTDB metadata files')
    if os.path.exists('ar122_metadata.tsv'):
        print('Already had the archaea metadata file ar122_metadata.tsv. Move this to a different folder if you want to download this again')
    else:
        os.system('wget -q https://data.gtdb.ecogenomic.org/releases/latest/ar122_metadata.tar.gz')
        os.system('tar -xf ar122_metadata.tar.gz')
        files = [f for f in os.listdir(os.getcwd()) if 'ar122_metadata' in f and 'tsv' in f][0]
        os.system('mv '+files+' ar122_metadata.tsv')
        print('Successfully downloaded and unzipped ar122_metadata.tar.gz and its now saved as ar122_metadata.tsv')
    if os.path.exists('bac120_metadata.tsv'):
        print('Already had the bacteria metadata file bac120_metadata.tsv. Move this to a different folder if you want to download this again\n')
    else:
        os.system('wget -q https://data.gtdb.ecogenomic.org/releases/latest/bac120_metadata.tar.gz -q')
        os.system('tar -xf bac120_metadata.tar.gz')
        files = [f for f in os.listdir(os.getcwd()) if 'bac120_metadata' in f and 'tsv' in f][0]
        os.system('mv '+files+' bac120_metadata.tsv')
        print('Successfully downloaded and unzipped bac120_metadata.tar.gz and its now saved as bac120_metadata.tsv\n')   
except:
    print("Failed to download GTDB metadata files. These should be saved as ar122_metadata.tar.gz and bac120_metadata.tar.gz. If these aren't downloaded, you can check whether the file paths are dead and modify these in the script. They should be here in on the GTDB website: https://data.gtdb.ecogenomic.org/releases/latest/\n")
    sys.exit()

print('Opening the metadata files')
bac, arc = pd.read_csv('bac120_metadata.tsv', header=0, sep='\t', low_memory=False), pd.read_csv('ar122_metadata.tsv', header=0, sep='\t', low_memory=False)
bac = bac.loc[:, ['gtdb_genome_representative', 'gtdb_taxonomy']]
arc = arc.loc[:, ['gtdb_genome_representative', 'gtdb_taxonomy']]
genomes = pd.concat([bac, arc]).set_index('gtdb_genome_representative')
genomes = genomes.groupby(by=genomes.index, axis=0).first()
new_acc_dict = {}
for acc in genomes.index.values: 
    if acc[2] == '_': new_acc_dict[acc] = acc[3:]
start = ['d__', 'k__', 'p__', 'c__', 'o__', 'f__', 'g__', 's__', 's1__']
genomes = genomes.rename(index=new_acc_dict).reset_index().rename(columns={'gtdb_genome_representative':'accession'})
for acc in genomes.index.values:
    tax = genomes.loc[acc, 'gtdb_taxonomy']
    tax = tax.split(';')
    this_dict = {}
    for t in range(len(tax)):
        this_dict[tax[t][:3]] = tax[t][3:]
    tax_new = ''
    for s in range(len(start)):
        if tax_new != '': tax_new += ';'
        if start[s] in this_dict:
            tax_new += start[s]+this_dict[start[s]]
        else:
            if s == 0: tax_new += start[s]+this_dict[start[s+1]]
            else: 
                tax_new += start[s]+this_dict[start[s-1]]
    genomes.loc[acc, 'gtdb_taxonomy'] = tax_new
                
genomes.to_csv('genomes_and_taxonomy.tsv', sep='\t', index=False, header=False)
print('Combined the metadata files to keep only the accession and taxonomy columns. This is saved as genomes_and_taxonomy.csv if you would like to look. It has '+str(genomes.shape[0])+' genomes\n')

if os.path.exists('genomes.tsv'):
    other_genomes = pd.read_csv('genomes.tsv', index_col=None, header=None, sep='\t').rename(columns={0:'accession', 1:'gtdb_taxonomy'})
    genomes = pd.concat([genomes, other_genomes])
    genomes.to_csv('genomes_and_taxonomy_added_custom.tsv', sep='\t', index=False, header=False)
    print('Your genomes have been added to the list of GTDB genomes. If you want to see the complete list then this is in genomes_and_taxonomy_added_custom.tsv')

print('Starting to make the taxonomy structure')
root, superkingdom, kingdom, phylum, classs, order, family, genus, species, strain = {'root':[1,1]}, {}, {}, {}, {}, {}, {}, {}, {}, {}
root_children, superkingdom_children, kingdom_children, phylum_children, classs_children, order_children, family_children, genus_children, species_children = {'root':[]}, {}, {}, {}, {}, {}, {}, {}, {}
genomes = genomes.set_index('accession')

genomes['taxid'] = ""
genomes['parent'] = ""

dict_list = [superkingdom, kingdom, phylum, classs, order, family, genus, species, strain]
children_list = [superkingdom_children, kingdom_children, phylum_children, classs_children, order_children, family_children, genus_children, species_children]
taxa = set(genomes.loc[:, 'gtdb_taxonomy'].values)
count = 2
rank_names = ['superkingdom', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'strain']
root_name1 = [1, '|', 'all', '|', '', '|', 'synonym', '|']
root_name2 = [1, '|', 'root', '|', '', '|', 'scientific name', '|']
root_node = [1, '|', 1, '|', 'no rank', '|', 'XX', '|', 0, '|', 0, '|', 11, '|', 1, '|', 1, '|', 0, '|', 0, '|', 0, '|']
names, nodes = [root_name1, root_name2], [root_node]
for tax in sorted(taxa):
  tax = tax.split(';')
  for a in range(len(tax)):
    if start[a] in tax[a]: #check that what we're looking at actually belongs to the level that we want to look at
      if tax[a] not in dict_list[a]: #if we haven't already added this level
        taxid = count #make taxid be the count
        if a > 0: #if we're looking at above the superkingdom level
          parent = dict_list[a-1][tax[a-1]][0] #get the taxid of the parent
        else:
          parent = 1 #otherwise, if we're looking at superkingdom level, the parent taxid must be 1
        dict_list[a][tax[a]] = [taxid, parent] #add the taxid and the parent taxid to our dictionary list
        if a < 8: #if we're not at the strain level
          children_list[a][tax[a]] = [] #add this taxa to the list, so we can add its children later
        if a > 0:  #if we're not at the superkingdom level
          children_list[a-1][tax[a-1]] = children_list[a-1][tax[a-1]]+[tax[a]] #add this child to its parent list
        else: #if we are at the superkingdom level
          root_children['root'] = root_children['root']+[tax[a]] #add this child to its parent list
        #node = [tax_id, parent tax_id, rank (phylum, class, ....), embl code=XX, division id=0, inherited div flag=0, genetic code id=11, inherited GC  flag=1, mitochondrial genetic code id=1, inherited MGC flag=0, GenBank hidden flag=0, hidden subtree root flag=0]
        node = [taxid, '|', parent, '|', rank_names[a], '|', 'XX', '|', 0, '|', 0, '|', 11, '|', 1, '|', 1, '|', 0, '|', 0, '|', 0, '|'] #make a list that fits the format of the nodes.dmp file
        #name = [tax_id, name_txt, unique name, name class]
        name = [taxid, '|', tax[a], '|','', '|', 'scientific name', '|'] #make a list that fits the format of the names.dmp file
        names.append(name) #add the name to the list
        nodes.append(node) #add the node to the list
        count += 1 #add to the count

for gen in genomes.index.values: #for each accession
  sp_name = genomes.loc[gen, 'gtdb_taxonomy'].split(';')[-1] #the species name is at the end of the GTDB taxonomy
  taxid, parent = strain[sp_name] #get the taxid and parent taxid for this species
  genomes.loc[gen, 'taxid'] = taxid #add the taxid to the right column of the dataframe
  genomes.loc[gen, 'parent'] = parent #add the parent taxid to the right column of the dataframe

genomes.to_csv('db_samples.tsv', sep='\t') #save the file that fits the format of the samples.txt file (we'll be able to use this to assign the taxid's to the fasta file headers)
print('Saved the genomes with taxonomy, taxid and parent taxid as db_samples.tsv\n')

#getting a list of all taxid's sorted by rank
superkingdom_children, kingodm_children, phylum_children, classs_children, order_children, family_children, genus_children, species_children = children_list
#taxid_list = ['taxID', 'parent', 'name', 'rank']
taxid_list = [[1, 1, 'root', 'no rank']]

#looping through each child within each parent to get these in order and adding them to a list
for a in root_children['root']:
  taxid_list.append(superkingdom[a]+[a, rank_names[0]])
  for b in superkingdom_children[a]:
    taxid_list.append(kingdom[b]+[b, rank_names[1]])
    for c in kingdom_children[b]:
      taxid_list.append(phylum[c]+[c, rank_names[3]])
      for d in phylum_children[c]:
        taxid_list.append(classs[d]+[d, rank_names[3]])
        for e in classs_children[d]:
          taxid_list.append(order[e]+[e, rank_names[4]])
          for f in order_children[e]:
            taxid_list.append(family[f]+[f, rank_names[5]])
            for g in family_children[f]:
              taxid_list.append(genus[g]+[g, rank_names[6]])
              for h in genus_children[g]:
                taxid_list.append(species[h]+[h, rank_names[7]])
                for i in species_children[h]:
                  taxid_list.append(strain[i]+[i, rank_names[8]])
              
taxid_df = pd.DataFrame(taxid_list, columns=['taxID', 'parent', 'name', 'rank']) #turn the list into a dataframe
taxid_df.to_csv('db_taxid_info.tsv', sep='\t', index=False) #save it to a file
print('Saved all taxid info (taxid, parent taxid, name and rank) as db_taxid_info.tsv\n')

names_df, nodes_df = pd.DataFrame(names), pd.DataFrame(nodes) #turn the names and nodes lists into a dataframe

names_df.to_csv('names.dmp', sep='\t', index=False, header=False) #save the file
nodes_df.to_csv('nodes.dmp', sep='\t', index=False, header=False) #save the file
print('Saved the names.dmp and nodes.dmp files\n')

#now make the taxid2accession files
#We need the accession, accession.version, taxid and gi (filled with NA)

taxid2accession = pd.DataFrame(genomes.loc[:, 'taxid'])
taxid2accession['accession.version'] = list(taxid2accession.index.values)
taxid2accession['gi'] = np.nan
taxid2accession = taxid2accession.loc[:, ['accession.version', 'taxid', 'gi']]

acc_rename = {}
for acc in taxid2accession.index.values:
    if '.' in acc:
        acc_rename[acc] = acc.split('.')[0]
  
taxid2accession = taxid2accession.rename(index=acc_rename)
taxid2accession.to_csv('nucl.accession2taxid', sep='\t')
print('Saved the nucl.accession2taxid file\n')

try:
    if download_genomes != True:
        print('Not downloading genomes, change the option you put in if you meant to download them\
              Either run python download_GTDB_latest.py or python download_GTDB_latest.py --download_genomes True if you want to download them\n')
    else:
        if os.path.exists('gtdb_genomes_reps/'):
            print('You already have a folder called gtdb_genomes_reps. Move this to a different folder if you want to download this again\n')
        else:
            print('Getting the representative GTDB genomes')
            os.system('wget -q https://data.gtdb.ecogenomic.org/releases/latest/genomic_files_reps/gtdb_genomes_reps.tar.gz')
            print('Downloaded the representative genomes')
            os.system('tar -xf gtdb_genomes_reps.tar.gz')
            files = [f for f in os.listdir(os.getcwd()) if 'gtdb_genomes_reps' in f and 'tar.gz' not in f and f != 'gtdb_genomes_reps']
            if len(files) == 0: 
                print('Something went wrong and with the unzipping. Try running this script again\n')
                sys.exit()
            files = files[0]
            output = str(subprocess.check_output("find "+files+" -type f -name '*.fna.gz'", shell=True)).replace("\\n", " ").split(" ")+str(subprocess.check_output("find "+files+" -type f -name '*.fna'", shell=True)).replace("\\n", " ").split(" ")
            os.system('mkdir gtdb_genomes_reps/')
            for out in output:
                print(out)
                os.system("mv "+out+" gtdb_genomes_reps/"+out.split('/')[-1])
            os.system("rm -r "+files)
            print('Unzipped the genomes folder and moved the fasta files from the subfolders that they were in to gtdb_genomes_reps\n')
except:
    print('Failed to download the representative GTDB genomes\n')
    sys.exit()
    

#Clean up the tar files from the directory
files = [f for f in os.listdir(os.getcwd()) if 'tar.gz' in f]
for f in files:
    os.system('rm '+f)
print('Removed the .tar.gz files that are no longer needed\n')
print('The script has now finished running. ')
