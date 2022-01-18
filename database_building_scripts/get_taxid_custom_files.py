#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 12 17:54:28 2021

@author: robynwright
"""
import pandas as pd
import numpy as np
import argparse
import sys

parser = argparse.ArgumentParser(description='This requires you to input a list of genomes and their taxonomy in a tab-delimited file and will output the nodes.dmp, names.dmp and nucl.taxid2accession files\n \
                                 The taxonomy should be a string separated by colons, i.e. d__;p__;c__;o__;f__;g__;s__ and the file should not have a header')
parser.add_argument('--genome_taxonomy', dest='genomes', default=None,
                    help='This should be a tab delimited .tsv file with two columns - one containing a list of genomes \n\
                        and the second with their taxonomy (as is output by download_GTDB_latest.py - genomes_and_taxonomy.tsv)\n\
                        The taxonomy should be a string separated by colons, i.e. d__;p__;c__;o__;f__;g__;s__ \n\
                        although no check is run to ensure that this is the case\n\
                        Note that this file should not have a header')
                                
args = parser.parse_args()
genomes = args.genomes
if genomes == None:
    print("You haven't given any input file.\nRun python get_taxid_custom_files.py --help to see how this should be run")
    sys.exit()
genomes = pd.read_csv(genomes, index_col=None, header=None, sep='\t').rename(columns={0:'accession', 1:'gtdb_taxonomy'})
print('Read in the file fine')

print('Starting to make the taxonomy structure')
root, superkingdom, kingdom, phylum, classs, order, family, genus, species, strain = {'root':[1,1]}, {}, {}, {}, {}, {}, {}, {}, {}, {}
root_children, superkingdom_children, kingdom_children, phylum_children, classs_children, order_children, family_children, genus_children, species_children = {'root':[]}, {}, {}, {}, {}, {}, {}, {}, {}

genomes = genomes.set_index('accession')
# for acc in genomes.index.values:
#   tax = genomes.loc[acc, 'gtdb_taxonomy']
#   genomes.loc[acc, 'gtdb_taxonomy'] = tax

genomes['taxid'] = ""
genomes['parent'] = ""

start = ['d__', 'k__', 'p__', 'c__', 'o__', 'f__', 'g__', 's__', 's1__']
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
  taxid, parent = species[sp_name] #get the taxid and parent taxid for this species
  genomes.loc[gen, 'taxid'] = taxid #add the taxid to the right column of the dataframe
  genomes.loc[gen, 'parent'] = parent #add the parent taxid to the right column of the dataframe

genomes.to_csv('db_samples.tsv', sep='\t') #save the file that fits the format of the samples.txt file (we'll be able to use this to assign the taxid's to the fasta file headers)
print('Saved the genomes with taxonomy, taxid and parent taxid as db_samples.tsv. This file is needed for the rename-fasta_headers.py script\n')

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
taxid2accession.to_csv('nucl.taxid2accession', sep='\t')
print('Saved the nucl.taxid2accession file\n')