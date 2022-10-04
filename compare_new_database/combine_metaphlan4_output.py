import os
import pandas as pd

results_dir = 'metaphlan_out/'
results = os.listdir(results_dir)

results = [f for f in results]
dfs = []

for fn in results:
  f = pd.read_csv(results_dir+fn, index_col=0, header=5, sep='\t')
  keeping = []
  rename_row = {}
  for row in f.index:
    if 's__' in row and 't__' not in row:
       keeping.append(row)
       rename_row[f.loc[row, 'clade_taxid']] = f.loc[row, 'clade_taxid'].split('|')[-1]
       if len(f.loc[row, 'clade_taxid'].split('|')[-1]) < 1:
         rename_row[f.loc[row, 'clade_taxid']] = f.loc[row, 'clade_taxid'].split('|')[-2]
  f = f.loc[keeping, :].set_index('clade_taxid')
  f = f.rename(index=rename_row)
  f = f.loc[:, ['estimated_number_of_reads_from_the_clade']].rename(columns={'estimated_number_of_reads_from_the_clade':fn.replace('.txt', '')})
  dfs.append(f)

metaphlan_out = pd.concat(dfs).fillna(value=0)
metaphlan_out = metaphlan_out.astype(float)
metaphlan_out = metaphlan_out.groupby(by=metaphlan_out.index, axis=0).sum()

metaphlan_out.to_csv('metaphlan_out.csv')
