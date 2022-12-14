import os
import pandas as pd

dfs = []
results_dir = 'kraken_test_results/'
results = os.listdir(results_dir)

all_names = [f.replace('.kreport', '').replace('.bracken', '').replace('_bracken_species', '') for f in results]
all_names = list(set(all_names))
all_names = [n for n in all_names if '.kraken' not in n and '_bracken' not in n]
results = [f for f in results if '.bracken' in f]

for fn in results:
  f = pd.read_csv(results_dir+fn, index_col=0, header=0, sep='\t')
  f = f.set_index('taxonomy_id').loc[:, ['new_est_reads']].rename(columns={'new_est_reads':fn.replace('.bracken', '')})
  if f.shape[0] == 0: continue
  dfs.append(f)

bracken_out = pd.concat(dfs).fillna(value=0)
bracken_out = bracken_out.groupby(by=bracken_out.index, axis=0).sum()
for f in all_names:
  if f not in bracken_out.columns:
    bracken_out[f] = 0

bracken_out.to_csv('bracken_out.csv')
