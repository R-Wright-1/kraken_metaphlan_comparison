import os
import pandas as pd

dfs = []
results_dir = 'test_results/'
results = os.listdir(results_dir)
for fn in results:
  f = pd.read_csv(results_dir+fn, index_col=0, header=0, sep='\t')
  f = f.set_index('taxonomy_id').loc[:, ['new_est_reads']].rename(columns={'new_est_reads':fn.replace('.bracken', '')})
  if f.shape[0] == 0: continue
  dfs.append(f)

bracken_out = pd.concat(dfs).fillna(value=0)
bracken_out = bracken_out.groupby(by=bracken_out.index, axis=0).sum()
bracken_out.to_csv('bracken_out.csv')
